require "bio"
require "./fasta_printer.rb"

class GcRunner
  WINDOW_DOMAIN = -5000..-1
  
  attr_reader :folders
  
  def initialize(directory)
    @folders = Dir[File.join(directory, "*")]
  end
  
  def load_from_folder(folder, suffix)
    file_path = Dir[File.join(folder, "*#{suffix}.txt")].first
    sequence  = Bio::Sequence::NA.new(Bio::FastaFormat.new(File.read(file_path).gsub(/\r/, "\n")).seq)
    
    sequence.empty? ? (raise "Empty sequence in folder #{folder} with suffix #{suffix}") : sequence
  end
  
  def window(folder)
    gc_sequence   = load_from_folder(folder, "env")
    window_domain = load_from_folder(folder, "erv")
    
    GcWindow.new(gc_sequence, window_domain[WINDOW_DOMAIN])
  end
  
  def write_best_windows!
    folders.map do |folder|
      window(folder).write_best_window!(folder)
    end
  end
end

class GcWindow
  WINDOW_SIZE = 1500
  
  attr_reader :gc_sequence, :window_domain
  
  def initialize(gc_sequence, window_domain)
    @gc_sequence, @window_domain = gc_sequence, window_domain
  end
  
  def target_gc_content
    @sequence_content ||= gc_sequence.gc_content.to_f
  end
  
  def best_window
    window_hashes.min { |a, b| a[:gc_difference] <=> b[:gc_difference] }
  end
  
  def window_hashes
    (0..(window_domain.length - WINDOW_SIZE)).map do |i|
      sequence = window_domain[i, WINDOW_SIZE]
      
      {
        index:         i + 1,
        sequence:      sequence,
        gc_difference: (target_gc_content - sequence.gc_content).abs
      }
    end
  end
  
  def write_best_window!(directory)
    GcPrinter.new(target_gc_content, best_window).write_fasta!(directory)
  end
end

class GcPrinter
  include FastaPrinter
  
  attr_reader :target_gc_content, :best_sequence_hash
  
  def initialize(target_gc_content, best_sequence_hash)
    @target_gc_content  = target_gc_content
    @best_sequence_hash = best_sequence_hash
  end
  
  def fasta_hash
    {
      comment: "> Deviated from GC content of %.2f by %.2f percent, starts in 3' genomic region at position %d (1-indexed)" % [
        target_gc_content * 100,
        best_sequence_hash[:gc_difference] * 100,
        best_sequence_hash[:index]
      ],
      sequence: best_sequence_hash[:sequence]
    }
  end
  
  def fasta_filename
    "Highest matching GC region.fasta"
  end
end

# path = "/Users/evansenter/Documents/School/BC/Rotations/Rotation 3 - Johnson/Kate's Stuff/GC Sequences"