raise "I DON'T THINK THIS WORKS"

require "bio"
require "./fasta_printer.rb"

class NormalizeGc
  attr_reader :folders, :orf_suffix, :window_size, :flanking_length
  
  def initialize(directory, orf_suffix = "env", window_size = 1_500, flanking_length = 5_000)
    @folders         = Dir[File.join(directory, "*")]
    @orf_suffix      = orf_suffix
    @window_size     = window_size
    @flanking_length = flanking_length
  end
  
  def write_best_windows!
    folders.map(&method(:write_best_window_by_gc_content))
  end
  
  def write_best_window_by_gc_content(folder)
    target_gc_content  = load_from_folder(folder, orf_suffix).gc_content.to_f
    erv_sequence       = load_from_folder(folder, "erv")[-5000..-1]
    best_sequence_hash = find_best_sequence(erv_sequence, target_gc_content)
    
    GcResult.new(target_gc_content, best_sequence_hash).write_fasta!(folder)
  end
  
  def find_best_sequence(flanking_region, target_gc_content)
    matches = (0..(flanking_length - window_size)).map do |i|
      sequence = flanking_region[i, window_size]
      
      {
        index:                 i,
        sequence:              sequence,
        gc_content_difference: (target_gc_content - sequence.gc_content.to_f).abs
      }
    end
      
    matches.min { |a, b| a[:gc_content_difference] <=> b[:gc_content_difference] }
  end
  
  def load_from_folder(folder, suffix)
    file_path = Dir[File.join(folder, "*#{suffix}.txt")].first
    
    Bio::Sequence::NA.new(Bio::FastaFormat.new(File.read(file_path)).seq)
  end
  
  class GcResult
    include FastaPrinter
    
    attr_reader :target_gc_content, :best_sequence_hash
    
    def initialize(target_gc_content, best_sequence_hash)
      @target_gc_content  = target_gc_content
      @best_sequence_hash = best_sequence_hash
    end
    
    def fasta_hash
      {
        comment: "> Deviated from GC content of %.2f by %.2f percent, starts in 3' genomic region at position %d (0-indexed)" % [
          target_gc_content * 100,
          best_sequence_hash[:gc_content_difference] * 100,
          best_sequence_hash[:index]
        ],
        sequence: best_sequence_hash[:sequence]
      }
    end
    
    def fasta_filename
      "Highest matching GC region.fasta"
    end
  end
end

# path = "/Users/evansenter/Documents/School/BC/Rotations/Rotation 3 - Johnson/Kate's Stuff/GC Sequences"