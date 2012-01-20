require "bio"
require "./fasta_parser.rb"
require "./printer.rb"

class Erv
  include FastaParser
  
  attr_reader :sequence, :annotations
  
  class << self
    def bootstrap_from_dirty_file(path)
      file     = File.read(path)
      sequence = file.gsub(/[^atcgn]/, "")
      new(Bio::Sequence::NA.new(sequence))
    end
    
    def bootstrap_from_bio_fasta(fasta)
      new(fasta.naseq)
    end
  end
  
  def initialize(sequence)
    @sequence    = sequence
    @annotations = []
  end
  
  def reading_frames
    @reading_frames ||= (1..3).map { |i| { i => sequence.translate(i) } }.inject(&:merge)
  end
  
  def annotate_orf(substring, name = nil)
    orf(substring).each do |orf_region|
      start_annotation_string = "[#{name || (substring.length > 10 ? substring[0, 7] + '...' : substring)} start (frame #{orf_region[:frame]}) =>"
      end_annotation_string   = "<= #{name || (substring.length > 10 ? substring[0, 7] + '...' : substring)} end (frame #{orf_region[:frame]})]"
      
      annotate(orf_region[:na_position].begin, start_annotation_string)
      annotate(orf_region[:na_position].end - end_annotation_string.length - 1, end_annotation_string)
    end
  end
  
  def orf(substring)
    locate_aa(substring).map do |frame, starting_aa| 
      matching_orf = reading_frames[frame][starting_aa..-1][/^.*?\*/]
      starting_na  = (3 * starting_aa + frame) # Switching from 0-indexed to 1-indexed
      
      { 
        sequence:    matching_orf,
        na_position: starting_na..(starting_na + matching_orf.length * 3 - 1),
        frame:       frame
      }
    end
  end
  
  def locate_aa(substring)
    reading_frames.map { |frame, sequence| { frame => sequence.index(substring.upcase) } }.reject { |hash| hash.values.compact.empty? }.inject(&:merge).to_a
  end
  
  def annotate(position, label)
    annotations.push([position, label]).sort!
  end
  
  def clear_annotations
    @annotations.clear
  end
  
  def printer(tens_of_aa = 6)
    Printer.new(self, tens_of_aa)
  end
  
  def print(tens_of_aa = 6)
    printer(tens_of_aa).print
  end
end

# erv = Erv.bootstrap_from_dirty_file("/Users/evansenter/Documents/School/BC/Rotations/Rotation 3 - Johnson/Canis Lupus Familiaris/Dog Fc Consensus.txt")
# erv.annotate_orf("mgtsqsk", "Gag")
# erv.annotate_orf("rspgsat", "ProPol")
# erv.annotate_orf("cmkgsgt", "Env")
# erv.print(7)