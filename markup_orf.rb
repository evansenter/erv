require "bio"

class Erv
  attr_reader :sequence, :annotations
  
  def self.bootstrap(path)
    raise "Needs to handle n in sequence and different cases, should just use Bio::Sequence"
    
    file     = File.read(path)
    sequence = file.gsub(/[^atcg]/, "")
    new(sequence)
  end
  
  def initialize(sequence)
    @sequence    = Bio::Sequence::NA.new(sequence)
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
  
  def print(tens_of_aa = 6)
    Printer.new(self, tens_of_aa).print
  end
  
  class Printer
    attr_reader :erv, :tens_of_aa
    
    def initialize(erv, tens_of_aa)
      @erv        = erv
      @tens_of_aa = tens_of_aa
    end
    
    def print
      formatted_indices.zip(formatted_pointers, *formatted_sequences, formatted_annotations).each do |lines| 
        lines.each(&method(:puts))
        puts
      end
      
      return nil
    end
    
    private
    
    def formatted_indices
      indices = (0...(erv.sequence.length / 30.0).ceil).map { |i| 10 * i + 1 }.map { |i| "#{((i - 1) * 3) + 1}/#{i}"}.map { |string| "%-30s" % string }
      indices.each_slice(tens_of_aa).map(&:join)
    end
    
    def formatted_pointers
      (erv.sequence.length / 30.0).ceil.times.map { ("|" + " " * 29) }.each_slice(tens_of_aa).map(&:join)
    end
    
    def formatted_sequences
      sequences = (1..3).map { |i| " " * (i - 1) + erv.reading_frames[i].gsub(/./) { |match| "#{match}  " } }.unshift(erv.sequence)
      sequences.map { |sequence| sequence.scan(/.{1,#{30 * tens_of_aa}}/) }
    end
    
    def formatted_annotations
      erv.annotations.inject(" " * erv.sequence.length) do |annotation_string, (position, label)|
        annotation_string.tap do
          annotation_string[(position - 1)..(position + label.length)] = label
        end
      end.scan(/.{1,#{30 * tens_of_aa}}/)
    end
  end
end

erv = Erv.bootstrap("/Users/evansenter/Documents/School/BC/Rotation 3 - Johnson/Dog Fc Consensus.txt")
erv.annotate_orf("mgtsqsk", "Gag")
erv.annotate_orf("rspgsat", "ProPol")
erv.annotate_orf("cmkgsgt", "Env")
erv.print