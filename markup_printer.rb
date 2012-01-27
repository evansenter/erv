class MarkupPrinter
  attr_reader :erv, :tens_of_aa
  
  def initialize(erv, tens_of_aa)
    @erv        = erv
    @tens_of_aa = tens_of_aa
  end
  
  def format
    formatted_indices.zip(formatted_pointers, *formatted_sequences, formatted_annotations).inject("") do |string, lines| 
      string + lines.join("\n") + "\n\n"
    end
  end
  
  def print
    super(format)
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