# Start modelling indels from http://mbe.oxfordjournals.org/content/26/7/1523.full.pdf+html paper using human / chimp numbers
# Run once using noncoding assumption, if % maintaining ORF is low can no longer maintain this assumption, and have to use coding numbers
# Insertion length roughly follows a geometric distribution ()

require "bio"
require "./ratio_calculator.rb"

class Mutaterv
  TRANSITION_HASH = {
    "a" => "g",
    "g" => "a",
    "c" => "t",
    "t" => "c"
  }
  
  TRANSVERSION_HASH = {
    "a" => %w(c t),
    "g" => %w(c t),
    "c" => %w(a g),
    "t" => %w(a g)
  }
  
  attr_reader :mutation_sequence, :current_sequence, :ancestral_sequence, :years, :sub_indel_ratio, :mutation_rate
  
  def initialize(mutation_sequence, current_sequence, ancestral_sequence, years, sub_indel_ratio = :sub_only, mutation_rate = 2.2e-9)
    @mutation_sequence  = Bio::Sequence::NA.new(mutation_sequence.downcase)
    @current_sequence   = Bio::Sequence::NA.new(current_sequence.downcase)
    @ancestral_sequence = Bio::Sequence::NA.new(ancestral_sequence.downcase)
    @years              = years
    @sub_indel_ratio    = sub_indel_ratio
    @mutation_rate      = mutation_rate
    
    raise "Sequences must be the same length" unless current_sequence.length == ancestral_sequence.length
  end
  
  def mutated_stops(iterations = 1000)
    iterations.times.inject(Hash.new { |hash, key| hash[key] = 0 }) do |hash, _|
      hash.tap do
        hash[count_stops(mutate)] += 1
      end
    end
  end
  
  def count_stops(sequence)
    translated_sequence = sequence.translate
    
    unless (number_of_stops = translated_sequence.count("*")) == 1
      number_of_stops
    else
      translated_sequence[0...15].include?(?M) && translated_sequence[-15..-1].include?(?*) ? -1 : 1
    end
  end
  
  def mutate
    mutation_ordering.inject(Bio::Sequence::NA.new(mutation_sequence)) do |string, action|
      send(:"#{action}_for", string)
    end
  end
  
  def number_of_mutations
    @number_of_mutations ||= (mutation_rate * ancestral_sequence.length * years).round
  end
  
  def number_of_indels
    @number_of_indels ||= (sub_indel_ratio == :sub_only ? 0 : (number_of_mutations / sub_indel_ratio).round)
  end
  
  def mutation_ratio
    @mutation_ratio ||= RatioCalculator.new(current_sequence, ancestral_sequence).ratio
  end
  
  def sub_for(string)
    string.tap do
      mutation_position         = rand(string.length)
      string[mutation_position] = generate_sub(string[mutation_position])
    end
  end
  
  def indel_for(string)
    if rand(2).zero?
      # Insert into sequence
      string.insert(rand(string.length + 1), generate_insert)
    else
      # Delete from sequence
      deletion_site = rand(string.length)
      string[0, deletion_site] + (string[(deletion_site + indel_length)..-1] || "")
    end
  end
  
  def generate_sub(character)
    if %w(a g c t).include?(character)
      rand < mutation_ratio[:transitions] ? TRANSITION_HASH[character] : TRANSVERSION_HASH[character][rand(2)]
    else
      %w(a g c t)[rand(4)]
    end
  end
  
  def generate_insert
    indel_length.times.map { %w(a g c t)[rand(4)] }.join
  end
  
  def indel_length
    # Closed form sum of geometric series through term n: 2 * (1 - 2 ** (-n - 1)) - 1
    # Randomly generating a sum, and solving for n to generate indices in a geometric series with r = 1/2
    (Math::log2(2 / (1 - rand)) - 1).ceil
  end
  
  def mutation_ordering
    (Array.new(number_of_mutations, :sub) + Array.new(number_of_indels, :indel)).shuffle
  end
  
  def count_differences(string_1, string_2)
    string_1.each_char.zip(string_2.each_char).select { |array| array.uniq.length > 1 }.size
  end
end