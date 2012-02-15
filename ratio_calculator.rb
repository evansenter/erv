class RatioCalculator
  include Enumerable
  
  TRANSITION_HASH = {
    "a" => "g",
    "g" => "a",
    "c" => "t",
    "t" => "c"
  }
  
  attr_reader :query, :target, :length
  
  def initialize(query, target)
    @query, @target = query.downcase, target.downcase
    
    if query.length != target.length
      raise "Both FASTA objects must have a sequence of the same length"
    else
      @length = query.length
    end
  end
  
  def ratio
    counter.tap do |hash|
      total_count = hash.values.inject(&:+)
      
      hash[:transitions]   /= total_count
      hash[:transversions] /= total_count
    end
  end
  
  def counter
    counter_hash = inject({ transitions: 0.0, transversions: 0.0 }) do |hash, (from, to)|
      hash.tap do
        if from != to
          hash[TRANSITION_HASH[from] == to ? :transitions : :transversions] += 1
        end
      end
    end
    
    if counter_hash[:transitions].zero? || counter_hash[:transversions].zero?
      { transitions: 1.0, transversions: 1.0 }
    else
      counter_hash
    end
  end
  
  def each(&block)
    query.each_char.zip(target.each_char).each(&block)
  end
end