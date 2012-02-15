require "benchmark"
require "./mutaterv.rb"
require "/Users/evansenter/Documents/School/BC/Rotations/Rotation 3 - Johnson/Ted's Stuff/Evolution Files/Hominid ERV-Fc1 env + ancestral.rb"

{ ancestral: [@acseq, @hseq], chimp: [@ccseq, @cseq], human: [@hcseq, @hseq] }.each do |sequence_type, (species, alignment)|
  { noncoding: 8.73, coding: 38.68 }.each do |sub_indel_type, sub_indel_ratio|
    puts sequence_type
    puts sub_indel_type
    puts "%.3f min" % (Benchmark.realtime { @mutations = Mutaterv.new(species, alignment, @aseq, 9e6, sub_indel_ratio).mutated_stops(1_000_000) } / 60)
    puts @mutations
    puts
  end
  
  print "\n" * 3
end