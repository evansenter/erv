require "./entrez.rb"
require "./fasta_printer.rb"

class PutativeErv
  include EntrezSequence
  include FastaPrinter
  
  attr_reader :ltrs
  
  def initialize(ltr_1, ltr_2)
    @ltrs = [ltr_1, ltr_2].sort_by(&:midpoint)
    
    raise "LTRs must be in the same orientation" unless ltrs.map(&:plus_strand?).uniq.length == 1
    raise "LTRs must have the same Hit ID" unless ltrs.map(&:hit_id).uniq.length == 1
  end
  
  def plus_strand?
    ltrs.all? { |ltr| ltr.plus_strand? }
  end
  
  def minus_strand?
    ltrs.all? { |ltr| ltr.minus_strand? }
  end
  
  def definition
    ltrs.first.definition
  end
  
  def accession
    ltrs.first.accession
  end
  
  def seq
    raw_sequence = na_sequence_from_entrez(ltrs.first.hit_id, ltrs.first.up_coord, 0...length)
    
    minus_strand? ? raw_sequence.complement : raw_sequence
  end
  
  def length
    ltrs.last.down_coord - ltrs.first.up_coord + 1
  end
  
  def coord_window(coord)
    (coord - 200)..(coord + 200)
  end
  
  def from
    plus_strand? ? ltrs.first.up_coord : ltrs.last.down_coord
  end
  
  def to
    minus_strand? ? ltrs.first.up_coord : ltrs.last.down_coord
  end
  
  def ==(other_erv)
    accession == other_erv.accession && coord_window(from).include?(other_erv.from) && coord_window(to).include?(other_erv.to)
  end
end