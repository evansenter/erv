require "./entrez.rb"
require "./fasta_printer.rb"

class Ltr
  include EntrezSequence
  include FastaPrinter
  
  attr_reader :hit, :hsp
  
  def initialize(hit, hsp)
    @hit = hit
    @hsp = hsp
  end
  
  def plus_strand?
    hsp.query_frame == hsp.hit_frame
  end
  
  def minus_strand?
    !plus_strand?
  end
  
  def definition
    hit.definition
  end
  
  def accession
    hit.accession
  end
  
  def hit_id
    # gi|224514737|ref|NT_009237.18|
    hit.hit_id.split("|")[1]
  end
  
  def seq
    raw_sequence = na_sequence_from_entrez(hit_id, up_coord, 0...length)
    
    minus_strand? ? raw_sequence.complement : raw_sequence
  end
  
  def length
    down_coord - up_coord + 1
  end
  
  def midpoint
    (hsp.hit_from + hsp.hit_to) / 2.0
  end
  
  def from
    hsp.hit_from
  end
  
  def to
    hsp.hit_to
  end
  
  def up_coord
    [hsp.hit_from, hsp.hit_to].min
  end
  
  def down_coord
    [hsp.hit_from, hsp.hit_to].max
  end
  
  def coord_window(coord)
    (coord - 200)..(coord + 200)
  end
  
  def ==(other_ltr)
    hit.accession == other_ltr.hit.accession && coord_window(from).include?(other_ltr.from) && coord_window(to).include?(to)
  end
  
  def inspect
    "<##{self.class}:#{object_id.to_s(8)} #{hsp.inspect}>"
  end
end