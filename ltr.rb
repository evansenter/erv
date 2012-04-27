require "./commune.rb"
require "./entrez.rb"
require "./fasta_printer.rb"

class Ltr
  include Commune
  include EntrezSequence
  include FastaPrinter
  
  attr_reader :hit, :hsp, :buffer_size
  
  def initialize(hit, hsp, buffer_size = 0)
    @hit         = hit
    @hsp         = hsp
    @buffer_size = buffer_size
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
  
  def raw_seq
    @raw_sequence ||= na_sequence_from_entrez(hit_id, up_coord, 0...length, buffer_size)
    
    minus_strand? ? @raw_sequence.complement : @raw_sequence
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
    accession == other_ltr.accession && coord_window(from).include?(other_ltr.from) && coord_window(to).include?(other_ltr.to)
  end
  
  def inspect
    "<##{self.class}:#{object_id.to_s(8)} #{hsp.inspect}>"
  end
end