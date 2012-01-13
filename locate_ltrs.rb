require "bio"
require "entrez"

class BlastParser
  FLANKING_DISTANCE = 5_000
  ERV_DISTANCE      = 4_000..10_000
  
  attr_reader :report

  def self.bootstrap(path)
    new(Bio::Blast::Report.new(File.read(path)))
  end
  
  def initialize(report)
    @report = report
  end
  
  def parse_xml
    report.hits.map do |hit|
      {
        definition: hit.definition,
        id:         hit.hit_id.split("|")[1],
        accession:  hit.accession,
        locations:  parse_coordinates(hit.hsps.select { |hsp| hsp.align_len > report.query_len * 0.9 })
      }
    end
  end
  
  def filter
    parse_xml.map do |hash|
      hash.merge({
        ltr_pairs: {
          plus:  filter_array(hash[:locations][:plus]),
          minus: filter_array(hash[:locations][:minus])
        }
      })
    end.reject do |hash|
      hash[:ltr_pairs][:plus].empty? && hash[:ltr_pairs][:minus].empty?
    end
  end
  
  def flanking_regions!
    filter.map do |hash|
      hash.merge(flanking_pairs: flanking_pairs(hash[:id], hash[:ltr_pairs]))
    end
  end
  
  private
  
  def flanking_pairs(id, pairs_hash)
    {
      plus: pairs_hash[:plus].map do |putative_erv|
        {
          putative_erv: putative_erv,
          flanking_5:  sequence_from_entrez(id, putative_erv.first.hit_from, -FLANKING_DISTANCE..-1),
          flanking_3:  sequence_from_entrez(id, putative_erv.last.hit_to, 1..FLANKING_DISTANCE)
        }
      end,
      minus: pairs_hash[:minus].map do |putative_erv|
        {
          putative_erv: putative_erv,
          flanking_5:  Bio::Sequence::NA.new(sequence_from_entrez(id, putative_erv.last.hit_to, 1..FLANKING_DISTANCE)).complement,
          flanking_3:  Bio::Sequence::NA.new(sequence_from_entrez(id, putative_erv.first.hit_from, -FLANKING_DISTANCE..-1)).complement
        }
      end
    }
  end
  
  def sequence_from_entrez(id, position, window)
    fasta = Entrez.EFetch("nuccore", {
      id:        id, 
      seq_start: position + window.begin, 
      seq_stop:  position + window.end, 
      retmode:   :fasta, 
      rettype:   :text
    }).response.body
    
    Bio::FastaFormat.new(fasta).seq
  end
  
  def filter_array(hsps)
    hsps.combination(2).select do |match_5, match_3|
      ERV_DISTANCE.include?(match_3.hit_from - match_5.hit_to)
    end
  end
  
  def parse_coordinates(hsps)
    partitioned_hsps = hsps.partition(&method(:same_strand?)).map do |hsps|
      hsps.sort_by(&method(:midpoint)) #{ |a, b| midpoint(a) <=> midpoint(b) * (same_strand?(a) ? 1 : -1) }
    end
    
    Hash[[:plus, :minus].zip(partitioned_hsps)]
  end
  
  def midpoint(hsp)
    (hsp.hit_from + hsp.hit_to) / 2.0
  end
  
  def same_strand?(hsp)
    hsp.query_frame == hsp.hit_frame
  end
end

parser             = BlastParser.bootstrap("/Users/evansenter/Documents/School/BC/Rotation 3 - Johnson/Canis Lupus Familiaris/GY1D8VT9016-Alignment.xml")
ltr_pairs          = parser.filter
flanking_sequences = parser.flanking_regions!.map { |hash| hash[:flanking_pairs] }