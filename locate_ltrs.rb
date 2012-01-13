# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=73950699&retmode=xml&rettype=text
# bio = Bio::FastaFormat.new(File.read("/Users/evansenter/Downloads/sequence.fasta"))

require "bio"
require "entrez"
require "nokogiri"

class BlastParser
  FLANKING_DISTANCE = 5_000
  ERV_DISTANCE      = 4_000..10_000
  
  attr_reader :xml, :query_length

  def self.bootstrap(path)
    new(Nokogiri::XML(File.open(path)))
  end
  
  def initialize(xml)
    @xml          = xml
    @query_length = xml.xpath("BlastOutput/BlastOutput_query-len").inner_text.to_i
  end
  
  def parse_xml
    xml.xpath("//BlastOutput_iterations//Iteration//Iteration_hits//Hit").map do |hit_group|
      {
        id:        hit_group.xpath("Hit_id").inner_text.split("|")[1],
        accession: hit_group.xpath("Hit_accession").inner_text,
        locations: parse_coordinates(hit_group)
      }
    end
  end
  
  def filter
    parse_xml.map do |hash|
      hash.merge(ltr_pairs: filter_array(hash[:locations][:plus]) + filter_array(hash[:locations][:minus]))
    end.reject do |hash|
      hash[:ltr_pairs].empty?
    end
  end
  
  def flanking_regions!
    filter.map do |hash|
      hash.merge(flanking_pairs: flanking_pairs(hash[:id], hash[:ltr_pairs]))
    end
  end
  
  def flanking_pairs(id, pairs_array)
    pairs_array.map do |putative_erv|
      {
        coordinates: putative_erv,
        flanking_5:  sequence_from_entrez(id, putative_erv.first.begin - FLANKING_DISTANCE, putative_erv.first.begin - 1),
        flanking_3:  sequence_from_entrez(id, putative_erv.last.end + 1, putative_erv.last.end + FLANKING_DISTANCE)
      }
    end
  end
  
  def sequence_from_entrez(id, start, stop)
    fasta = Entrez.EFetch("nuccore", {
      id:        id, 
      seq_start: start, 
      seq_stop:  stop, 
      retmode:   :fasta, 
      rettype:   :text
    }).response.body
    
    Bio::FastaFormat.new(fasta).seq
  end
  
  private
  
  def filter_array(coordinate_array)
    coordinate_array.combination(2).select do |lower, upper|
      ERV_DISTANCE.include?(upper.begin - lower.end)
    end
  end
  
  def parse_coordinates(hit_group)
    partitioned_coordinates = hit_group.xpath("Hit_hsps/Hsp").map do |hit|
      [hit.xpath("Hsp_hit-from").inner_text.to_i, hit.xpath("Hsp_hit-to").inner_text.to_i]
    end.select do |coordinates|
      (coordinates.last - coordinates.first).abs > query_length * 0.9
    end.partition do |coordinates|
      coordinates.last > coordinates.first
    end.map do |array|
      array.map { |coordinates| Range.new(*coordinates.sort) }.sort_by(&method(:midpoint))
    end
    
    Hash[[:plus, :minus].zip(partitioned_coordinates)]
  end
  
  def midpoint(range)
    (range.begin + range.end) / 2.0
  end
end

parser    = BlastParser.bootstrap("/Users/evansenter/Documents/School/BC/Rotation 3 - Johnson/Canis Lupus Familiaris/GY1D8VT9016-Alignment.xml")
ltr_pairs = parser.filter