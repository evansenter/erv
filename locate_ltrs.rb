# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=73950699&retmode=xml&rettype=text
# bio = Bio::FastaFormat.new(File.read("/Users/evansenter/Downloads/sequence.fasta"))

require "nokogiri"

class BlastParser
  ERV_DISTANCE = 4_000..10_000
  
  attr_reader :xml, :query_length

  def self.bootstrap(path)
    new(Nokogiri::XML(File.open(path)))
  end
  
  def initialize(xml)
    @xml          = xml
    @query_length = xml.xpath("BlastOutput/BlastOutput_query-len").inner_text.to_i
  end
  
  def parse_xml
    @parsed_xml ||= xml.xpath("//BlastOutput_iterations//Iteration//Iteration_hits//Hit").map do |hit_group|
      {
        id:        hit_group.xpath("Hit_id").inner_text.split("|")[1],
        accession: hit_group.xpath("Hit_accession").inner_text,
        locations: parse_coordinates(hit_group)
      }
    end
  end
  
  def filter
    parse_xml.map do |hash|
      {
        id:        hash[:id],
        accession: hash[:accession],
        ltr_pairs: filter_array(hash[:locations][:plus]) + filter_array(hash[:locations][:minus])
      }
    end.reject do |hash|
      hash[:ltr_pairs].empty?
    end
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