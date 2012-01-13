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
    end.reject do |hash|
      hash[:flanking_pairs][:plus].empty? && hash[:flanking_pairs][:minus].empty?
    end
  end
  
  def flanking_regions_as_fasta!(path)
    raise "File can't already exist" if File.exists?(path)
    
    fasta_content = formatted_flanking_regions
    
    File.open(path, "w") do |file|
      fasta_content.each do |fasta_sequence|
        file.write(fasta_sequence[:comment])
        file.write("\n")
        file.write(fasta_sequence[:sequence].scan(/.{1,120}/).join("\n"))
        file.write("\n\n")
      end
    end
  end
  
  private
  
  def formatted_flanking_regions
    flanking_regions!.inject([]) do |fasta_sequences, hash|
      plus_regions = hash[:flanking_pairs][:plus].inject([]) do |array, putative_erv|
        array.concat([
          formatted_flanking_region(hash[:definition], "5' genomic region", "plus", putative_erv[:putative_erv].first.hit_from, -FLANKING_DISTANCE..-1, putative_erv[:flanking_5])
          formatted_flanking_region(hash[:definition], "3' genomic region", "plus", putative_erv[:putative_erv].first.hit_to, 1..FLANKING_DISTANCE, putative_erv[:flanking_3])
        ])
      end
        
      minus_regions = hash[:flanking_pairs][:minus].inject([]) do |array, putative_erv|
        array.concat([
          formatted_flanking_region(hash[:definition], "5' genomic region", "minus", putative_erv[:putative_erv].first.hit_to, 1..FLANKING_DISTANCE, putative_erv[:flanking_5])
          formatted_flanking_region(hash[:definition], "3' genomic region", "minus", putative_erv[:putative_erv].first.hit_from, -FLANKING_DISTANCE..-1, putative_erv[:flanking_3])
        ])
      end
      
      fasta_sequences.concat(plus_regions).concat(minus_regions)
    end
  end
  
  def formatted_flanking_region(definition, region, strand, position, window, sequence)
    {
      comment: "> %s (%s on %s strand: %s - coordinates %s to %s" % [
        definition,
        region,
        strand,
        position + window.begin,
        position + window.end
      ],
      sequence: sequence
    }
  end
  
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
      hsps.sort_by(&method(:midpoint))
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

# parser = BlastParser.bootstrap("/Users/evansenter/Documents/School/BC/Rotation 3 - Johnson/Canis Lupus Familiaris/GY1D8VT9016-Alignment.xml")
# parser.flanking_regions_as_fasta!("/Users/evansenter/Desktop/genomic_regions.fasta")