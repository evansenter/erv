require "bio"
require "entrez"

class BlastParser
  ERV_DISTANCE = 4_000..12_000
  
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
        locations:  parse_coordinates(hit.hsps.select { |hsp| hsp.align_len > report.query_len * 0.8 })
      }
    end
  end
  
  def filter
    parse_xml.map do |hash|
      plus_ervs, plus_solo_ltrs   = filter_array(hash[:locations][:plus])
      minus_ervs, minus_solo_ltrs = filter_array(hash[:locations][:minus])
      
      hash.merge({
        ltr_pairs: {
          plus:  plus_ervs,
          minus: minus_ervs
        },
        solo_ltrs: {
          plus:  plus_solo_ltrs,
          minus: minus_solo_ltrs
        }
      })
    end
  end
  
  def flanking_regions!
    filter.map do |hash|
      hash.merge(flanking_pairs: flanking_pairs(hash[:id], hash[:ltr_pairs]))
    end.reject do |hash|
      hash[:flanking_pairs][:plus].empty? && hash[:flanking_pairs][:minus].empty?
    end
  end
  
  def write_fasta!(path, fasta_content)
    raise "File can't already exist" if File.exists?(path)
    
    File.open(path, "w") do |file|
      fasta_content.each do |fasta_sequence|
        file.write(fasta_sequence[:comment])
        file.write("\n")
        file.write(fasta_sequence[:sequence].scan(/.{1,120}/).join("\n"))
        file.write("\n\n")
      end
    end
  end
  
  def formatted_erv_regions
    filter.inject([]) do |fasta_sequences, hash|
      plus_regions = hash[:ltr_pairs][:plus].inject([]) do |array, putative_erv|
        start = putative_erv.first.hit_from
        stop  = putative_erv.last.hit_to
        range = 0..(stop - start)
        
        array.concat([formatted_region(hash[:definition], "putative ERV", "plus", hash[:accession], start, range, sequence_from_entrez(hash[:id], start, range))])
      end
        
      minus_regions = hash[:ltr_pairs][:minus].inject([]) do |array, putative_erv|
        start = putative_erv.first.hit_to
        stop  = putative_erv.last.hit_from
        range = 0..(stop - start)
        
        array.concat([formatted_region(hash[:definition], "putative ERV", "minus", hash[:accession], start, range, complement(sequence_from_entrez(hash[:id], start, range)))])
      end
      
      fasta_sequences.concat(plus_regions).concat(minus_regions)
    end
  end
  
  def formatted_genomic_regions(distance)
    flanking_regions!.inject([]) do |fasta_sequences, hash|
      plus_regions = hash[:flanking_pairs][:plus].inject([]) do |array, putative_erv|
        array.concat([
          formatted_region(hash[:definition], "5' genomic region", "plus", hash[:accession], putative_erv[:putative_erv].first.hit_from, -distance..-1, putative_erv[:flanking_5]),
          formatted_region(hash[:definition], "3' genomic region", "plus", hash[:accession], putative_erv[:putative_erv].first.hit_to, 1..distance, putative_erv[:flanking_3])
        ])
      end
        
      minus_regions = hash[:flanking_pairs][:minus].inject([]) do |array, putative_erv|
        array.concat([
          formatted_region(hash[:definition], "5' genomic region", "minus", hash[:accession], putative_erv[:putative_erv].first.hit_from, 1..distance, putative_erv[:flanking_5]),
          formatted_region(hash[:definition], "3' genomic region", "minus", hash[:accession], putative_erv[:putative_erv].first.hit_to, -distance..-1, putative_erv[:flanking_3])
        ])
      end
      
      fasta_sequences.concat(plus_regions).concat(minus_regions)
    end
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
  
  def complement(string)
    Bio::Sequence::NA.new(string).complement
  end
  
  def print_solo_ltrs(cell = "%-30s")
    headers = ["Accession number", "Strand", "From", "To", "Identity"]
    
    printf (cell * headers.length + "\n") % headers
    
    filter.each do |hsps|
      print_solo_ltr_group(
        cell,
        headers, 
        (hsps[:solo_ltrs][:plus] + hsps[:solo_ltrs][:minus]).sort_by(&:identity).reverse, 
        hsps[:accession]
      )
    end
    
    return nil
  end
  
  private
  
  def print_solo_ltr_group(cell, headers, solo_ltrs, accession)
    solo_ltrs.each do |solo_ltr|
      printf (cell * headers.length + "\n") % [
        accession,
        same_strand?(solo_ltr) ? "Plus" : "Minus",
        solo_ltr.hit_from, 
        solo_ltr.hit_to,
        solo_ltr.identity
      ]
    end
  end
  
  def formatted_region(definition, region, strand, accession, position, window, sequence)
    {
      comment: "> %s (%s on %s strand: %s - coordinates %s to %s" % [
        definition,
        region,
        strand,
        accession,
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
          flanking_5:  complement(sequence_from_entrez(id, putative_erv.last.hit_from, 1..FLANKING_DISTANCE)),
          flanking_3:  complement(sequence_from_entrez(id, putative_erv.first.hit_to, -FLANKING_DISTANCE..-1))
        }
      end
    }
  end
  
  def filter_array(hsps)
    ervs = hsps.combination(2).select do |match_5, match_3|
      ERV_DISTANCE.include?(match_3.hit_from - match_5.hit_to)
    end
    
    [ervs, hsps.reject { |hsp| ervs.flatten.include?(hsp) }]
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

# parser = BlastParser.bootstrap("/Users/evansenter/Documents/School/BC/Rotations/Rotation 3 - Johnson/Canis Lupus Familiaris/Candidate Sequences/LTR against CanFam3.1.xml")
# parser.flanking_regions_as_fasta!("/Users/evansenter/Desktop/genomic_regions.fasta")
# parser.print_solo_ltrs