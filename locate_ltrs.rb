require "bio"
require "./putative_erv.rb"
require "./ltr.rb"

class BlastParser
  ERV_DISTANCE = 4_000..12_000
  
  attr_reader :report
  
  class << self
    def bootstrap(path)
      new(Bio::Blast::Report.new(File.read(path)))
    end

    def solo_ltrs_from_batch(files, ignore_regex = nil)
      # Dir["/Users/evansenter/Documents/School/BC/Rotations/Rotation 3 - Johnson/Kate's Stuff/*.xml"]

      parsers   = files.map(&method(:bootstrap))
      solo_ltrs = parsers.map { |parser| parser.solo_ltrs(ignore_regex) }.inject(&:concat)

      curate_similar_ltrs(solo_ltrs)
    end

    def write_solo_ltr_batch!(solo_ltrs, directory)
      solo_ltrs.each do |solo_ltr|
        write_fasta!(File.join(directory, solo_ltr.fasta_filename), solo_ltr.fasta_hash)
      end
    end
    
    def curate_similar_ltrs(solo_ltrs)
      grouped_solo_ltrs = solo_ltrs.inject([]) do |grouped_ltrs, ltr|
        grouped_ltrs.tap do
          if matching_ltr_group = grouped_ltrs.find { |ltr_group| ltr_group.first == ltr }
            matching_ltr_group << ltr
          else
            grouped_ltrs << [ltr]
          end
        end
      end

      grouped_solo_ltrs.map { |ltr_group| ltr_group.sort_by(&:length).last }
    end

    private

    def write_fasta!(path, fasta_content)
      unless File.exists?(path)
        File.open(path, "w") do |file|
          [fasta_content].flatten.each do |fasta_sequence|
            file.write(fasta_sequence[:comment])
            file.write("\n")
            file.write(fasta_sequence[:sequence])
            file.write("\n\n")
          end
        end
      end
    end
  end
  
  def initialize(report)
    @report = report
  end
  
  def filtered_ltrs(ignore_regex = nil)
    ltrs = report.hits.reject { |hit| hit.definition =~ ignore_regex if ignore_regex }.map do |hit|
      parse_hsps(hit)
    end.flatten
    
    self.class.curate_similar_ltrs(ltrs)
  end
  
  def putative_ervs
    report.hits.map do |hit|
      parse_hsps(hit).partition(&:plus_strand?).map do |ltrs|
        find_ervs(ltrs)
      end
    end.flatten
  end
  
  def solo_ltrs(ignore_regex = nil)
    ervs = putative_ervs.map(&:ltrs).inject(&:concat)
    
    filtered_ltrs(ignore_regex).reject { |ltr| ervs.include?(ltr) }
  end
  
  private
  
  def find_ervs(ltrs)
    ltrs.combination(2).select do |match_5, match_3|
      ERV_DISTANCE.include?(match_3.up_coord - match_5.down_coord)
    end.map do |ltr_pair|
      PutativeErv.new(ltr_pair)
    end
  end
  
  def parse_hsps(hit)
    hit.hsps.select do |hsp| 
      hsp.align_len > report.query_len * 0.8
    end.map do |hsp|
      Ltr.new(hit, hsp)
    end.sort_by(&:midpoint)
  end
end

# parser = BlastParser.bootstrap("/Users/evansenter/Documents/School/BC/Rotations/Rotation 3 - Johnson/Canis Lupus Familiaris/Candidate Sequences/LTR against CanFam3.1.xml")
# parser = BlastParser.bootstrap("/Users/evansenter/Documents/School/BC/Rotations/Rotation 3 - Johnson/Kate's Stuff/J34ZR93001R-Alignment.xml")
# parser.flanking_regions_as_fasta!("/Users/evansenter/Desktop/genomic_regions.fasta")
# parser.print_solo_ltrs