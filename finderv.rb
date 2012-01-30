require "bio"
require "./putative_erv.rb"
require "./ltr.rb"

class Finderv
  ERV_DISTANCE = 200..12_000
  
  attr_reader :report
  
  class << self
    def bootstrap(path)
      new(Bio::Blast::Report.new(File.read(path)))
    end

    def solo_ltrs_from_batch(files, ignore_regex = nil)
      parsers   = files.map(&method(:bootstrap))
      solo_ltrs = parsers.map { |parser| parser.solo_ltrs(ignore_regex) }.inject(&:concat)

      curate_similar(solo_ltrs)
    end

    def write_solo_ltr_batch!(solo_ltrs, directory)
      solo_ltrs.each do |solo_ltr|
        solo_ltr.write_fasta!(directory)
      end
    end
    
    def putative_ervs_from_batch(files)
      parsers = files.map(&method(:bootstrap))
      ervs    = parsers.map(&:putative_ervs).inject(&:concat)
      
      curate_similar(ervs)
    end
    
    def write_putative_erv_batch!(putative_ervs, directory)
      putative_ervs.each do |putative_erv|
        putative_erv.write_fasta!(directory)
      end
    end
    
    def curate_similar(elements)
      grouped_elements = elements.inject([]) do |grouped_elements, element|
        grouped_elements.tap do
          if matching_element_group = grouped_elements.find { |element_group| element_group.first == element }
            matching_element_group << element
          else
            grouped_elements << [element]
          end
        end
      end

      grouped_elements.map { |element_group| element_group.sort_by(&:length).last }
    end
  end
  
  def initialize(report)
    @report = report
  end
  
  def filtered_ltrs(ignore_regex = nil)
    ltrs = report.hits.reject { |hit| hit.definition =~ ignore_regex if ignore_regex }.map do |hit|
      parse_hsps(hit)
    end.flatten
    
    curate_similar(ltrs)
  end
  
  def putative_ervs
    ervs = report.hits.map do |hit|
      parse_hsps(hit).partition(&:plus_strand?).map do |ltrs|
        find_ervs(ltrs)
      end
    end.flatten
    
    curate_similar(ervs)
  end
  
  def solo_ltrs(ignore_regex = nil)
    ervs = putative_ervs.map(&:ltrs).inject(&:concat)
    
    filtered_ltrs(ignore_regex).reject { |ltr| ervs.include?(ltr) }
  end
  
  private
  
  def curate_similar(elements)
    self.class.curate_similar(elements)
  end
  
  def find_ervs(ltrs)
    ltrs.combination(2).select do |match_5, match_3|
      ERV_DISTANCE.include?(match_3.up_coord - match_5.down_coord)
    end.map do |ltr_pair|
      PutativeErv.new(*ltr_pair)
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

# files  = Dir["/Users/evansenter/Documents/School/BC/Rotations/Rotation 3 - Johnson/Kate's Stuff/*.xml"]
# parser = Finderv.bootstrap("/Users/evansenter/Documents/School/BC/Rotations/Rotation 3 - Johnson/Canis Lupus Familiaris/Candidate Sequences/LTR against CanFam3.1.xml")
# parser = Finderv.bootstrap("/Users/evansenter/Documents/School/BC/Rotations/Rotation 3 - Johnson/Kate's Stuff/J34ZR93001R-Alignment.xml")