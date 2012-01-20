require "bio"

module FastaParser
  def self.included(base)
    base.extend(ClassMethods)
  end
  
  module ClassMethods
    def parse_fasta_file(path)
      data = File.read(path)
      data.split(/\n\n/).map { |fasta_sequence| Bio::FastaFormat.new(fasta_sequence) }
    end
  end
end

# path = "/Users/evansenter/Documents/School/BC/Rotations/Rotation 3 - Johnson/Canis Lupus Familiaris/putative_ervs.fasta"