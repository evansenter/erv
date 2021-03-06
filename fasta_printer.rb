require "bio"
require "active_support/core_ext/kernel/reporting"

module FastaPrinter
  def self.included(base)
    base.send(:include, InstanceMethods)
  end
  
  module InstanceMethods
    def fasta_hash
      {
        comment:  fasta_comment,
        sequence: seq
      }
    end

    def fasta_comment
      "> %s ((%s) %s nt. long) %s %s %s" % [
        definition,
        plus_strand? ? "+" : "-",
        length,
        accession,
        from,
        to
      ]
    end

    def fasta_filename
      "%s %s %s %s %s %s.fasta" % [
        definition.split(/,/).first,
        accession,
        plus_strand? ? "plus" : "minus",
        from,
        to,
        type
      ]
    end
    
    def print_fasta
      silence_stream(STDERR) do
        puts fasta_comment
        puts seq
        puts
      end
    end
    
    def write_fasta!(directory)
      path = File.join(directory, fasta_filename)
      
      unless File.exists?(path)
        File.open(path, "w") do |file|
          [fasta_hash].flatten.each do |fasta_sequence|
            file.write(fasta_sequence[:comment])
            file.write("\n")
            file.write(fasta_sequence[:sequence])
            file.write("\n\n")
          end
        end
      end
    end
  end
end