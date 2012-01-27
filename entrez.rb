require "entrez"

module EntrezSequence
  def self.included(base)
    base.send(:include, InstanceMethods)
  end
  
  module InstanceMethods
    def na_sequence_from_entrez(id, position, window)
      Bio::Sequence::NA.new(sequence_from_entrez(id, position, window).seq)
    end
    
    def aa_sequence_from_entrez(id, position, window)
      Bio::Sequence::AA.new(sequence_from_entrez(id, position, window).seq)
    end
    
    def sequence_from_entrez(id, position, window)
      fasta = Entrez.EFetch("nuccore", {
        id:        id, 
        seq_start: position + window.min, 
        seq_stop:  position + window.max, 
        retmode:   :fasta, 
        rettype:   :text
      }).response.body

      Bio::FastaFormat.new(fasta)
    end
  end
end