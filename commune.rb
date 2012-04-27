require "active_support/inflector"

module Commune
  def self.included(base)
    base.send(:include, InstanceMethods)
  end
  
  module InstanceMethods
    TSD_LENGTH = 5
    
    def type
      "#{self.class.name.titleize.downcase} " + if has_tsd?
        "(pared to TSD #{tsd_5.gsub(/\W/, '_')})"
      else
        "(TSD mismatch 5' #{tsd_5.gsub(/\W/, '_')}, 3' #{tsd_3.gsub(/\W/, '_')})"
      end
    end
    
    def seq
      try_paring_to_tsd
    end
    
    def try_paring_to_tsd
      # In order to identify TSD (target site duplication, insertion points for the retrovirus), we buffer 75nt on either end of
      # the putative ERV (in na_sequence_from_entrez call), and pare it down if the closest flanking 5'/3' 5nt are identical.
    
      if has_tsd?
        raw_seq[(buffer_size - TSD_LENGTH)...(-buffer_size + TSD_LENGTH)]
      else
        raw_seq
      end
    end
    
    def tsd_5
      raw_seq[buffer_size - TSD_LENGTH, TSD_LENGTH]
    end
    
    def tsd_3
      raw_seq[-buffer_size, TSD_LENGTH]
    end
  
    def has_tsd?
      tsd_5 == tsd_3
    end
  end
end