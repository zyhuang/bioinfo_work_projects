'''This module deals with variant normalization and conversion between various
formats.

VCF and Annovar variant notations are not unique. A variant can have multiple
representations and/or not left-normalized. This module converts these into a
unique representation, which is parsimonious and lossless
(chrom[nochr]:pos[09]:ref[1]:alt[+-ACGT]).

Annovar data format has two forms:

- 5 columns (chrom,pstart,pend,ref,alt tab-delimited). pstart and pend are
  1-based inclusive. ref and alt may be "-" in case of indel.

  Examples:
    hg19_avsnp150.txt
    hg19_esp6500siv2_all.txt
    hg19_exac03nontcga.txt
    hg19_gnomad_exome.txt
    hg19_gnomad_genome.txt

- 4 columns (chrom,pos,ref,alt tab-delimited). pos is 1-based inclusive.

  Examples:
    hg19_ALL.sites.2015_08.txt

'''

import sys
from . import fasta


class Variant:

    def __init__(self):

    	self.chrom = None
    	self.pos = None
    	self.ref = None
    	self.alt = None
    	self.info = {}
        self.normalized = False
        self.varkey = None


    def format_varkey(self):
        '''Format variant into colon-seperated normalized format.

        Update:
            self.normalized = True
            set varkey member
        '''

    	self.varkey = ('{}:{}:{}:{}'
    		       .format(self.chrom, str(self.pos).zfill(9), self.ref, self.alt))


    def normalize(self):
        pass


    def parse_annovar_format_5col(chrom, pstart, pend, ref, alt):

















def parse_vcf_record(vcf_line):
