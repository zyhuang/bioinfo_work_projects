'''This module deals with variant normalization and conversion between various
formats.

VCF and Annovar variant notations are not unique. A variant can have multiple
representations and/or not left-normalized. This module converts these into a
unique representation, which is parsimonious and lossless
(chrom[nochr]:pos[09s]:ref[ACGT]+:alt[ACGT]+).

Annovar data format has two forms:

- 4 columns (chrom,pos,ref,alt tab-delimited). pos is 1-based inclusive.

  Examples:
    hg19_ALL.sites.2015_08.txt


- 5 columns (chrom,pstart,pend,ref,alt tab-delimited). pstart and pend are
  1-based inclusive. ref and alt may be "-" in case of indel.

  Examples:
    hg19_avsnp150.txt
    hg19_esp6500siv2_all.txt
    hg19_exac03nontcga.txt
    hg19_gnomad_exome.txt
    hg19_gnomad_genome.txt

TODO:

    add fasta reference base query
'''

import sys
import re
# from . import fasta


class Variant:
    '''This module deals with variant notation conversion and left-normalization.

    '''

    def __init__(self):

        self.chrom = None
        self.pos = None
        self.ref = None
        self.alt = None
        self.left_norm = False


    def __str__(self):
        '''Format variant into colon-seperated normalized format.

        Update:
            self.normalized = True
            set varkey member
        '''

        if self.chrom == None and self.pos == None:
            return str(None)
        else:
            return (
                '{}:{}:{}:{}'.format(self.chrom, str(self.pos).zfill(9),
                                     self.ref, self.alt))


    def left_normalize(self):
        '''Remove redundant bases among REF and ALT and left-align ALT allele.

        '''
        if len(self.ref) == 1 or len(self.alt) == 1:
            self.left_norm = True
            return

        if self.ref[0] != self.alt[0] and self.ref[-1] != self.alt[-1]:
            self.left_norm = True
            return

        print(self.ref)
        print(self.alt)

        if self.ref[-1] == self.alt[-1]:
            i = 1
            while (i != min(len(self.ref), len(self.alt)) and
                   self.ref[-i] == self.alt[-i]):
                i += 1
            if i != 1:
                self.ref = self.ref[:-i+1]
                self.alt = self.alt[:-i+1]

        print(self.ref)
        print(self.alt)

        if len(self.ref) == 1 or len(self.alt) == 1:
            self.left_norm = True
            return


        if self.ref[0] == self.alt[0]:
            i = 0
            while (i != min(len(self.ref), len(self.alt)) - 1 and
                   self.ref[i] == self.alt[i]):
                i += 1
            self.ref = self.ref[i:]
            self.alt = self.alt[i:]
            self.pos += i

        print(self.ref)
        print(self.alt)

        self.left_norm = True



    def var_type(self):

        if not self.left_norm:
            self.left_nromalize()
        # if SNV, INDEL, MNV...


    def parse_vcf(self, chrom, pos, ref, alt, left_norm=False):

        # deal with ALT == '*'
        parse_annovar(chrom, pstart, ref, alt, pend=pstart,
                      left_norm=left_norm)
        pass


    def parse_annovar(self, chrom, pstart, ref, alt, pend=None, left_norm=False):
        '''parse first 4 or 5 cols from annovar table.

        Args:
            chrom:  chromosome name (chr will be stripped)
            pstart (pend):  1-based start and end positions in the reference
            ref:  reference allele (can be ACGT, -)
            alt:  alternate allele (can contain <XXX>, -, nnnACGT, nnn)

        Attributes:
            member data: chrom, pos, ref, alt are updated

        Raises:
            ValueError: if ALT is CNV or SV (in form of <XXX>) or multiallelic,
                        or in 4-col format, ALT is not [0-9]+[ACGT]+
        '''

        if alt.startswith('<') or alt.endswith('>'):
            raise ValueError(
                'CNV/SV is not supported (variant = {}:{}:{}:{})'
                .format(chrom, pstart, ref, alt))

        if ',' in alt:
            raise ValueError(
                'Multiallelic ALT is not supported (variant = {}:{}:{}:{})'
                .format(chrom, pstart, ref, alt))

        # preparation
        if chrom[:3] in {'chr', 'Chr'}:
            self.chrom = chrom[3:]
        else:
            self.chrom = chrom

        try:
            pos = int(pstart)
        except:
            raise ValueError(
                'input position must be integer (variant = {}:{}:{}:{})'
                .format(chrom, pstart, ref, alt))

        # for the most of cases (SNV)
        if len(ref) == 1 and len(alt) == 1:
            self.pos = pos
            self.ref = ref
            self.alt = alt
            if left_norm:
                self.left_norm = True
            return

        # 4col version, alt with ([0-9]+) ([ACGT]+) or both notation
        if pend == None:

            if alt[-1].isdigit():
                # [0-9]+  deletion, 1 255403 AGGG 4
                self.pos = pos - 1
                self.ref = 'R' + ref
                self.alt = 'R'

            elif alt[0] in "ACGT":
                # [ACGT]+  indel or MNV
                self.pos = pos
                self.ref = ref
                self.alt = alt

            else:
                # [0-9]+[ACGT]+, indel or MNV
                token = re.findall(r'([0-9]+)([ACGT]+)', alt)
                if len(token) != 1:
                    raise ValueError(
                        'input ALT is malformatted, expect ([0-9]+)([ACGT]+) '
                        '(variant = {}:{}:{}:{})'
                        .format(chrom, pstart, ref, alt))

                if token[0][0] == '0':
                    # insertion, A:0C => A:AC
                    self.pos = pos
                    self.ref = ref
                    self.alt = ref + token[0][1]

                elif token[0][0] == '1':
                    # insertion, A:1AC => A:AC
                    self.pos = pos
                    self.ref = ref
                    self.alt = token[0][1]

                else:
                    # deletion TTAA:4T => TTAA:T
                    self.pos = pos
                    self.ref = ref
                    self.alt = token[0][1]

        # 5col version, ref or alt with - notation
        else:
            if ref != '-' and alt != '-':
                # MNV or Indel
                self.pos = pos
                self.ref = ref
                self.alt = alt

            elif ref == '-':
                # insertion
                self.pos = pos
                self.ref = 'R'
                self.alt = 'R' + alt

            elif alt == '-':
                # deletion
                self.pos = pos - 1
                self.ref = 'R' + ref
                self.alt = 'R'

        if left_norm:
            self.left_normalize()
