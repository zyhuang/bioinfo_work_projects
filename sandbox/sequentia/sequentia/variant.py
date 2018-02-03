'''This module deals with variant normalization and conversion between various
formats.

VCF and Annovar variant notations are not unique. A variant can have multiple
representations and/or not left-normalized. This module converts these into a
unique representation, which is parsimonious and lossless
(chrom[nochr]:pos[09s]:ref[ACGT]+:alt[ACGT]+).

Annovar data format has two forms:

- 4 columns (chrom,pos,ref,alt tab-delimited). pos is 1-based inclusive.
  REF are [ACGT]+, ALT are [0-9]+, or [ACGT]+, or [0-9]+[ACGT]+

  Examples:
    hg19_ALL.sites.2015_08.txt


- 5 columns (chrom,pstart,pend,ref,alt tab-delimited). pstart and pend are
  1-based inclusive. REF and ALT are [ACGT]+ or -

  Examples:
    hg19_avsnp150.txt
    hg19_esp6500siv2_all.txt
    hg19_exac03nontcga.txt
    hg19_gnomad_exome.txt
    hg19_gnomad_genome.txt


TODO:

'''

import sys
import re


class Variant:
    '''This module deals with variant notation conversion and left-normalization.

    Attributes:
        chrom (str):  chromosome name, without chr or Chr
        pos (int):  reference position
        ref (str):  reference allele [ACGT]+
        alt (str):  alternate allele [ACGT]+
        left_norm (bool):  if the variant is left-normalized
        type (str):   variant type: SNV, Insertion, Deletion, or MNV

    Methods:
        __str__():  convert the variant into string (varkey)
        left_normalize():  apply left-normalization (left-alignment and remove
            redundant bases on both ends
        var_type():  determine the variant type (SNV/Insertion/Deletion/MNV)
        parse():  parse input information of a variant (chrom, pos, ref, alt)

    '''

    def __init__(self):
        '''Initiliaze a variant object.'''

        self.chrom = None
        self.pos = None
        self.ref = None
        self.alt = None
        self.lnorm = False


    def __str__(self):
        '''Format variant into colon-seperated normalized format.

        Update:
            self.normalized = True
            set varkey member
        '''

        if (self.chrom == None and self.pos == None and self.ref == None and
            self.alt == None):
            return str(None)
        else:
            return (
                '{}:{}:{}:{}'.format(self.chrom, str(self.pos).zfill(9),
                                     self.ref, self.alt))


    def left_normalize(self):
        '''Remove redundant bases among REF and ALT and left-align ALT allele.

        '''
        # no need to left_norm SNP
        if len(self.ref) == 1 or len(self.alt) == 1:
            self.lnorm = True
            return

        # already left normalized
        if self.ref[0] != self.alt[0] and self.ref[-1] != self.alt[-1]:
            self.lnorm = True
            return

        # trim redundant bases on the right-end
        if self.ref[-1] == self.alt[-1]:
            i = 1
            while (i != min(len(self.ref), len(self.alt)) and
                   self.ref[-i] == self.alt[-i]):
                i += 1
            if i > 1:
                self.ref = self.ref[:-i+1]
                self.alt = self.alt[:-i+1]

        # if left-normalized after trim
        if len(self.ref) == 1 or len(self.alt) == 1:
            self.lnorm = True
            return

        # trim redundant bases on the left-end (except for the last redundant one)
        if self.ref[0] == self.alt[0]:
            i = 0
            while (i != min(len(self.ref), len(self.alt)) - 1 and
                   self.ref[i] == self.alt[i]):
                i += 1
            self.ref = self.ref[i:]
            self.alt = self.alt[i:]
            self.pos += i

        self.lnorm = True



    def var_type(self):
        '''Determine variant type: SNV, Insertion, Deletion, MNV

        Returns:
           variant type (str):  SNV/Insertion/Deletion/MNV
        '''

        if not self.lnorm:
            self.left_normalize()

        if len(self.ref) == 1 and len(self.alt) == 1:
            return 'SNV'

        elif len(self.ref) == 1:
            if self.ref != self.alt[0]:
                return 'MNV'
            else:
                return 'Insertion'

        elif len(self.alt) == 1:
            if self.alt != self.ref[0]:
                return 'MNV'
            else:
                return 'Deletion'

        else:
            return 'MNV'


    def parse(self, chrom, pos, ref, alt, refgen_fa, left_norm=False):
        '''parse variant given chrom, pos, ref, alt.

        Args:
            chrom:  chromosome name (chr will be stripped)
            pos:  1-based position in the reference
            ref:  reference allele (can be -, [ACGT]+)
            alt:  alternate allele (can be -, [ACGT]+, [0-9]+, [0-9]+[ACGT]+)
            refgen_fa:  Fasta object, see Fasta module

        Attributes:
            member data: chrom, pos, ref, alt are updated
            member attribute left_norm can optionally be set

        Raises:
            ValueError: when
               ALT is CNV or SV (in form of <XXX>) or multiallelic,
               REF is none of "-", or [ACGT]+
               ALT is none of "-", [0-9]+, [ACGT]+ or [0-9]+[ACGT]+
        '''

        # trim "chr" if exists
        if chrom[:3] in {'chr', 'Chr'}:
            self.chrom = chrom[3:]
        else:
            self.chrom = chrom

        # convert pos to integer
        try:
            pos = int(pos)
        except:
            raise ValueError(
                'Input position must be integer (variant = {}:{}:{}:{})'
                .format(chrom, pos, ref, alt))

        # for the most of cases (SNV)
        if len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT":
            self.pos = pos
            self.ref = ref
            self.alt = alt
            if left_norm:
                self.lnorm = True
            return

        # check illegal REF or ALT
        if alt.startswith('<') or alt.endswith('>'):
            raise ValueError(
                'CNV/SV is not supported (variant = {}:{}:{}:{})'
                .format(chrom, pos, ref, alt))

        if ',' in alt:
            raise ValueError(
                'Multiallelic ALT is not supported (variant = {}:{}:{}:{})'
                .format(chrom, pos, ref, alt))

        if (re.sub(r'[ACGT-]+', '', ref) != '' or
            re.sub(r'[0-9ACGT-]+', '', alt) != ''):
            raise ValueError(
                'REF should be "[ACGT]+" or "-", ALT should be "[0-9]+[ACGT]+"'
                ' or "-" (variant = {}:{}:{}:{})'
                .format(chrom, pos, ref, alt))

        # Indel and MNV
        if ref == '-' and alt != '-':
            # insertion, e.g.  -:ACG => R:RACG
            self.pos = pos
            self.ref = refgen_fa.query(self.chrom, pstart=self.pos,
                                       pend=self.pos)
            self.alt = self.ref + alt

        elif alt == '-' and ref != '-':
            # deletion, e.g. ACG:- => RACG:R
            self.pos = pos - 1
            self.alt = refgen_fa.query(self.chrom, pstart=self.pos,
                                       pend=self.pos)
            self.ref = self.alt + ref

        elif alt[0] in "ACGT":
            # [ACGT]+ indel or MNV
            self.pos = pos
            self.ref = ref
            self.alt = alt

        elif alt[-1].isdigit():
            # [0-9]+, alt is a number, =deletion, e.g. AGGG:4 => RAGGG:R
            self.pos = pos - 1
            self.alt = refgen_fa.query(self.chrom, pstart=self.pos,
                                       pend=self.pos)
            self.ref = self.alt + ref

        elif alt[0].isdigit() and alt[-1] in "ACGT":
            # [0-9]+[ACGT]+, indel or MNV
            token = re.findall(r'([0-9]+)([ACGT]+)', alt)
            if len(token) != 1:
                raise ValueError(
                    'input ALT is malformatted, expect ([0-9]+)([ACGT]+) '
                    '(variant = {}:{}:{}:{})'
                    .format(chrom, pos, ref, alt))

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

        if left_norm:
            self.left_normalize()


def sort_varkey(varkey):

    chrom_hash = {'X':23, 'Y':24, 'MT':25}

    chrom, pos, ref, alt = varkey.split(':')
    if chrom in chrom_hash:
        chrom = chrom_hash[chrom]
    else:
        chrom = int(chrom)
    pos = int(pos)

    return chrom, pos, ref, alt
