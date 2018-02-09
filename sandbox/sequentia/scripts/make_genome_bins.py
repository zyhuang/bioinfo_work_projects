'''This script prints a list of all 1MB genomic bin ids (chrom.nnn).

e.g. 1.000 ... 1.249, 2.000 ... 2.051 ... X.000 ... Y.000 ... MT.000

For hg19, in total, there are 3114 bins,
  3053 bins on chr1-22 and chrX,
    60 bins on chrY
     1 bin on chrMT.
'''

import sys

def make_genome_bins(refgen_fa):
    '''print a list of 1MB-bin ids to stdout.

    Args:
        refgen_fa (str)  reference genome fasta file name

    Return:
        stdout
    '''

    chrom_set = set([str(i) for i in range(1,23)] + ['X', 'Y', 'MT'])
    for line in open(refgen_fa + '.fai'):
        chrom, length = line.rstrip().split('\t')[:2]
        if chrom not in chrom_set:
            continue
        length_prefix = length.zfill(9)[:3]
        for i in range(int(length_prefix)+1):
            print('{}.{}'.format(chrom, str(i).zfill(3)))


if __name__ == '__main__':

    refgen_fa = sys.argv[1]
    make_genome_bins(refgen_fa)
