'''This module deals with reference fasta file.


TODO:
- remove testing part
- add command line query interface
'''

class Fasta:


    def __init__(self, fasta_name):
        '''Create a fasta object.'''

        self.fasta_name = fasta_name
        self.fasta_index_name = fasta_name + '.fai'
        self.index = self.load_index()


    def load_index(self):
        '''Load fasta index (.fai file)

        Returns:
           a dictionary with chromosome name (str) as key, and a dict as value.
           The value dict has keys:
              chrom_len:  chromosome length (number of bases)
              offset:  number of bytes offset from the file beginning
              nbase:  number of bases per line
              nbyte:  number of bytes per line
        '''

        index = {}
        for line in open(self.fasta_index_name):
            chrom, chrom_len, offset, nbase, nbyte = line.strip().split('\t')
            index[chrom] = {
                'chrom_len':int(chrom_len),
                'offset':int(offset),
                'nbase':int(nbase),
                'nbyte':int(nbyte)
            }
        return index


    def query(self, chrom, pstart=None, pend=None):
        '''Query fasta sequence by chrom(, pstart, pend).

        Args:
            chrom (str):  chromosome name
            pstart (int):  starting position of query region (1-based inclusive)
            pend (int): ending position of query region (1-based inclusive)

        Returns:
            sequence (str)

        Raises:
            ValueError  if input range is invalid
        '''

        if chrom not in self.index:
            raise ValueError('query chromosome {} is not found'.format(chrom))

        if pstart and pstart > self.index[chrom]['chrom_len']:
            raise ValueError('query region outside of chromosome {}:{}-'
                             .format(chrom, pstart))

        if pend and pend < 1:
            raise ValueError('query region outside of chromosome {}:-{}'
                             .format(chrom, pend))

        if pstart and pend and pstart > pend:
            raise ValueError('query region is illegal {}:{}-{}'
                             .format(chrom, pstart, pend))

        # set default start/end position if not specified
        if not pstart or pstart < 1:
            pstart = 1
        if not pend or pend > self.index[chrom]['chrom_len']:
            pend = self.index[chrom]['chrom_len']

        bstart = self._calc_byte_offset(chrom, pstart)  # 0-based inclusive
        bend = self._calc_byte_offset(chrom, pend)  # 0-based inclusive
        sequence = ''
        with open(self.fasta_name) as f:
            f.seek(bstart)
            sequence = f.read(bend - bstart + 1)
        sequence = sequence.replace('\n','')

        return sequence


    def query_region(self, region):
        '''Query fasta sequence by a region string

        Args:
            region (str):  region of format chrom(:pstart-pend)

        Returns:
            sequence (str)

        Raises:
            ValueError  if input range is invalid
        '''

        if ':' not in region:
            return self.query(chrom)
        chrom, prange = region.split(':')
        pstart, pend = map(int, prange.split('-'))
        return self.query(chrom, pstart, pend)


    def _calc_byte_offset(self, chrom, pos):
        '''Calculate byte offset (0-based) of a genomic position in fasta file.

        Args:
            chrom (str):  query chromosome name
            pos (int):  query base position (1-based)

        Returns:
            offset (int):  offset in fasta file (0-based), to be used by seek()
        '''

        nbase = self.index[chrom]['nbase']
        nbyte = self.index[chrom]['nbyte']
        nline = (pos-1) // nbase
        residual = (pos-1) % nbase
        offset = self.index[chrom]['offset'] + nline * nbyte + residual

        return offset


def samtools_query(fasta_name, region):
    '''Samtools fasta query wrapper.

    Args:
        fasta_name (str):  query fasta filename
        region (str): query region

    Returns:
        sequence (str):  sequence in the fasta
    '''

    from subprocess import Popen, PIPE

    cmd = 'samtools faidx {} {}'.format(fasta_name, region)
    proc = Popen(cmd, stdout=PIPE, shell=True, universal_newlines=True)
    result = []
    for line in proc.stdout:
        result.append(line.rstrip())
    proc.communicate()
    sequence = ''.join(result[1:])

    return sequence


def test(region, fasta_name):

    fasta_name = fasta_name
    fa = Fasta(fasta_name)

    seq1 = samtools_query(fasta_name, region)
    seq2 = fa.query_region(region)

    if seq1 != seq2:
        print(region, 'failed', flush=True)
    # else:
    #     print(region, 'success!')


if __name__ == '__main__':

    import sys
    import random

    fasta_name = sys.argv[1]
    test('20:11234133-11234133', fasta_name)

    chrom_set = [str(i) for i in range(1,23)] + ['X','Y','MT']

    for i in range(100):
        chrom = random.choice(chrom_set)
        chrom = '2'
        pstart = random.randrange(1,100000000)
        pend = pstart + random.randrange(1,1000)
        region = '{}:{}-{}'.format(chrom, pstart, pend)
        if (i+1) % 10 == 0:
            print('> tested {} cases'.format(i+1), flush=True)
        test(region, fasta_name)
