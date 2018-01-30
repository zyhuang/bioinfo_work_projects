class Fasta:


    def __init__(self, fasta_name):
        self.fasta_name = fasta_name
        self.fasta_index_name = fasta_name + '.fai'
        self.index = self.load_index()


    def load_index(self):
        index = {}
        for line in open(self.fasta_index_name):
            chrom, chrom_len, offset, nbase, nbyte = line.strip().split('\t')
            index[chrom] = {'chrom_len':int(chrom_len), 'offset':int(offset),
                            'nbase':int(nbase), 'nbyte':int(nbyte)}
        return index


    def query(self, chrom, pstart=None, pend=None):
        ''' pstart, pend are 1-based inclusive '''
        if chrom not in self.index:
            return '', ''

        if pstart and pstart > self.index[chrom]['chrom_len']:
            return '', ''

        if pend and pend < 1:
            return '', ''

        if pstart and pend and pstart > pend:
            return '', ''

        if not pstart or pstart < 1:
            pstart = 1
        if not pend or pend > self.index[chrom]['chrom_len']:
            pend = self.index[chrom]['chrom_len']

        region = '{}:{}-{}'.format(chrom, pstart, pend)

        bstart = self._calc_byte_offset(chrom, pstart)  # 0-based inclusive
        bend = self._calc_byte_offset(chrom, pend)  # 0-based inclusive
        sequence = ''
        with open(self.fasta_name) as f:
            f.seek(bstart)
            sequence = f.read(bend - bstart + 1)
        sequence = sequence.replace('\n','')

        return sequence, region


    def query_region(self, region):
        if ':' not in region:
            return self.query(chrom)
        chrom, prange = region.split(':')
        pstart, pend = map(int, prange.split('-'))
        return self.query(chrom, pstart, pend)


    def _calc_byte_offset(self, chrom, pos):

        nbase = self.index[chrom]['nbase']
        nbyte = self.index[chrom]['nbyte']
        nline = (pos-1) // nbase
        residual = (pos-1) % nbase
        offset = self.index[chrom]['offset'] + nline * nbyte + residual

        return offset


def samtools_query(fasta_name, region):

    from subprocess import Popen, PIPE

    cmd = 'samtools faidx {} {}'.format(fasta_name, region)
    proc = Popen(cmd, stdout=PIPE, shell=True, universal_newlines=True)
    result = []
    for line in proc.stdout:
        result.append(line.rstrip())
    proc.communicate()
    region = result[0][1:]
    sequence = ''.join(result[1:])

    return sequence, region


def test(region):

    fasta_name = '/stornext/snfs4/1000-gen/DATA/REF/human/GRCh37/human.g1k.v37.fa'
    fa = Fasta(fasta_name)

    seq1, reg1 = samtools_query(fasta_name, region)
    seq2, reg2 = fa.query_region(region)

    if seq1 != seq2:
        print(region, 'failed', flush=True)
    # else:
    #     print(region, 'success!')


if __name__ == '__main__':

    test('20:11234133-11234133')

    import random
    chrom_set = [str(i) for i in range(1,23)] + ['X','Y','MT']

    for i in range(100):
        chrom = random.choice(chrom_set)
        pstart = random.randrange(1,100000000)
        pend = pstart + random.randrange(1,100000)
        region = '{}:{}-{}'.format(chrom, pstart, pend)
        if (i+1) % 10 == 0:
            print('> tested {} cases'.format(i+1), flush=True)
        test(region)
