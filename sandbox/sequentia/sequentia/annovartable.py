'''This module deals with query of annovar table.

'''
import sys

class AnnovarTable:

    def __init__(self, table_name):
        '''Create a table object.'''

        self.table_name = table_name
        self.table_index_name = table_name + '.idx'
        self.index = {}
        self.bin_size = 0
        self.file_size = 0
        self.load_index()


    def load_index(self):
        '''Load table index (.idx file)

        Set Attributes:
            self.index (dict):  stores offset in byte of each bin position
                {chrom: {pos: {offset_start: offset, offset_end: offset}}}
                offset_start is 0-based inclusive (used by seek())
                offset_end is 0-based non-inclusive
                (used in nbyte = offset_end - offset_start)
            self.binsize (int):  number of bases per bin
            self.filesize (int):  total number of bytes in the table file
        '''

        self.index = {}
        print('loading ' + self.table_index_name, file=sys.stderr)
        for line in open(self.table_index_name):
            if line.startswith('#'):
                foo, binsize, filesize = line.rstrip().split('\t')
                self.bin_size = int(binsize)
                self.file_size = int(filesize)
                continue
            chrom, pos, bstart, bend = line.rstrip().split('\t')
            pos = int(pos)
            if chrom not in self.index:
                self.index[chrom] = {}
            self.index[chrom][pos] = {
                'offset_start': int(bstart),
                'offset_end': int(bend),
            }


    def query_region(self, region=None):
        '''Query annovar table by a region string

        Args:
            region (str):  region of format chrom(:pstart-pend)

        Returns:
           generator of lines overlapping with query region.
        '''

        if region == None:
            return self.query()

        if ':' not in region:
            return self.query(chrom=region)

        chrom, prange = region.split(':')
        pstart, pend = map(int, prange.split('-'))
        return self.query(chrom=chrom, pstart=pstart, pend=pend)


    def query(self, chrom=None, pstart=None, pend=None):
        '''Query table by chrom(, pstart, pend).

        Args:
            chrom (str):  chromosome name
            pstart (int):  starting position of query region (1-based inclusive)
            pend (int): ending position of query region (1-based inclusive)

        Returns:
           generator of lines overlapping with query region.
        '''

        if chrom == None:
            with open(self.table_name) as f:
                for line in f:
                    yield line.rstrip()
            return

        chrom_len = sorted(self.index[chrom])[-1] + self.bin_size - 1

        if chrom not in self.index:
            print('*WARNING* query chromosome {} is not found'.format(chrom),
                  file=sys.stderr)
            yield from []
            return

        if pstart and pstart > chrom_len:
            print('*WARNING* query region outside of chromosome {}:{}-'
                  .format(chrom, pstart), file=sys.stderr)
            yield from []
            return

        if pend and pend < 1:
            print('*WARNING* query region outside of chromosome {}:-{}'
                  .format(chrom, pend), file=sys.stderr)
            yield from []
            return

        if pstart and pend and pstart > pend:
            print('*WARNING* query region is illegal {}:{}-{}'
                  .format(chrom, pstart, pend), file=sys.stderr)
            yield from []
            return

        # set default start/end position if not specified
        if not pstart or pstart < 1:
            pstart = 1
        if not pend or pend > chrom_len:
            pend = chrom_len

        # find overlapping bins
        bin_start = self._find_bin(chrom, pstart, forward=True)
        bin_end = self._find_bin(chrom, pend, forward=False)

        # specified query range is either before the first bin or after the last bin
        if bin_start == None or bin_end == None:
            yield from []
            return

        offset_start = self.index[chrom][bin_start]['offset_start']
        offset_end = self.index[chrom][bin_end]['offset_end']

        with open(self.table_name) as f:
            f.seek(offset_start)
            data = f.read(offset_end-offset_start).split('\n')
            data.pop()
            for line in data:
                lpos = int(line.split('\t')[1])
                if lpos >= pstart and lpos <= pend:
                    yield line

        return


    def _find_bin(self, chrom, pos, forward=True):
        '''Find the bin position of query position.

        Example:
        1..99 => 0, 100..199 => 100, 200..299 => 200

        if bin is not found, then find the next existing bin
        (diretion='forward') or prvious existing bin (direction='backward')

        if next existing bin is after the last bin, return None
        if previous existing bin is before the before bin, return None
        '''

        pos_floor = pos // self.bin_size * self.bin_size
        if pos_floor in self.index[chrom]:
            return pos_floor

        bin_list = sorted(self.index[chrom])

        if forward:
            if pos_floor > bin_list[-1]:
                return None
            for bin_pos in range(pos_floor, bin_list[-1]+1, self.bin_size):
                if bin_pos in self.index[chrom]:
                    return bin_pos

        else:  # backward
            if pos_floor < bin_list[0]:
                return None
            for bin_pos in range(pos_floor, bin_list[0]-1, -self.bin_size):
                if bin_pos in self.index[chrom]:
                    return bin_pos

def test():

    table_name = '../scripts/humandb/hg19_esp6500siv2_aa.txt'
    table = AnnovarTable(table_name)
#    for line in table.query(chrom='X', pstart=154001557, pend=154001557+100000):
    for line in table.query(chrom='X', pstart=155239822, pend=155239900):
    # for line in table.query(chrom='X', pstart=155239899):
    # for line in table.query():
        print(line)


if __name__ == '__main__':

    table_name, region_str = sys.argv[1:3]
    table = AnnovarTable(table_name)
    for line in table.query_region(region_str):
        print(line)
