'''This module deals with query of annovar table.

'''

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
                {chrom: {pos: {offset: offset, nbyte: nbyte}}}
            self.binsize (int):  number of bases per bin
            self.filesize (int):  total number of bytes in the table file
        '''

        self.index = {}
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
                'offset': int(bstart),
                'nbyte': int(bend) - int(bstart)
            }


    def query_region(self, region):
        '''Query annovar table by a region string

        Args:
            region (str):  region of format chrom(:pstart-pend)

        Returns:
            iterator

        Raises:
            ValueError  if input range is invalid
        '''

        if ':' not in region:
            return self.query(chrom)
        chrom, prange = region.split(':')
        pstart, pend = map(int, prange.split('-'))
        return self.query(chrom, pstart, pend)


    def query(self, chrom=None, pstart=None, pend=None):
        '''Query table by chrom(, pstart, pend).

        Args:
            chrom (str):  chromosome name
            pstart (int):  starting position of query region (1-based inclusive)
            pend (int): ending position of query region (1-based inclusive)

        Returns:
           iterator of lines overlapping with query region.

        Raises:
            ValueError  if input range is invalid
        '''

        if chrom == None:
            with open(self.table_name) as f:
                for line in f:
                    yield line.rstrip()
            return

        chrom_len = sort(self.index[chrom])[-1] + self.bin_size

        if chrom not in self.index:
            raise ValueError('query chromosome {} is not found'.format(chrom))

        if pstart and pstart > chrom_len:
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
        if not pend or pend > chrom_len:
            pend = chrom_len

        # find overlapping bins
