'''This program maps Ensembl Gene IDs to Gene Names.

It takes two tables (xref.txt.gz and gene.txt.gz) which can be found at
ftp://ftp.ensembl.org/pub/release-XX/mysql/homo_sapiens_funcgen_XX_37/
(XX is version number)

The output (ensg_gene.list) is a list of mapping between ENSG IDs and gene name.
'''

import sys
import gzip


def map_ensg_gene():

    xref_gzip = 'xref.txt.gz'
    gene_gzip = 'gene.txt.gz'
    out_list = 'ensg_gene.list'

    xref_dict = dict()
    print('reading ' + xref_gzip, file=sys.stderr)
    for line in gzip.open(xref_gzip):
        data = line.decode('ascii').rstrip().split('\t')
        xref = data[0]
        gene = data[3]
        xref_dict[xref] = gene

    ensg_dict = dict()
    print('reading ' + gene_gzip, file=sys.stderr)
    for line in gzip.open(gene_gzip):
        data = line.decode('ascii').rstrip().split('\t')
        xref = data[7]
        ensg = data[13]
        ensg_dict[ensg] = xref_dict[xref]

    print('writing ' + out_list, file=sys.stderr)
    fout = open(out_list, 'w')
    for ensg, gene in sorted(ensg_dict.items()):
        print(ensg, gene, sep='\t', file=fout)
    fout.close()

map_ensg_gene()
