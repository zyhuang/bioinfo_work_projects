'''This program creates a list of ENSG-ID, type and gene name.

It merges two Ensemble table gene.txt.gz and gene_attrib.txt.gz downloaded from
wget ftp://ftp.ensembl.org/pub/current_mysql/homo_sapiens_core_91_38/

The output file is named "ensid.type.name.txt". If ENSG-ID is not present in
"gene_attrib.txt.gz", the gene name is set as "Null".
'''

import os
import sys
import gzip


def get_gene_name():

    gene_data = 'gene.txt.gz'
    gene_attrib_data = 'gene_attrib.txt.gz'
    fout_name = 'ensid.type.name.txt'

    for fname in [gene_data, gene_attrib_data]:
        if not os.path.isfile(fname):
            os.system('wget ftp://ftp.ensembl.org/pub/current_mysql/homo_sapiens_core_91_38/' + fname)

    ensid_dict = dict()
    print('reading ' + gene_data)
    for line in gzip.open(gene_data):
        data = line.decode('ascii').rstrip().split('\t')
        key, gtype, ensid = data[0], data[1], data[-4]
        ensid_dict[key] = {'ensid': ensid, 'type': gtype, 'name': 'Null'}

    ensid_gene_dict = dict()
    print('reading ' + gene_attrib_data)
    for line in gzip.open(gene_attrib_data):
        data = line.decode('ascii').rstrip().split('\t')
        key, number, name = data
        if key in ensid_dict and number == '4':
            ensid_dict[key]['name'] = name

    fout = open(fout_name, 'w')
    print('writing output list ' + fout_name)
    for key,value in sorted(ensid_dict.items(), key=lambda x: x[1]['ensid']):
        print(value['ensid'], value['type'], value['name'], sep='\t', file=fout)
    fout.close()

if __name__ == '__main__':

    get_gene_name()
