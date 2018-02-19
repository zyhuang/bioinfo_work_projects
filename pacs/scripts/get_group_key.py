'''This program surveys the gene keys in pdist lists

skip intergenic and up/down-stream regions when examing gene/ncRNA.
exonic, splicing and intronic, UTR3/5 are considered for a gene/ncRNA to be
surveyed.

The output saves the name og gene/ncRNA in each binkey.

The "_all" version merges all gene/ncRNA from all bins and output names of bins
and the type (ncrna/gene) for each gene/ncRNA.
'''

import sys
import json

#         exonic:xxx, splicing, intronic, UTR5/UTR3,
#         ncRNA_exonic, ncRNA_splicing, ncRNA_intronic,
#         upstream/downstream, intergenic


def get_gene_binkey(binkey):

    in_list_name = '../af_pdist_list/{}.pdist.list'.format(binkey)
    out_list_name = 'temp/{}.gene.list'.format(binkey)

    gene_dict = {'gene':set(), 'ncrna':set()}
    print('reading ' + in_list_name, file=sys.stderr)
    for line in open(in_list_name):
        varkey, vardata = line.rstrip().split('\t')
        vardata = json.loads(vardata)
        for func in vardata['group']['func']:
            if 'x3b' in func:
                print('*WARNING*: found \x3b in {} of {}'.format(varkey, in_list_name))
            gene, anno = func.split(':')
            # skip intergenic and up/down-stream (including exonic, splicing, intronic, UTR
            if 'upstream/downstream' in anno or 'intergenic' in anno:
                continue
            if anno.startswith('ncRNA_'):
                gene_dict['ncrna'].add(gene)
            else:
                gene_dict['gene'].add(gene)

    gene_dict['gene'] = list(gene_dict['gene'])
    gene_dict['ncrna'] = list(gene_dict['ncrna'])
    fout = open(out_list_name, 'w')
    print('writing ' + out_list_name, file=sys.stderr)
    print(json.dumps(gene_dict, sort_keys=True), file=fout)
    fout.close()


def get_gene_binkey_all():

    gene_dict = {'gene': {}, 'ncrna': {}}
    for binkey in open('../binkey.list'):
        binkey = binkey.rstrip()
        for line in open('temp/{}.gene.list'.format(binkey)):
            data = json.loads(line.rstrip())
            for gtype in ['gene', 'ncrna']:
                for gene in data[gtype]:
                    if gene not in gene_dict[gtype]:
                        gene_dict[gtype][gene] = set()
                    gene_dict[gtype][gene].add(binkey)

    fout = open('gene.binkey.list', 'w')
    for gtype in ['gene', 'ncrna']:
        for gene in sorted(gene_dict[gtype]):
            # skip gene classified as both ncRNA and gene
            if gtype == 'gene' and gene in gene_dict['ncrna']:
                continue
            print(gene, gtype, ','.join(sorted(list(gene_dict[gtype][gene]))),
                  sep='\t', file=fout)
    fout.close()


if __name__ == '__main__':

    if len(sys.argv) == 2:
        binkey = sys.argv[1]
        get_gene_binkey(binkey)
    else:
        get_gene_binkey_all()
