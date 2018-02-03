'''This module calculate AF in each population of 1000G in a given input VCF.

'''
import sys
from os.path import dirname, abspath
root_dir = dirname(dirname(abspath(__file__)))
sys.path.append('{}/sequentia'.format(root_dir))


import variant
import fasta
import annovartable
import json


def get_rsid_annovar(avsnp_name, refgen_name, out_list_name, region_str=None):
    '''
    If region_str is None, query the whole table.
    '''

    refgen_fa = fasta.Fasta(refgen_name)
    table = annovartable.AnnovarTable(avsnp_name)

    var_dict = {}
    print('reading ' + avsnp_name, file=sys.stderr)
    for line in table.query_region(region_str):
        chrom, pstart, pend, ref, alt, rsid = line.split('\t')
        var = variant.Variant()
        try:
            var.parse(chrom, pstart, ref, alt, refgen_fa, left_norm=True)
        except ValueError as e:
            print('*WARNING* ' + str(e), file=sys.stderr)
            continue
        varkey = str(var)
        var_dict[varkey] = rsid

    fout = open(out_list_name, 'w')
    print('writing ' + out_list_name, file=sys.stderr)
    for varkey in sorted(var_dict, key=lambda x: variant.sort_varkey(x)):
        print(varkey, var_dict[varkey], sep='\t', file=fout)
    fout.close()


if __name__ == '__main__':

    avsnp_name, refgen_name, out_list_name, region_str = sys.argv[1:5]
    get_rsid_annovar(avsnp_name, refgen_name, out_list_name,
                     region_str=region_str)
