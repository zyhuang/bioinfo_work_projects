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


def get_rsid_annovar(gnomad_name, refgen_name, out_list_name, region_str=None):
    '''
    If region_str is None, query the whole table.
    '''

    if region_str != None:
        query_chrom = region_str.split(':')[0]
    else:
        query_chrom = None

    refgen_fa = fasta.Fasta(refgen_name)
    table = annovartable.AnnovarTable(gnomad_name, query_chrom=query_chrom)

    var_dict = {}
    col_pop = {}
    for i,x in enumerate(table.header):
        if i < 5:
            continue
        col_pop[i] = x.split('_')[-1]

    print('reading ' + gnomad_name, file=sys.stderr)
    for line in table.query_region(region_str):
        data = line.rstrip().split('\t')
        chrom, pstart, pend, ref, alt = data[:5]
        var = variant.Variant()
        try:
            var.parse(chrom, pstart, ref, alt, refgen_fa, left_norm=True)
        except ValueError as e:
            print('*WARNING* ' + str(e), file=sys.stderr)
            continue
        varkey = str(var)
        popaf_dict = {}
        for i,x in enumerate(data):
            if i < 5:
                continue
            popaf_dict[col_pop[i]] = float(x) if x != '.' else -1
        var_dict[varkey] = popaf_dict
    table.close()

    fout = open(out_list_name, 'w')
    print('writing ' + out_list_name, file=sys.stderr)
    for varkey in sorted(var_dict, key=lambda x: variant.sort_varkey(x)):
        print(varkey, json.dumps(var_dict[varkey], sort_keys=True),
              sep='\t', file=fout)
        var_dict.pop(varkey)
    fout.close()



if __name__ == '__main__':

    gnomad_name, refgen_name, out_list_name, region_str = sys.argv[1:5]
    get_rsid_annovar(gnomad_name, refgen_name, out_list_name,
                     region_str=region_str)
