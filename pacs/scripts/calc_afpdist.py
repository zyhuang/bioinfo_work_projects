''' This program calculates pair-wise distance of AF between two provinces for
any group of variants.

For each variant, assign a region group key, and functional group keys,
calculate pair-wise distance of AF between any unique pair of different
provinces.
'''

import sys
import json


def get_reg_key(varkey):
    ''' Only 10kb group. 100kb and 1Mb can be further grouped based on this.'''

    chrom, pos = varkey.split(':')[:2]
    return '{}.{}'.format(chrom, pos[:5])


def get_anno_keys(vardata):
    ''' Only consider refGene here.

    Input func key-value: ExonicFunc.refGene, Func.refGene, Gene.refGene

    Output func keys in following set:
        exonic:xxx, splicing, ncRNA_exonic, ncRNA_splicing, intronic,
        ncRNA_intronic, UTR5/UTR3, upstream/downstream, intergenic
    '''

    rename_rule = {
        'exonic\x3bsplicing': ['exonic', 'splicing'],
        'ncRNA_exonic\x3bsplicing': ['ncRNA_exonic', 'ncRNA_splicing'],
        'UTR5': ['UTR5/UTR3'],
        'UTR3': ['UTR5/UTR3'],
        'UTR5\\x3bUTR3': ['UTR5/UTR3'],
        'upstream': ['upstream/downstream'],
        'downstream': ['upstream/downstream'],
        'upstream\\x3bdownstream': ['upstream/downstream'],
    }

    if vardata['Func.refGene'] in rename_rule:
        func = rename_rule[vardata['Func.refGene']]
    else:
        func = [vardata['Func.refGene']]

    if func[0] == 'exonic':
        func[0] += '.' + vardata['ExonicFunc.refGene']

    genes = vardata['Gene.refGene'].replace('\\x3b',',').split(',')

    anno_keys = []
    for f in func:
        for g in genes:
            anno_keys.append(g + ':' + f)

    return anno_keys


def calc_pdist(vardata):
    '''Calculate pair-wise AF distance between all provinces.
    only for binom(n,2) half off-diagnal values
    '''

    pdist_dict = {}
    prov_list = sorted(vardata)
    for i,p1 in enumerate(prov_list):
        for j,p2 in enumerate(prov_list):
            if j <= i:
                continue
            # if insufficient sampling
            if vardata[p1]['ns_all'] <= 0 or vardata[p2]['ns_all'] <= 0:
                continue
            key = '{}:{}'.format(p1, p2)
            df = vardata[p1]['af_mle'] - vardata[p2]['af_mle']
            pdist_dict[key] = {
                'n': 1 , 'df': df, 'df2': df*df,
            }

    return pdist_dict


def calc_pdist_af(binkey):

    af_list = 'af_prov_list/{}.af.list'.format(binkey)
    anno_list = 'anno_list/{}.anno.list'.format(binkey)
    out_pdist_list = 'af_pdist_list/{}.pdist.list'.format(binkey)

    group_dict = {}
    print('reading ' + anno_list)
    for line in open(anno_list):
        varkey, vardata = line.rstrip().split('\t')
        vardata = json.loads(vardata)
        anno_keys = get_anno_keys(vardata)
        region_key = get_reg_key(varkey)
        group_dict[varkey] = {
            'group': {
                'func': anno_keys,
                'region': region_key,
            },
            'anno': {
                'rsid': vardata['avsnp147'],
                'cadd': vardata['CADD_Phred'],
            },
        }

    fout = open(out_pdist_list, 'w')
    print('reading ' + af_list)
    print('writing ' + out_pdist_list)
    for line in open(af_list):
        varkey, vardata = line.rstrip().split('\t')
        vardata = json.loads(vardata)
        pdist_dict = calc_pdist(vardata)
        group_dict[varkey]['pdist'] = pdist_dict
        print(varkey, json.dumps(group_dict[varkey], sort_keys=True),
              sep='\t', file=fout)
    fout.close()


if __name__ == '__main__':

    binkey = sys.argv[1]
    calc_pdist_af(binkey)
