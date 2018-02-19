'''This program aggregates the afpdist stats by region/func group.

output each group afpdist
'''

import sys
import json


def calc_afpdist_byregion(binkey):

    afpdist_list = 'af_pdist_list/{}.pdist.list'.format(binkey)
    out_afpdist_list = 'af_pdist_region_list/{}.pdist.list'.format(binkey)

    region_pdist = {}
    print('reading ' + afpdist_list, file=sys.stderr)
    for line in open(afpdist_list):
        varkey, vardata = line.rstrip().split('\t')
        vardata = json.loads(vardata)
        region_10kb = vardata['group']['region']
        region_100kb = region_10kb[:-1]
        region_1mb = region_10kb[:-2]
        for reg in [region_10kb, region_100kb, region_1mb]:
            if reg not in region_pdist:
                region_pdist[reg] = {'pdist': {}}
            for pair in vardata['pdist']:
                if pair not in region_pdist[reg]['pdist']:
                    region_pdist[reg]['pdist'][pair] = {
                        'n': 0, 'df': 0, 'df2': 0,
                    }
                for x in ['n', 'df', 'df2']:
                    region_pdist[reg]['pdist'][pair][x] \
                        += vardata['pdist'][pair][x]

    print('writing ' + out_afpdist_list, file=sys.stderr)
    fout = open(out_afpdist_list, 'w')
    for reg in sorted(region_pdist):
        print(reg, json.dumps(region_pdist[reg], sort_keys=True),
              sep='\t', file=fout)
    fout.close()


def calc_afpdist_bygene(qgene):

    gene_list = 'gene.binkey.list'

    qtype = None
    qbinkeys = []
    print('reading ' + gene_list, file=sys.stderr)
    for line in open(gene_list):
        gene, gtype, binlist = line.rstrip().split('\t')
        if qgene != gene:
            continue
        qtype = gtype
        qbinkeys = binlist.split(',')

    # qgene may contain '/'
    out_gene_list = 'af_pdist_{}_list/{}.pdist.list'.format(qtype, qgene.replace('/','_'))

    pdist_dict = {'pdist': {}}
    for binkey in qbinkeys:
        af_list = 'af_pdist_list/{}.pdist.list'.format(binkey)
        print('reading ' + af_list, file=sys.stderr)
        for line in open(af_list):
            varkey, vardata = line.rstrip().split('\t')
            vardata = json.loads(vardata)
            for func in vardata['group']['func']:
                gene, anno = func.split(':')
                if gene != qgene:
                    continue
                if anno in ['upstream/downstream', 'intergenic']:
                    continue

                for pair in vardata['pdist']:
                    if pair not in pdist_dict['pdist']:
                        pdist_dict['pdist'][pair] = {
                            'n': 0, 'df': 0, 'df2': 0,
                        }
                    for x in ['n', 'df', 'df2']:
                        pdist_dict['pdist'][pair][x] \
                            += vardata['pdist'][pair][x]

                break   # each variant is counted once

    print('writing ' + out_gene_list, file=sys.stderr)
    fout = open(out_gene_list, 'w')
    print(qgene, json.dumps(pdist_dict, sort_keys=True),
          sep='\t', file=fout)
    fout.close()




if __name__ == '__main__':

    group = sys.argv[1]
    if group == 'region':
        binkey = sys.argv[2]
        calc_afpdist_byregion(binkey)
    elif group == 'gene':
        qgene = sys.argv[2]
        calc_afpdist_bygene(qgene)
