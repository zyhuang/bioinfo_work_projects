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


def calc_afpdist_bygene(binkey):







if __name__ == '__main__':

    group, binkey = sys.argv[1:3]
    if group == 'region':
        calc_afpdist_byregion(binkey)
