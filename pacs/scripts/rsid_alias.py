'''This program extracts pdist line of variants provided in rsid.list.

The output is in af_pdist_rsid_list dir.
'''

import sys
import json


def rsid_alias():

    for line in open('rsid.list'):
        rsid, qvarkey, others = line.rstrip().split('\t')
        chrom, pos = qvarkey.split(':')[:2]
        binkey = '{}.{}'.format(chrom, pos[:3])
        out_list = 'af_pdist_rsid_list/{}.pdist.list'.format(rsid)
        print('writing '+ out_list, file=sys.stderr)
        fout = open(out_list, 'w')
        for line in open('af_pdist_list/{}.pdist.list'.format(binkey)):
            varkey, vardata = line.rstrip().split('\t')
            if varkey == qvarkey:
                print(line.rstrip(), file=fout)
                break
        fout.close()


rsid_alias()
