import sys
import json


def get_anno(binkey):

    anno_dir = ('/stornext/snfs2/1000GENOMES/zhuoyih/projects/221E/work-20161007-caller/work/nipt.214k/v1.1/anno/function/out')

    in_af_list = 'af_prov_list/{}.af.list'.format(binkey)
    in_anno_list = '{}/{}.anno.list'.format(anno_dir, binkey)
    out_anno_list = 'anno_list/{}.anno.list'.format(binkey)

    varkey_set = set()

    print('reading {} and {}'.format(in_af_list, in_anno_list), file=sys.stderr)
    print('writing {}'.format(out_anno_list), file=sys.stderr)

    fout = open(out_anno_list, 'w')

    for line in open(in_af_list):
        varkey, vardata = line.rstrip().split('\t')
        varkey_set.add(varkey)

    for line in open(in_anno_list):
        varkey, vardata = line.rstrip().split('\t')
        if varkey in varkey_set:
            print(line.rstrip(), file=fout)

    fout.close()


binkey = sys.argv[1]
get_anno(binkey)
