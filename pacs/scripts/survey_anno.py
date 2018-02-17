import sys
import json


def survey_binkey(binkey, annokey_stat):

    af_list = '../anno_list/{}.anno.list'.format(binkey)
    print('reading {}'.format(af_list), file=sys.stderr)
    for line in open(af_list):
        varkey, vardata = line.rstrip().split('\t')
        vardata = json.loads(vardata)
        annokey = [
            vardata['ExonicFunc.refGene'],
            vardata['ExonicFunc.wgEncodeGencodeBasicV19'],
            vardata['Func.refGene'],
            vardata['Func.wgEncodeGencodeBasicV19'],
        ]
        annokey = '|'.join(annokey)
        if annokey not in annokey_stat:
            annokey_stat[annokey] = 0
        annokey_stat[annokey] += 1



if __name__ == '__main__':


    annokey_stat = {}
    nbin = 0
    for binkey in open('../binkey.list'):
        nbin += 1
        binkey = binkey.rstrip()
        survey_binkey(binkey, annokey_stat)
        # if nbin == 100:
        #     break

    fout_name = 'survey.out'
    print('writing {}'.format(fout_name), file=sys.stderr)
    fout = open(fout_name, 'w')
    for annokey, count in sorted(annokey_stat.items()):
        print(count, annokey, sep='\t', file=fout)
    fout.close()
