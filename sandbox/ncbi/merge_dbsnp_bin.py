'''Merge bins on GRCh37 from different chromosomes on GRCh38.

The downloaded chromosomal partition of dbsnp-150 json data is based on
GRCh38. Some variants are on different chromosome on assembly GRCh37.
Therefore all variants in the same bin on GRCh37 are collected, merged and
sorted using this script.

'''
import sys
import os


def merge(merge_list, out_dir):

    bin_dict = {}
    for line in open(merge_list):
        line = line.rstrip()
        key = line.split('/')[-1].rsplit('.', 4)[0]
        if key not in bin_dict:
            bin_dict[key] = []
        bin_dict[key].append(line)

    for key in bin_dict:
        file_list = ' '.join(bin_dict[key])
        file_out = '{}/{}.dbsnp.json'.format(out_dir, key)
        cmd = 'cat {} | sort -gk2,2 > {}'.format(file_list, file_out)
        print(':: ' + cmd)
        os.system(cmd)

merge('merge_dbsnp_bin.list', '../snp')
