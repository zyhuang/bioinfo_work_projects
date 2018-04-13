import sys
import json
import os
import re


def make_bin_key(var_key):
    chrom, pos = var_key.split(':')[:2]
    bin_key = '{}.{}'.format(chrom, pos[:3])
    return bin_key


def sort_json(json_name):

    print('sorting ' + json_name, file=sys.stderr)
    sort_cmd = 'sort -t\: -gk2,2 {} -o {}'.format(json_name, json_name)
    print('> ' + sort_cmd, file=sys.stderr)
    os.system(sort_cmd)


def split_dbsnp_bin(json_name, out_dir):

    fout_dict = {}

    chrom = re.findall(r'refsnp-chr(.*?).json', json_name)[0]
    print('reading ' + json_name, file=sys.stderr)
    for line in open(json_name):
        line = line.rstrip()
        var_key, value = line.split('\t')
        bin_key = make_bin_key(var_key)
        out_json_name = '{}/{}.dbsnp.from.{}.json'.format(out_dir, bin_key, chrom)
        if out_json_name not in fout_dict:
            fout_dict[out_json_name] = open(out_json_name, 'w')
        print(line, file=fout_dict[out_json_name])

    for out_json_name in fout_dict:
        fout_dict[out_json_name].close()
        sort_json(out_json_name)


if __name__ == '__main__':

    json_name, out_dir = sys.argv[1:3]
    split_dbsnp_bin(json_name, out_dir)
