'''This module calculate AF in each population of 1000G in a given input VCF.

'''
import sys
from os.path import dirname, abspath
root_dir = dirname(dirname(abspath(__file__)))
sys.path.append('{}/sequentia'.format(root_dir))

import variant
import fasta
from subprocess import PIPE, Popen
import json



def load_sample_pop(panel_name):

    sample_pop_hash = {}
    for line in open(panel_name):
        if line.startswith('sample'):
            continue
        sample, pop1, pop2, gender = line.rstrip().split('\t')
        sample_pop_hash[sample] = {
            'subpop': pop1, 'suppop': pop2,
        }
    return sample_pop_hash


def load_sample_header(vcf_line, panel_name):

    sample_pop_hash = load_sample_pop(panel_name)
    sample_header_hash = {}
    for i,s in enumerate(vcf_line.rstrip().split('\t')):
        if i < 9:
            continue
        sample_header_hash[i] = sample_pop_hash[s]
    return sample_header_hash


def load_gt_stat(vcf_data, sample_header_hash):
    '''
    return gt_stat = {'ALL': {"0/1": cnt, "0/0": cnt, ...}, "EAS": {...}, ...}
    '''

    gt_stat = {}
    for i, gt in enumerate(vcf_data):
        if i < 9:
            continue
        subpop = sample_header_hash[i]['subpop']
        suppop = sample_header_hash[i]['suppop']

        for pop in [subpop, suppop, 'ALL']:
            if pop not in gt_stat:
                gt_stat[pop] = {}
            if gt not in gt_stat[pop]:
                gt_stat[pop][gt] = 0
            gt_stat[pop][gt] += 1

    return gt_stat


def calc_alt_freq(gt_stat, alt_index):
    '''
    alt_index is 0-based in the alt list, alt = 0 is the 1 of "1/2" in genotype
    return af_stat = {'ALL': {'ns': ns, 'af': af, 'ns1': ns_male, 'af1': af_male,
                      'ns2': ns_female, 'af2': af_female, ...}, 'EAS': {...}, ...}
    '''

    af_stat = {}
    for pop, gt_dict in gt_stat.items():
        ns1 = 0   # 1 chromosome
        ac1 = 0
        ns2 = 0   # 2 chromosomes
        ac2 = 0
        for gt, count in gt_dict.items():
            # two chromosomes
            if '|' in gt:  # | for phased genotype in 1000G
                ns2 += count
                if gt == '{}|{}'.format(alt_index+1, alt_index+1):
                    ac2 += count*2
                elif (gt.startswith('{}|'.format(alt_index+1)) or
                      gt.endswith('|{}'.format(alt_index+1))):
                    ac2 += count
            # one chromosome
            else:
                 ns1 += count
                 if gt == '{}'.format(alt_index+1):
                     ac1 += count

        ns = ns1 + ns2
        af1 = ac1 / ns1 if ns1 else None
        af2 = ac2 / (2*ns2) if ns2 else None
        af = (ac1 + ac2) / (ns1 + 2*ns2) if (ns1 + 2*ns2) else None

        af_stat[pop] = {
            'af': af, 'ns': ns,
            'af1': af1, 'ns1': ns1,
            'af2': af2, 'ns2': ns2,
        }

    return af_stat


def parse_info_af(vcf_data):

    info_af = {}
    for x in vcf_data[7].split(';'):
        if '=' not in x:
            continue
        key, values = x.split('=')
        if key.endswith('AF'):
            value_list = [float(v) for v in values.split(',')]
            info_af[key] = value_list

    return info_af



def subpop_af(vcf_name, panel_name, refgen_name, out_list_name,
              region_str=None):
    '''

    Note: about RsID.
    '''

    refgen_fa = fasta.Fasta(refgen_name)
    fout = open(out_list_name, 'w')
    varkey_set = set()

    if region_str:
        cmd = 'tabix -h {} {}'.format(vcf_name, region_str)
    else:
        cmd = 'zcat {}'.format(vcf_name)

    proc = Popen(cmd, shell=True, stdout=PIPE, universal_newlines=True)

    print('> ' + cmd, file=sys.stderr)
    print('writing ' + out_list_name, file=sys.stderr)

    sample_header_hash = {}
    for line in proc.stdout:
        if line.startswith('##'):
            continue
        if line.startswith('#'):
            sample_header_hash = load_sample_header(line.rstrip(), panel_name)
            continue

        vcf_data = line.rstrip().split('\t')

        chrom, pos, rsids, ref, alts = vcf_data[:5]
        alt_list = alts.split(',')
        # rsid_list = rsids.split(';')
        # if len(alt_list) != len(rsid_list):
        #     rsid_list = [rsids for alt in alt_list]

        # ----- debug -----
        # if not vcf_data[4].startswith('<'):
        #     continue
        # if len(alt_list) == 1:
        #     continue

        # for validation only
        info_af = parse_info_af(vcf_data)

        # loop through multiallelic variants
        gt_stat = None
        for i, alt in enumerate(alt_list):
            var = variant.Variant()
            try:
                var.parse(chrom, pos, ref, alt, refgen_fa, left_norm=True)
            except ValueError as e:
                print('*WARNING* ' + str(e), file=sys.stderr)
                continue
            # if var.var_type() != 'SNV':
            #     continue

            varkey = str(var)
            if not gt_stat:
                gt_stat = load_gt_stat(vcf_data, sample_header_hash)
            af_stat = calc_alt_freq(gt_stat, i)

            # to prevent duplicated variants (appearing on chr8, 12, 14 and X)
            if varkey not in varkey_set:
                print(varkey, json.dumps(af_stat, sort_keys=True),
                      sep='\t', file=fout)
                varkey_set.add(varkey)

            # ----- debug -----
            # print(vcf_data[:8])
            # print(varkey)
            # print(af_stat)
            # print(gt_stat)
            # print()

            # debugging print, check if the calculated AF the same as those in INFO
            for p in ['AFR', 'AMR', 'EUR', 'EAS', 'SAS', 'ALL']:
                q = p+'_AF' if p != 'ALL' else 'AF'
                if abs(af_stat[p]['af'] - info_af[q][i]) > 1e-3:
                    print(vcf_data[:8], file=sys.stderr)
                    print(p, q, af_stat[p], info_af[q][i], file=sys.stderr)
                    print(varkey, file=sys.stderr)
                    print(af_stat, file=sys.stderr)
                    print(gt_stat, file=sys.stderr)
                    print(file=sys.stderr)

    proc.communicate()
    fout.close()


def test():

    subpop_af(
        'test/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',
        'test/integrated_call_samples_v3.20130502.ALL.panel',
        '/data1/work/niptcaller/bundle/refgen/bwa/human.g1k.v37.fa',
        'test.list', region_str='20:1000000-2000000')

# test()


if __name__ == '__main__':

    vcf_name, panel_name, refgen_name, out_list_name = sys.argv[1:5]
    if len(sys.argv) == 6:
        region_str = sys.argv[5]
    else:
        region_str = None

    subpop_af(vcf_name, panel_name, refgen_name, out_list_name,
              region_str=region_str)
