import sys
from subprocess import PIPE, Popen


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


def has_snp(ref, alt_list):

    if len(ref) > 1:
        return False

    for alt in alt_list:
        if len(alt) == 1:
            return True

    return False


def load_gt_stat(vcf_data, sample_header_hash):

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

    af_stat = {}
    for pop,gt_dict in gt_stat.items():
        ns = 0
        ac = 0
        for gt,count in gt_dict.items():
            ns += count
            if gt == '{}|{}'.format(alt_index+1, alt_index+1):
                ac += 2*count
            elif (gt.startswith('{}|'.format(alt_index+1)) or
                  gt.endswith('|{}'.format(alt_index+1))):
                ac += count
        af_stat[pop] = ac/ns/2

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



def subpop_af(vcf_name, panel_name):

    sample_header_hash = {}
    proc = Popen('zcat {}'.format(vcf_name), shell=True, stdout=PIPE,
                 universal_newlines=True)
    for line in proc.stdout:
        if line.startswith('##'):
            continue
        if line.startswith('#'):
            sample_header_hash = load_sample_header(line.rstrip(), panel_name)
            continue

        vcf_data = line.rstrip().split('\t')
        chrom, pos, rsid, ref, alts = vcf_data[:5]
        alt_list = alts.split(',')

        # skip non-SNP
        if not has_snp(ref, alt_list):
            continue

        # load stats of GT by populations
        gt_stat = load_gt_stat(vcf_data, sample_header_hash)

        # load info data
        info_af = parse_info_af(vcf_data)


        for i, alt in enumerate(alt_list):
            if len(alt) == 1:
                varkey = '{}:{}:{}:{}'.format(chrom, pos.zfill(9), ref, alt)
                af_stat = calc_alt_freq(gt_stat, i)
                # checking
                for p in ['AFR', 'AMR', 'EUR', 'EAS', 'SAS', 'ALL']:
                    q = p+'_AF' if p != 'ALL' else 'AF'
                    if abs(af_stat[p] - info_af[q][i]) > 1e-3:
                        print(vcf_data[:8])
                        print(p, q, af_stat[p], info_af[q][i])
                        print(varkey)
                        print(af_stat)
                        print(gt_stat)
                        print()



    proc.communicate()


def test():

    subpop_af(
        'test/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',
        'test/integrated_call_samples_v3.20130502.ALL.panel')


test()
