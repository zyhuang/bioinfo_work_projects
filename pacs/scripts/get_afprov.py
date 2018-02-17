import sys
import json

def get_afprov(binkey, min_maf=0.01, min_ns=100):

    af_dir = ('/stornext/snfs2/1000GENOMES/zhuoyih/projects/221E/work-20161007-caller/work/nipt.214k/v1.1/analysis1/af_pdf_breve/data')

    # load varkey
    in_prov_list = 'prov.list'
    out_af_list = 'af_prov_list/{}.af.list'.format(binkey)

    prov_set = set()
    for line in open(in_prov_list):
        prov_set.add(line.rstrip())

    afprov_dict = {}
    in_af_list = '{}/{}/{}.af.list'.format(af_dir, 'All', binkey)
    print('reading {}'.format(in_af_list), file=sys.stderr)
    for line in open(in_af_list):
        varkey, vardata = line.rstrip().split('\t')
        vardata = json.loads(vardata)

        # filter sites with insufficient samples at global level
        if vardata['var_union']['0']['ns'] < min_ns:
            continue

        # filter sites with global small MAF
        af = vardata['af_pdf']['mean']['af']
        if af < min_maf or 1-af < min_maf:
            continue

        afprov_dict[varkey] = {}
        for prov in prov_set:
            if prov != 'All':
                afprov_dict[varkey][prov] = {
                    'af_mle': -1, 'af_max': -1, 'af_min': -1, 'af_raw': -1,
                    'ns_all': -1, 'ns_ra': -1, 'ns_rr': -1, 'ns_aa': -1,
                }
            else:
                afprov_dict[varkey][prov]= {
                    'af_mle': vardata['af_pdf']['mean']['af'],
                    'af_min': vardata['af_pdf']['ci_mean_95%'][0],
                    'af_max': vardata['af_pdf']['ci_mean_95%'][1],
                    'af_raw': vardata['var_union']['0']['af'],
                    'ns_all': vardata['var_union']['0']['ns'],
                    'ns_rr': vardata['var_union']['0']['rr'],
                    'ns_ra': vardata['var_union']['0']['ra'],
                    'ns_aa': vardata['var_union']['0']['aa'],
                }


    for prov in sorted(prov_set):
        if prov == 'All':
            continue
        in_afprov_list = '{}/{}/{}.af.list'.format(af_dir, prov, binkey)
        print('reading {}'.format(in_afprov_list), file=sys.stderr)
        for line in open(in_afprov_list):
            varkey, vardata = line.rstrip().split('\t')
            vardata = json.loads(vardata)

            # filter sites whose global MAF is small
            if varkey not in afprov_dict:
                continue

            # fitler sites with insufficient samples at province level
            if vardata['var_union']['0']['ns'] < min_ns:
                continue

            afprov_dict[varkey][prov]= {
                'af_mle': vardata['af_pdf']['mean']['af'],
                'af_min': vardata['af_pdf']['ci_mean_95%'][0],
                'af_max': vardata['af_pdf']['ci_mean_95%'][1],
                'af_raw': vardata['var_union']['0']['af'],
                'ns_all': vardata['var_union']['0']['ns'],
                'ns_rr': vardata['var_union']['0']['rr'],
                'ns_ra': vardata['var_union']['0']['ra'],
                'ns_aa': vardata['var_union']['0']['aa'],
            }

    fout = open(out_af_list, 'w')
    print('writing {}'.format(out_af_list), file=sys.stderr)
    for varkey in sorted(afprov_dict):
        print(varkey, json.dumps(afprov_dict[varkey], sort_keys=True),
              sep='\t', file=fout)
    fout.close()


binkey = sys.argv[1]
get_afprov(binkey)
