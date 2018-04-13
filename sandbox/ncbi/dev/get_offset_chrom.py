import os
import sys
import json
import requests


def get_chrom_offset(rsid):

    url = 'https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/' + rsid
    print('get ' + url)
    r = requests.get(url)
    data = json.loads(r.text)
    fout = open('foo.json', 'w')
    print(json.dumps(data, sort_keys=True, indent=2), file=fout)
    fout.close()

    allele_list = data['primary_snapshot_data']['placements_with_allele']
    hgvs_list = []
    for allele in allele_list:
        assembly = False
        for x in allele['placement_annot']['seq_id_traits_by_assembly']:
            if ((x['assembly_name'].startswith('GRCh37.') or
                 x['assembly_name'].startswith('GRCh38.'))
                 and x['is_chromosome'] and not x['is_patch']):
                assembly = x['assembly_name']
                break

        if assembly:
            for a in allele['alleles']:
                if (a['allele']['spdi']['deleted_sequence'] !=
                    a['allele']['spdi']['inserted_sequence']):
                    hgvs_list.append(assembly + '/' + a['hgvs'])

    # print(rsid, sorted(hgvs_list))
    return sorted(hgvs_list)


def loop():

    fout = open('chrom_offset.list', 'w')
    for fname in sorted(os.listdir('../note')):
        fname = '../note/' + fname
        for line in open(fname):
            key, value = line.rstrip().split('\t')
            value = json.loads(value)
            rsid = value['meta']['rsid']
            hgvs_list = get_chrom_offset(rsid)
            out = {
                'varkey': key,
                'hgvs': hgvs_list
            }
            print(rsid, out, sep='\t', file=fout)
    fout.close()

if __name__ == '__main__':

    loop()
