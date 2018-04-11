'''This progrom parses the dbsnp json data and save the results in a line-based
key-value json list file.

The JSON data is downloaded from the following FTP
ftp://ftp.ncbi.nlm.nih.gov/snp/.redesign/latest_release/JSON/

The key is of format '{chrom}:{pos09}:{ref}:{alt}'.
The variants are left-normalized, ref or alt may be * for indels.

The extracted information about each article are...

'''

# list can be empty



import sys
import json
import gzip


def get_value(data, keypath_str, default=None):

    try:
        key_list = keypath_str.split('.')
        tmp_data = data
        for key in key_list:
            tmp_data = tmp_data[key]
        value = tmp_data
    except:
        value = default

    return value


def get_meta(data):

    meta = {}
    meta['rsid'] = get_value(data, 'refsnp_id')
    meta['type'] = get_value(data, 'primary_snapshot_data.variant_type')
    meta['date_create'] = get_value(data, 'create_date')
    meta['date_update'] = get_value(data, 'last_update_date')
    meta['build_id'] = get_value(data, 'last_update_build_id')
    meta['citations'] = get_value(data, 'citations')

    return meta


def get_alleles(data, assembly_keyword='GRCh37'):
    '''Get chromosome variant (NC) keys and their corresponding HGVS notation.'''

    assembly_keyword += '.'
    allele_list = []

    chrom_map = {
        '23': 'X', '24': 'Y', '12920': 'MT',
    }
    for placement in data['primary_snapshot_data']['placements_with_allele']:
        if not placement['seq_id'].startswith('NC_'):
            continue
        seq_id = placement['seq_id']
        chrom = seq_id.split('.')[0].split('_')[1].lstrip('0')
        if chrom in chrom_map:
            chrom = chrom_map[chrom]
        for assembly in placement['placement_annot']['seq_id_traits_by_assembly']:
            if not (assembly['assembly_name'].startswith(assembly_keyword) and
                    assembly['is_chromosome']):
                continue
            for allele in placement['alleles']:
                hgvs = allele['hgvs']
                pos = allele['allele']['spdi']['position'] + 1
                ref = allele['allele']['spdi']['deleted_sequence']
                alt = allele['allele']['spdi']['inserted_sequence']
                tmp = {
                    'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
                    'hgvs': hgvs,
                }
                allele_list.append(tmp)

    return allele_list


def get_annotations(data):

    anno_list = []

    # collect per allele annotation data
    # each allele has annotation, clinical, frequency and submissions
    for annotation in data['primary_snapshot_data']['allele_annotations']:

        allele = {}
        allele['submission'] = len(annotation['submissions'])

        allele['frequency'] = []
        for freq in annotation['frequency']:
            tmp = {
                'study': '{}.v{}'.format(freq['study_name'],
                                         freq['study_version']),
                'ac': freq['allele_count'],
                'an': freq['total_count'],
                'af': freq['allele_count']/freq['total_count'],
            }
            allele['frequency'].append(tmp)

        allele['clinical'] = []
        for clin in annotation['clinical']:
            tmp = {
                'citations': clin['citations'],
                'significance': clin['clinical_significances'],
                'collection': clin['collection_method'],
                'disease': clin['disease_names'],
                'review': clin['review_status'],
                'date_create': clin['create_date'],
                'date_update': clin['update_date'],
                'date_evalue': get_value(clin, 'last_evaluated_date'),
            }
            allele['clinical'].append(tmp)

        allele['annotation'] = {'genes': []}
        for anno in annotation['assembly_annotation']:
            # gene also means at DNA level here (c.f. rna/protein level below)
            for gene in anno['genes']:
                tmp_gene = {
                    'anno_release': anno['annotation_release'],
                    'locus': gene['locus'],
                    'name': gene['name'],
                    'is_pseudo': gene['is_pseudo'],
                    'is_minus': bool(gene['orientation'] == 'minus'),
                    'so_term': [x['name'] for x in gene['sequence_ontology']],
                    'rnas': [],
                }
                for rna in gene['rnas']:
                    tmp_rna = {
                        'so_term': [x['name'] for x in rna['sequence_ontology']],
                        'variant': {},
                        'protein': {},
                    }
                    if 'transcript_change' in rna:
                        tmp_rna['variant'] = {
                            'seq': rna['transcript_change']['seq_id'],
                            'pos': rna['transcript_change']['position'],
                            'ref': rna['transcript_change']['deleted_sequence'],
                            'alt': rna['transcript_change']['inserted_sequence'],
                        }

                    if 'protein' in rna:
                        if 'spdi' in rna['protein']['variant']:
                            x = rna['protein']['variant']['spdi']
                        else:
                            x = rna['protein']['variant']['frameshift']
                        tmp_protein = {
                            'so_term': [x['name'] for x in rna['protein']['sequence_ontology']],
                            'variant': {
                                'seq': x.get('seq_id'),
                                'pos': x.get('position'),
                                'ref': x.get('deleted_sequence'),
                                'alt': x.get('inserted_sequence'),
                            },
                        }
                        tmp_rna['protein'] = tmp_protein
                    tmp_gene['rnas'].append(tmp_rna)
                allele['annotation']['genes'].append(tmp_gene)

        anno_list.append(allele)

    return anno_list


def get_snp_info(data):

    # skip variants merged to another one
    if 'merge_snapshot_data' in data:
        return {}

    meta = get_meta(data)
    allele_list = get_alleles(data)

    # skip variant only present in GRCh38
    if len(allele_list) == 0:
        return {}

    anno_list = get_annotations(data)
    if len(anno_list) != len(allele_list):
        return {}

    out = {
        'meta': meta,
        'alleles': [],
    }

    for allele,anno in zip(allele_list, anno_list):
        if allele['alt'] == allele['ref']:
            continue
        tmp = {}
        tmp['variant'] = allele
        for k in anno:
            tmp[k] = anno[k]
        out['alleles'].append(tmp)

    return out


def make_varkey(chrom, pos, ref, alt):

    if alt == '':
        alt = '*'
    if ref == '':
        ref = '*'

    if len(ref) > 1 and len(alt) > 1:
        # left normalization
        tail = ''
        for i in range(min(len(ref),len(alt))-1):
            if ref[-i-1] == alt[-i-1]:
                tail += ref[-i-1]
            else:
                break
        if tail:
            tail = tail[::-1]
            ref = ref.rsplit(tail, 1)[0]
            alt = alt.rsplit(tail, 1)[0]

        head = ''
        for i in range(min(len(ref),len(alt))-1):
            if ref[i] == alt[i]:
                head += ref[i]
            else:
                break
        if head:
            pos += len(head)
            ref = ref.split(head, 1)[-1]
            alt = alt.split(head, 1)[-1]

    varkey ='{}:{}:{}:{}'.format(chrom, pos, ref, alt)

    return varkey


def parse_dbsnp_json(json_name, out_json_name):

    print(':: reading ' + json_name, file=sys.stderr, flush=True)
    if json_name.endswith('.gz'):
        fin = gzip.open(json_name)
    else:
        fin = open(json_name)

    print(':: writing ' + out_json_name, file=sys.stderr, flush=True)
    fout = open(out_json_name, 'w')
    for line in fin:
        if json_name.endswith('.gz'):
            line = line.decode('utf-8')
        line = line.rstrip()
        data = json.loads(line)

        try:
            out= get_snp_info(data)
        except Exception as e:
            print('*ERROR* can not process' + data['refsnp_id'],
                  file=sys.stderr)
            raise

        if not out:
            continue

        for allele in out['alleles']:
            # add meta data shared between alleles at the same locus
            allele['meta'] = out['meta']

            chrom = allele['variant']['chrom']
            pos = str(allele['variant']['pos']).zfill(9)
            ref = allele['variant']['ref']
            alt = allele['variant']['alt']
            varkey = make_varkey(chrom, pos, ref, alt)
            print(varkey, json.dumps(allele, sort_keys=True),
                  sep='\t', file=fout)

            # a,b,c,d = chrom, pos, ref, alt
            # if len(c) == 1 or len(d) == 1:
            #     continue
            # rr, aa = varkey.split(':')[2:]
            # if len(rr) > 1 and len(aa) > 1:
            # # if varkey.split(':')[1] != b:
                # print(varkey, a,b,c,d, out['meta']['rsid'])

    fin.close()
    fout.close()


if __name__ == '__main__':

    json_name, out_json_name = sys.argv[1:3]
    # json_name, out_json_name = 'JSON/refsnp-chr22.json.gz', 'foo.json'
    parse_dbsnp_json(json_name, out_json_name)
