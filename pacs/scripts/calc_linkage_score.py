'''This program calculate the clustering linkage vector for each line of input
based on the AF pdist.
'''

import sys
import json
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet


def prov_pair_to_prov_count(prov_pair):
    '''Get a list of unique provinces.'''

    prov_count = {}
    for pair in prov_pair:
        p1, p2 = pair.split(':')
        for p in [p1, p2]:
            if p not in prov_count:
                prov_count[p] = 0
            prov_count[p] += 1

    return prov_count


def prov_count_to_prov_pair(prov_count):
    '''create a unique pair of provinces.'''

    prov_list = sorted(prov_count.keys())
    prov_pair = []
    for i,p1 in enumerate(prov_list):
        for j,p2 in enumerate(prov_list):
            if j <= i:
                continue
            prov_pair.append('{}:{}'.format(p1,p2))

    return prov_pair


def get_pairs(pdist_dict):
    '''Get a list of province with valid pair-wise distance.

    The output list of N provinces matches with (N-1)N/2 pair-distance.
    '''

    exclude_prov = {'All', 'Unknown', 'Korea', 'Hainan'}
    prov_pair = []
    for pair in pdist_dict:
        p1, p2 = pair.split(':')
        if p1 in exclude_prov or p2 in exclude_prov:
            continue
        prov_pair.append(pair)

    prov_count = {}
    while True:
        prov_count = prov_pair_to_prov_count(prov_pair)
        n = len(prov_count)
        m = len(prov_pair)
        # print(n, n*(n-1)/2, m)
        if m == (n-1)*n/2:
            break

        # drop prov with min pair-wise distance
        prov_min = sorted(prov_count.items(), key=lambda x:x[1])[0][0]

        # print('pop prov_min {}: {}'.format(prov_min, prov_count[prov_min]))
        prov_count.pop(prov_min)

        # drop prov in prov_pair
        new_prov_pair = []
        for pair in prov_pair:
            p1, p2 = pair.split(':')
            if p1 == prov_min or p2 == prov_min:
                continue
            new_prov_pair.append(pair)
        prov_pair = new_prov_pair


    prov_list = sorted(prov_count.keys())
    prov_pair = prov_count_to_prov_pair(prov_count)
    # print(prov_list)

    return prov_pair, prov_list


def calc_linkage_score(pdist_list_name, out_linkage_name):

    print('reading ' + pdist_list_name)
    print('writing ' + out_linkage_name)
    fout = open(out_linkage_name, 'w')
    for line in open(pdist_list_name):
        key, data = line.rstrip().split('\t')
        data = json.loads(data)
        prov_pair, prov_list = get_pairs(data['pdist'])
        if len(prov_list) <= 2:
            continue

        prov_dist = []
        for pair in prov_pair:
            dist = np.sqrt(data['pdist'][pair]['df2']/data['pdist'][pair]['n'])
            prov_dist.append(dist)

        linkage_out = {}
        for method in ['single', 'complete', 'average', 'weighted', 'centroid',
                       'median', 'ward']:
            linkage_matrix = linkage(prov_dist, method, optimal_ordering=False)
            cophenet_coeff, cophenet_dist = cophenet(linkage_matrix, prov_dist)
            linkage_matrix = linkage_matrix.tolist()
            linkage_out[method] = {
                'cophenet_coeff': cophenet_coeff,
                'linkage_coeff': [
                    linkage_matrix[-1][-2],
                    linkage_matrix[-1][-1],   # last level height and #node
                    linkage_matrix[-2][-2],
                    linkage_matrix[-2][-1],   # 2nd-last height and #node
                ],
                # 'linkage_matrix': linkage_matrix,
            }

        out_dict = {
            'prov_list': prov_list,
            'prov_dist': prov_dist,
            'linkage_out': linkage_out,
        }

        print(key, json.dumps(out_dict, sort_keys=True), sep='\t', file=fout)

    fout.close()


if __name__ == '__main__':

    pdist_list_name, out_linkage_name = sys.argv[1:3]
    calc_linkage_score(pdist_list_name, out_linkage_name)
