import json
import sys
import numpy as np


def calc_pop_angular_distance(pop1, pop2, af_list_name, nfile_max=-1):

    af1_list = []
    af2_list = []

    nfile = 0
    nvar = 0
    for list_name in open(af_list_name):
        list_name = list_name.rstrip()

        nfile += 1
        for line in open(list_name):
            varkey, vardata = line.rstrip().split('\t')
            vardata = json.loads(vardata)

            ns1 = vardata[pop1]['ns_all']
            ns2 = vardata[pop2]['ns_all']

            if ns1 < 0 or ns2 < 0:
                continue

            af1 = vardata[pop1]['af_mle']
            af2 = vardata[pop2]['af_mle']
            af1_list.append(af1)
            af2_list.append(af2)
            nvar += 1

        if nfile == nfile_max:
            break

    if nvar:
        cov_matrix = np.cov(af1_list, af2_list)
        print(pop1, pop2, nvar,
              cov_matrix[0][0], cov_matrix[1][1], cov_matrix[0][1],
              cov_matrix[0][1]/np.sqrt(cov_matrix[0][0]*cov_matrix[1][1]),
              sep='\t')
    else:
        print(pop1, pop2, nvar, 0, 0, 0, 0, sep='\t')


def calc_pop_euclidean_distance(pop1, pop2, af_list_name, nfile_max=-1):

    df_list = []
    nfile = 0
    for list_name in open(af_list_name):
        list_name = list_name.rstrip()

        nfile += 1
        for line in open(list_name):
            varkey, vardata = line.rstrip().split('\t')
            vardata = json.loads(vardata)

            ns1 = vardata[pop1]['ns_all']
            ns2 = vardata[pop2]['ns_all']

            if ns1 < 0 or ns2 < 0:
                continue

            af1 = vardata[pop1]['af_mle']
            af2 = vardata[pop2]['af_mle']
            df_list.append(af1 - af2)

        if nfile == nfile_max:
            break

    nvar = len(df_list)
    if nvar:
        print(pop1, pop2, nvar, np.mean(df_list), np.std(df_list), sep='\t')
    else:
        print(pop1, pop2, nvar, -1, -1, sep='\t')


pop1, pop2, af_list_name, nfile_max = sys.argv[1:5]
nfile_max = int(nfile_max)
calc_pop_euclidean_distance(pop1, pop2, af_list_name, nfile_max)
