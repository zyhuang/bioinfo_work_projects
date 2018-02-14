import json
import sys
import numpy as np


def calc_pop_angular_distance(pop1, pop2, af_list_name, maf_min=0.01,
                              nfile_max=-1):

    af1_list = []
    af2_list = []

    nfile = 0
    nvar = 0
    for list_name in open(af_list_name):
        list_name = list_name.rstrip()
        # print('reading ' + list_name, file=sys.stderr)

        nfile += 1
        for line in open(list_name):
            varkey, rsid, vardata = line.rstrip().split('\t')
            vardata = json.loads(vardata)

            if vardata['ALL'] < maf_min or 1-vardata['ALL'] < maf_min:
                continue

            af1_list.append(vardata[pop1])
            af2_list.append(vardata[pop2])
            nvar += 1

        if nfile == nfile_max:
            break

    cov_matrix = np.cov(af1_list, af2_list)

    print(pop1, pop2, nvar,
          cov_matrix[0][0], cov_matrix[1][1], cov_matrix[0][1],
          cov_matrix[0][1]/np.sqrt(cov_matrix[0][0]*cov_matrix[1][1]),
          sep='\t')


def calc_pop_euclidean_distance(pop1, pop2, af_list_name, maf_min=0.01,
                                nfile_max=-1):

    df_list = []
    nfile = 0
    for list_name in open(af_list_name):
        list_name = list_name.rstrip()

        nfile += 1
        for line in open(list_name):
            varkey, rsid, vardata = line.rstrip().split('\t')
            vardata = json.loads(vardata)

            if vardata['ALL'] < maf_min or 1-vardata['ALL'] < maf_min:
                continue

            df_list.append(vardata[pop1] - vardata[pop2])

        if nfile == nfile_max:
            break

    nvar = len(df_list)
    print(pop1, pop2, nvar, np.mean(df_list), np.std(df_list), sep='\t')



pop1, pop2, af_list_name, maf_min, nfile_max = sys.argv[1:6]
maf_min = float(maf_min)
nfile_max = int(nfile_max)
calc_pop_euclidean_distance(pop1, pop2, af_list_name, maf_min, nfile_max)
