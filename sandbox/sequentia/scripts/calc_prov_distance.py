import json
import sys
import numpy as np


def calc_pop_distance(pop1, pop2, af_list_name, maf_min=0.01, nfile_max=-1,
                      ac_min=3):

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

            # drop MAF < maf_min
            if (vardata['All']['af_mle'] < maf_min or
                1-vardata['All']['af_mle'] < maf_min):
                continue

            ns1 = vardata[pop1]['ns']
            ns2 = vardata[pop2]['ns']

            # drop sites with insufficient sample size
            if ns1*maf_min < ac_min or ns2*maf_min < ac_min:
                continue

            af1 = max(vardata[pop1]['af_mle'], 0)
            af2 = max(vardata[pop2]['af_mle'], 0)
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


pop1, pop2, af_list_name, maf_min, nfile_max = sys.argv[1:6]
maf_min = float(maf_min)
nfile_max = int(nfile_max)
calc_pop_distance(pop1, pop2, af_list_name, maf_min, nfile_max)
