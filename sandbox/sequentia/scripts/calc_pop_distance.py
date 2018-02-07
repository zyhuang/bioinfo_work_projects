import json
import sys
import numpy as np


def calc_pop_distance(pop1, pop2, maf_min=0.01):

    af1_list = []
    af2_list = []

    for binkey in open('test.binkey.list'):
        list_name = '../../af_list/{}.af.list'.format(binkey.rstrip())

        for line in open(list_name):
            varkey, rsid, vardata = line.rstrip().split('\t')
            vardata = json.loads(vardata)

            if vardata['ALL'] < maf_min or 1-vardata['ALL'] < maf_min:
                continue

            af1_list.append(vardata[pop1])
            af2_list.append(vardata[pop2])

        break


    cov_matrix = np.cov(af1_list, af2_list)

    print(pop1, pop2, cov_matrix[0][0], cov_matrix[1][1], cov_matrix[0][1],
          cov_matrix[0][1]/np.sqrt(cov_matrix[0][0]*cov_matrix[1][1]))


pop1, pop2 = sys.argv[1:3]
calc_pop_distance(pop1, pop2)
