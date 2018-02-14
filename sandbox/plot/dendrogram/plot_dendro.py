import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage


def proc_data(data_file_name, data_label_name):

    data_label = []
    for line in open(data_label_name):
        if line.startswith('#'):
            continue
        data_label.append(line.rstrip())

    nitem = len(data_label)
    data_matrix = np.zeros([nitem, nitem])
    label_index = {x:i for i,x in enumerate(data_label)}

    for line in open(data_file_name):
        data = line.rstrip().split('\t')
        item1, item2, dist = data[0], data[1], float(data[-1])
        if item1 not in label_index or item2 not in label_index:
            continue
        idx1 = label_index[item1]
        idx2 = label_index[item2]
        data_matrix[idx1][idx2] = dist
        data_matrix[idx2][idx1] = dist

    return data_matrix, data_label


def plot_data(data_matrix, data_label, out_pdf_name, out_order_name):

    font_size = 7
    matplotlib.rcParams.update({'font.size': font_size})

    fig, ax = plt.subplots()

    distVec = squareform(data_matrix)
    max_dist = (1-max(distVec)) * 1.1
    mylinkage = linkage(distVec, 'ward')
    r = dendrogram(mylinkage, labels=data_label, leaf_font_size=font_size,
               leaf_rotation=90)

    fout = open(out_order_name, 'w')
    for x in r['ivl']:
        print(x, file=fout)
    fout.close()

    plt.savefig(out_pdf_name, bbox_inches='tight')


def plot_dendro():

    data_file_name, data_label_name, out_pdf_name, out_order_name = sys.argv[1:5]

    data_matrix, data_label = proc_data(data_file_name, data_label_name)
    plot_data(data_matrix, data_label, out_pdf_name, out_order_name)


if __name__ == '__main__':

    plot_dendro()
