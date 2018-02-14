import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


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


def plot_data(data_matrix, data_label, out_pdf_name):

    matplotlib.rcParams.update({'font.size': 7})

    fig, ax = plt.subplots()

    cplot = plt.imshow(data_matrix, cmap='Spectral_r',
                       interpolation='nearest')
    cbar = fig.colorbar(cplot)
    cbar.set_label('similarity')
    ax.set_xticks(np.arange(len(data_label)))
    ax.set_yticks(np.arange(len(data_label)))
    ax.set_xticklabels(data_label, minor=False)
    ax.set_yticklabels(data_label, minor=False)
    plt.xticks(rotation=90)

    # plt.show()
    plt.savefig(out_pdf_name, bbox_inches='tight')


def plot_heat():

    data_file_name, data_label_name, out_pdf_name = sys.argv[1:4]
    # data_file_name = 'distance.list'
    # data_label_name = 'prov_order.list'
    data_matrix, data_label = proc_data(data_file_name, data_label_name)
    plot_data(data_matrix, data_label, out_pdf_name)


if __name__ == '__main__':

    plot_heat()
