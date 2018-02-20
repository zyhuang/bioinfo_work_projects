import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet
import json


def make_dist_matrix(prov_list, prov_dist, prov_reorder_list):

    dist_dict = {}
    k = 0
    for i,p1 in enumerate(prov_list):
        for j,p2 in enumerate(prov_list):
            if j <= i:
                continue
            key = '{}:{}'.format(p1, p2)
            dist_dict[key] = prov_dist[k]
            k += 1

    nprov = len(prov_list)
    dist_matrix = np.zeros([nprov, nprov])
    for i,p1 in enumerate(prov_reorder_list):
        for j,p2 in enumerate(prov_reorder_list):
            if j == i:
                continue
            key = '{}:{}'.format(p1, p2)
            if key not in dist_dict:
                key = '{}:{}'.format(p2, p1)
            dist_matrix[i][j] = dist_dict[key]
            dist_matrix[j][i] = dist_dict[key]

    return dist_matrix


def plot_dendro_heat(varkey, linkage_list_name, out_dendro_name, out_heat_name,
                     order=False, method='ward'):

    font_size = 7
    matplotlib.rcParams.update({'font.size': font_size})

    prov_list = []
    prov_dist = []
    print('reading ' + linkage_list_name, file=sys.stderr)
    for line in open(linkage_list_name):
        key, data = line.rstrip().split('\t')
        data = json.loads(data)
        if key == varkey:
            prov_list = data['prov_list']
            prov_dist = data['prov_dist']
            break

    print('making dendrogram plot ' + out_dendro_name, file=sys.stderr)
    linkage_matrix = linkage(prov_dist, method=method,
                             optimal_ordering=order)
    fig, ax = plt.subplots()
    res = dendrogram(linkage_matrix, labels=prov_list,
                     leaf_font_size=font_size, leaf_rotation=90)
    plt.title(varkey)

    plt.savefig(out_dendro_name, bbox_inches='tight')

    print('making heat plot ' + out_heat_name, file=sys.stderr)
    prov_reorder_list = res['ivl']

    dist_matrix = make_dist_matrix(prov_list, prov_dist, prov_reorder_list)

    fig, ax = plt.subplots()
    cplot = plt.imshow(dist_matrix, cmap='Spectral', interpolation='nearest')
    cbar = fig.colorbar(cplot)
    cbar.set_label('distance')
    ax.set_xticks(np.arange(len(prov_reorder_list)))
    ax.set_yticks(np.arange(len(prov_reorder_list)))
    ax.set_xticklabels(prov_reorder_list, minor=False)
    ax.set_yticklabels(prov_reorder_list, minor=False)
    plt.xticks(rotation=90)
    plt.title(varkey)

    plt.savefig(out_heat_name, bbox_inches='tight')


def test():

    # varkey, linkage_list_name, out_dendro_name, out_heat_name = sys.argv[1:5]
    varkey, linkage_list_name, out_dendro_name, out_heat_name \
        = 'ALDH2', 'foo.list', 'plot.dendro.pdf', 'plot.heat.pdf'
    plot_dendro_heat(varkey, linkage_list_name, out_dendro_name, out_heat_name,
                     order=True, method='ward')


if __name__ == '__main__':

    varkey, linkage_list_name, out_dendro_name, out_heat_name, order, method \
        = sys.argv[1:7]

    order = True if order == 'True' else False
    plot_dendro_heat(varkey, linkage_list_name, out_dendro_name, out_heat_name,
                     order=order, method=method)
