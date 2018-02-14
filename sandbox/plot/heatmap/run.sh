#!/bin/bash


function plot_prov () {
    nprov=$1
    echo "python3 plot_heat.py ../dendrogram/distance_pacs.list ../dendrogram/dendro.pacs.$nprov.list heat.pacs.$nprov.pdf"
    python3 plot_heat.py ../dendrogram/distance_pacs.list ../dendrogram/dendro.pacs.$nprov.list heat.pacs.$nprov.pdf
}


function plot_pop () {
    npop=$1
    echo "python3 plot_heat.py ../dendrogram/distance_1000g.list ../dendrogram/dendro.1000g.$npop.list heat.1000g.$npop.pdf"
    python3 plot_heat.py ../dendrogram/distance_1000g.list ../dendrogram/dendro.1000g.$npop.list heat.1000g.$npop.pdf
}


plot_prov 24
plot_prov 30
plot_pop 26




