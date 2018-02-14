#!/bin/bash

function plot_prov () {
    nprov=$1
    echo python3 plot_dendro.py distance_pacs.list order_pacs.$nprov.list dendro.pacs.$nprov.pdf dendro.pacs.$nprov.list
    python3 plot_dendro.py distance_pacs.list order_pacs.$nprov.list dendro.pacs.$nprov.pdf dendro.pacs.$nprov.list
}

function plot_pop () {
    npop=$1
    echo python3 plot_dendro.py distance_1000g.list order_1000g.$npop.list dendro.1000g.$npop.pdf dendro.1000g.$npop.list
    python3 plot_dendro.py distance_1000g.list order_1000g.$npop.list dendro.1000g.$npop.pdf dendro.1000g.$npop.list
}

plot_prov 24
plot_prov 30
plot_pop 26




