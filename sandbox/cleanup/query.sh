#!/bin/bash

#PBS -q scavenger
#PBS -A scavenger
#PBS -o log
#PBS -V
#PBS -j oe
#PBS -N query
#PBS -l nodes=1:ppn=1,mem=5G
#PBS -t 1-3

#-q analysis
#-A proj-fy0006
#-l walltime=13:00:00:00

#cd $PBS_O_WORKDIR

#
PBS_ARRAYID=1

t1=$(date +%s)

incr=3
istart=$[(PBS_ARRAYID-1)*incr+1]
iend=$[PBS_ARRAYID*incr]

timestamp=$(date +%Y%m%d)
mkdir -p $timestamp
for list in $(ls sfns?.query.txt | awk 'NR>='$istart' && NR<='$iend''); do
    cp $list $timestamp/
    echo "search.pl -o -f $timestamp/$list"
    search.pl -o -f $timestamp/$list
done



t2=$(date +%s)
dt=$[t2-t1]
echo "% Job $PBS_ARRAYID is done in $(date -ud @$dt +%j:%T | awk -F\: '{printf("%03d:%02d:%02d:%02d\n",$1-1,$2,$3,$4)}')"
