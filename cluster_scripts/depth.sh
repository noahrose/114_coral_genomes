#!/bin/bash

samtools depth -r $1 -Q20 *.bam |\
    awk -v OFS="\t" '{for(i=3;i<=NF;i++) t+=$i; print $1,$2,t; t=0}' |\
    awk -v OFS="\t" '$3>113 && $4<457 {print $1,$2-1,$2}' |\
    bedtools coverage -nonamecheck -b - -a ref/ahy_10kb_intervals.unsort.bed |\
    awk -v chr=$1 '$1 == chr' > depth/${1}_depth.txt
