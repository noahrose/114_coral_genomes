#!/bin/bash
#SBATCH -N1 -c1 -t 2:00:00

realSFS fst print mafs/$1.$2.fst.idx | awk '{print $1,$2-1,$2,$3,$4}' OFS="\t" | bedtools map -nonamecheck -a ref/ahy_10kb_intervals.bed -b - -c 4,5,5 -o sum,sum,count | sort -k1,1 -k2,2n > processed_data/$1.$2.fst.components.txt
