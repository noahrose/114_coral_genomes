#!/bin/bash
#SBATCH -N1 -c2 -t 2:00:00

BASE=$(basename $1 .vcf.gz)
bcftools index $1
bcftools consensus -f $2 $1 > ref/$BASE.fna
cd ref
bwa index $BASE.fna
