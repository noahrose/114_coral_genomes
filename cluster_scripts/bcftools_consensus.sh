#!/bin/bash
#SBATCH -N1 -c2 -t 2:00:00

bcftools index $1
bcftools consensus -f $2 $1 > ref/$3.fna
cd ref
bwa index $3.fna
