#!/bin/bash
#SBATCH -t 48:00:00 -N 1 -c 3

bcftools mpileup -q20 -f $2 -b $1 -Ou -BI -a AD,DP,INFO/AD | bcftools call -vmOu  | bcftools view -v snps -q 0.5:alt1 | bcftools norm -Ou -m  - | bcftools norm -Oz -d snps -  > $(basename $1 .txt)_consensus.vcf.gz


