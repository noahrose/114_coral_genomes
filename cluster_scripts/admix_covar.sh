#!/bin/bash
#SBATCH -N1 -c16 -t 8:00:00

angsd -doCounts 1 -setMinDepth 114 -setMaxDepth 456 -minMapQ 20 -GL 1 -out genolike -nThreads 20 -doGlf 2 -doMajorMinor 1 -minMaf 0.05 -SNP_pval 1e-6 -doMaf 1 -bam $1
NGSadmix -likes genolike.beagle.gz -K $2 -P 16 -o admix -minMaf 0.05
python ~/pcangsd/pcangsd.py -beagle genolike.beagle.gz -o pcangsd -threads 10
