#!/bin/bash
#SBATCH -N1 -c1 -t 2:00:00

module load anaconda

angsd -r $2 -doCounts 1 -setMinDepth 114 -setMaxDepth 456 -minMapQ 20 -GL 1 -out $2 -nThreads 1 -doGlf 2 -doMajorMinor 1 -minMaf 0.05 -SNP_pval 1e-6 -doMaf 1 -bam $1
angsd -r $2 -minMapQ 20 -GL 1 -out $2_HA -nThreads 1 -doGlf 2 -doMajorMinor 3 -doMaf 1 -bam HA.txt -doGeno 2 -doPost 1 -postCutoff 0.8 -sites genolike.sites
angsd -r $2 -minMapQ 20 -GL 1 -out $2_HC -nThreads 1 -doGlf 2 -doMajorMinor 3 -doMaf 1 -bam HC.txt -doGeno 2 -doPost 1 -postCutoff 0.8 -sites genolike.sites
angsd -r $2 -minMapQ 20 -GL 1 -out $2_HD -nThreads 1 -doGlf 2 -doMajorMinor 3 -doMaf 1 -bam HD.txt -doGeno 2 -doPost 1 -postCutoff 0.8 -sites genolike.sites
angsd -r $2 -minMapQ 20 -GL 1 -out $2_HE -nThreads 1 -doGlf 2 -doMajorMinor 3 -doMaf 1 -bam HE.txt -doGeno 2 -doPost 1 -postCutoff 0.8 -sites genolike.sites
python ~/pcangsd/pcangsd.py -beagle $2.beagle.gz -o $2
