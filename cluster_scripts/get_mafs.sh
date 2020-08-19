#!/bin/bash
#SBATCH -N1 -c4 -t 8:00:00

angsd -minMapQ 20 -anc ref/Amilv2_ahy_consensus.fna -sites $2 -gl 1 -doMajorMinor 3 -doMaf 3 -doSaf 1  -bam $1 -out mafs/$(basename $1 .txt) -nThreads 4
realSFS mafs/$(basename $1 .txt).saf.idx -P 4 > mafs/$(basename $1 .txt).sfs
angsd -minMapQ 20 -anc ref/Amilv2_ahy_consensus.fna -sites $2 -gl 1 -pest mafs/$(basename $1 .txt).sfs -doThetas 1 -doSaf 1 -bam $1 -out mafs/$(basename $1 .txt) -nThreads 4
thetaStat do_stat mafs/$(basename $1 .txt).thetas.idx -type 2 -win 9999 -step 10000 -outnames mafs/$(basename $1 .txt)


