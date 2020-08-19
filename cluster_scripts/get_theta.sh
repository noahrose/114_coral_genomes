#!/bin/bash
#SBATCH -N1 -c16 -t 8:00:00

angsd -doCounts 1 -setMinDepth $2 -setMaxDepth $3 -minMapQ 20 -anc ref/Amilv2_ahy_consensus.fna -gl 1 -doSaf 1  -bam $1 -out mafs/allsites_$(basename $1 .txt) -nThreads 16
realSFS mafs/$(basename $1 .txt).saf.idx -P 16 > mafs/allsites_$(basename $1 .txt).sfs
angsd -doCounts 1 -setMinDepth $2 -setMaxDepth $3 -minMapQ 20 -anc ref/Amilv2_ahy_consensus.fna -gl 1 -pest mafs/allsites_$(basename $1 .txt).sfs -doThetas 1 -doSaf 1 -bam $1 -out mafs/allsites_$(basename $1 .txt) -nThreads 16
thetaStat do_stat mafs/allsites_$(basename $1 .txt).thetas.idx -type 2 -win 9999 -step 10000 -outnames mafs/allsites_$(basename $1 .txt)


