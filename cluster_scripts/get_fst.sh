#!/bin/bash
#SBATCH -N1 -c4 -t 12:00:00

pop1=$1
pop2=$2
realSFS mafs/$pop1.saf.idx mafs/$pop2.saf.idx -P 4 > mafs/$pop1.$pop2.ml
#prepare the fst for easy window analysis etc
realSFS fst index mafs/$pop1.saf.idx mafs/$pop2.saf.idx -sfs mafs/$pop1.$pop2.ml -fstout mafs/$pop1.$pop2 -P 4
#get the global estimate
realSFS fst stats mafs/$pop1.$pop2.fst.idx 
realSFS fst stats2 mafs/$pop1.$pop2.fst.idx -win 9999 -step 10000 -type 2 > processed_data/$pop1.$pop2.10kb.fst
echo "getting dxy"
getDxy.pl --pop1maf <(zcat mafs/$pop1.mafs.gz ) --pop2maf <(zcat mafs/$pop2.mafs.gz ) --minInd 1 | awk 'NR > 1{print $1,$2-1,$2,$3}' OFS="\t" | bedtools map -a ref/ahy_10kb_intervals.bed -b - -c 4 -o sum > processed_data/$pop1.$pop2.dxy
