#!/bin/bash

>processed_data/global_fst.txt
realSFS fst stats mafs/HA.HC.fst.idx >> processed_data/global_fst.txt
realSFS fst stats mafs/HA.HD.fst.idx >> processed_data/global_fst.txt
realSFS fst stats mafs/HA.HE.fst.idx >> processed_data/global_fst.txt
realSFS fst stats mafs/HC.HD.fst.idx >> processed_data/global_fst.txt
realSFS fst stats mafs/HC.HE.fst.idx >> processed_data/global_fst.txt
realSFS fst stats mafs/HD.HE.fst.idx >> processed_data/global_fst.txt

realSFS fst print mafs/HA.HC.fst.idx | awk '{print $1":"$2"\t"$3/$4}' | sort -k1,1 > HA.HC.snpfst.txt
realSFS fst print mafs/HA.HD.fst.idx | awk '{print $1":"$2"\t"$3/$4}' | sort -k1,1 > HA.HD.snpfst.txt
realSFS fst print mafs/HA.HE.fst.idx | awk '{print $1":"$2"\t"$3/$4}' | sort -k1,1 > HA.HE.snpfst.txt
realSFS fst print mafs/HC.HD.fst.idx | awk '{print $1":"$2"\t"$3/$4}' | sort -k1,1 > HC.HD.snpfst.txt
realSFS fst print mafs/HC.HE.fst.idx | awk '{print $1":"$2"\t"$3/$4}' | sort -k1,1 > HC.HE.snpfst.txt
realSFS fst print mafs/HD.HE.fst.idx | awk '{print $1":"$2"\t"$3/$4}' | sort -k1,1 > HD.HE.snpfst.txt

multijoin.sh *snpfst.txt snpEff.sort.txt | grep "missense_variant\|synonymous_variant" > processed_data/snpfst_sns.txt
multijoin.sh *snpfst.txt | gzip -c > processed_data/snpfst.txt.gz
