#!/bin/bash
#SBATCH -N1 -c2
#SBATCH -t 16:00:00

NAME=$(basename $1 _1.txt.gz)
echo $NAME
bwa mem -R "@RG\tID:$NAME\tSM:$NAME" -t 2 $2 ${NAME}_1.txt.gz ${NAME}_2.txt.gz | samtools sort -@2 -o $3/$NAME.sort.bam -
samtools rmdup $3/$NAME.sort.bam $3/$NAME.bam
rm $3/$NAME.sort.bam
samtools index $3/$NAME.bam

