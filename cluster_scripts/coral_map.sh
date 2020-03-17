#!/bin/bash
#SBATCH -N1 -c2
#SBATCH -t 12:00:00

NAME=$(basename $1 _1.txt.gz)

bwa mem -R "@RG\tID:$NAME\tSM:$NAME" -t 2 $2 ${NAME}_1.txt.gz ${NAME}_2.txt.gz | samtools sort -@2 -o $NAME.sort.bam -
samtools rmdup $NAME.sort.bam $NAME.bam
rm $NAME.sort.bam
samtools index $NAME.bam

java -Xmx4000m -jar ~/bin/picard.jar CollectAlignmentSummaryMetrics \
          R=$2 \
          I=$NAME.bam \
          O=$NAME.picard.txt

