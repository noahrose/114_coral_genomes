#!/bin/bash
>$1.faa
for i in Amil.all.maker.proteins.fasta H*.faa
    do echo ">$i" >> $1.faa && samtools faidx $i $1"-RA" | grep -v '>' >> $1.faa
done
