#!/bin/bash 
hisat2-build B_anynana_v2.fa B_anynana_v2_index
for i in `ls -1 *_fwd.fq | sed 's/_fwd.fq//'`
do
hisat2 -p 24 -x B_anynana_v2_index -1 $i\_fwd.fq -2 $i\_rev.fq --dta | samtools view - -Sb | samtools sort - -@24 -o $i.bam
done
