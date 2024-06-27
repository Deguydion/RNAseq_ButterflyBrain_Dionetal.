#!/bin/bash 
for i in 'ls -1 *_P.fq | sed 's/_P.fq//'`
do
bbsplit.sh in1=$i\_1P.fq in2=$i\_2P.fq ref='/home/nus/Desktop/NGS/bbmap/bacterial_genome/bacteria_genome.fa' basename=bbclean$i_%.fq outu1=$i\_1.fq outu2=$i\_2.fq
done
##--- END HERE ---
