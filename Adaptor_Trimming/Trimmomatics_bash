#!/bin/bash 
for i in `ls -1 *_1.fq.gz | sed 's/_1.fq.gz//'` 
do 
echo "Run trimmomatic for ${i}"
java -jar /path/to/trimmomatic.jar PE -threads 24 $i\_1.fq.gz $i\_2.fq.gz $i\_1P.fq.gz $i\_1U.fq.gz $i\_2P.fq.gz $i\_2U.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:true MAXINFO:40:0.2 
done
##--- END HERE --- 
