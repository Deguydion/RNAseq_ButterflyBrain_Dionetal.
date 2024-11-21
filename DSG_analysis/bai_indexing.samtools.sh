#!/bin/bash

for bam_file in *.bam; do
    base_name="${bam_file%.bam}"
    bai_file="${base_name}.bam.bai"
    samtools index "$bam_file" "$bai_file"
    
done
