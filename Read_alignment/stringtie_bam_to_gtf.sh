#!/bin/bash

THREADS=24

# Loop through all BAM files in the directory
for BAM_FILE in *.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
    OUTPUT_GTF="${SAMPLE_NAME}.gtf"
    stringtie "$BAM_FILE" -l ass_ -p "$THREADS" -o "$OUTPUT_GTF"
done
