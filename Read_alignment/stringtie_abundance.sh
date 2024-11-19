#!/bin/bash
# Number of threads and reference annotation file
THREADS=24
ANNOTATION="Bicyclus_anynana_BaGv2.gff3"

# Loop through all BAM files in the directory
for BAM_FILE in *.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
    OUTPUT_GTF="abund_${SAMPLE_NAME}.gtf"
    stringtie -e -B -p "$THREADS" -G "$ANNOTATION" -o "$OUTPUT_GTF" "$BAM_FILE"
done
