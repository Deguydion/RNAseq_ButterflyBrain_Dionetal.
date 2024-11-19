#!/bin/bash

# Define the number of threads to use
THREADS=24

# Loop through all BAM files in the directory
for BAM_FILE in *.bam; do
    # Extract the sample name by removing the ".bam" extension
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
    
    # Define the output GTF file name
    OUTPUT_GTF="${SAMPLE_NAME}.gtf"
    
    # Run StringTie
    echo "Processing $BAM_FILE..."
    stringtie "$BAM_FILE" -l ass_ -p "$THREADS" -o "$OUTPUT_GTF"
    
    echo "Generated $OUTPUT_GTF"
done
