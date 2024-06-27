#!/bin/bash

# Define paths to the SortMeRNA database files and the input directory
DATABASE_DIR="/AMLab_Ext1/Emilie/AB_transcr_2020/sortmerna/sortmerna-4.2.0/data/rRNA_databases"
INPUT_DIR="/AMLab_Ext1/Emilie/AB_transcr_2020/sortmerna"
OUTPUT_DIR="/AMLab_Ext1/Emilie/AB_transcr_2020/sortmerna/rRNA"

# Loop through each pair of files in the input directory
for INPUT_FILE_R1 in $INPUT_DIR/*_R1.fastq; do
    # Get the corresponding R2 file
    INPUT_FILE_R2="${INPUT_FILE_R1/_R1.fastq/_R2.fastq}"

    # Get the base name of the file (without the directory and _R1/_R2.fastq extension)
    BASE_NAME=$(basename "$INPUT_FILE_R1" _R1.fastq)

    # Define output files
    OUTPUT_RRNA_R1="$OUTPUT_DIR/${BASE_NAME}_rRNA_R1.fastq"
    OUTPUT_RRNA_R2="$OUTPUT_DIR/${BASE_NAME}_rRNA_R2t.fastq"
    OUTPUT_NON_RRNA_R1="$OUTPUT_DIR/${BASE_NAME}_sortme_R1.fastq"
    OUTPUT_NON_RRNA_R2="$OUTPUT_DIR/${BASE_NAME}_sortme_R2.fastq"

    # Run SortMeRNA
    sortmerna --ref "$DATABASE_DIR/silva-bac-16s-id90.fasta" \
              --reads "$INPUT_FILE_R1" \
              --reads "$INPUT_FILE_R2" \
              --aligned "$OUTPUT_RRNA_R1" "$OUTPUT_RRNA_R2" \
              --other "$OUTPUT_NON_RRNA_R1" "$OUTPUT_NON_RRNA_R2" \
              --paired_in \
              --fastx \
              --log -v

    echo "Processed $INPUT_FILE_R1 and $INPUT_FILE_R2"
done

echo "All samples processed."
