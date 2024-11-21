## Purpose:

These scripts find alternative splicing events in the transcriptome and identify differentially spliced genes (DSGs) between two treatments.
We used rMATS https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md

## Process:
  
1. Run rmats using "statoff" to calculate the percentage of alternative splicing without any apriori on the groups. Used for building the heatmap.
Inputs are a list of samples (here called brain_samples.txt) with file names separated by commas; and the merged gtf asssembly that contains transcripts.
2. Run rmats with "task both" to detect splicing events and compare the treatments to get DSG.
3. Perform the comparisons between the different groups. A folder is created for each comparison and we code the treatments by calling them groups 1 and 2 etc. Each reps within treatments are called by a number corresponding to the order they appear in the brain_samples.txt file.
4. Perform statistics on each comparison using "task stats".
5. The remaining filtering and processing steps after the statoff and after task stats are run in R (see Splicing.R).
6. To build the sashimi plots in IGV (https://igv.org/), we need to index each .bam library into a.bai file using samtools (http://www.htslib.org/). Run bai_indexing.samtools.sh to index all files, and load the .bai files into IGV.
