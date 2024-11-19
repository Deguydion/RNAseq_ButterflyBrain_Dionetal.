### Purpose: 
These scripts find alternative splicing events in the transcriptome and identify differentially spliced genes (DSGs) between two treatments.
We used rMATS https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md

### Steps:
1. run rmats using "statoff" to calculate the percentage of alternative splicing without any apriori on the groups. Used for building the heatmap
Inputs: list of samples (here called brain_samples.txt), as a list with file names separated by commas, merged gtf asssembly from all the samples
2. calculate the percenage of alternative splicing for all genes in all samples and do treatment comparisons
3: perform the comparisons between the different groups. 
A folder is created for each comparison and we code the treatments by calling them groups 1 and 2, and each reps within treatments are called by a number corresponding to the order they appear in the all_samples.txt file
4: perform statistics on each comparisons
The remaining filtering and procesing steps are run in R (see DSG.R). 
