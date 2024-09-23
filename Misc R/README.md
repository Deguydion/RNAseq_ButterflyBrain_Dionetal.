## This folder contains R scripts used for misc purposes

### Extract gene and genome ids from the gtf assembly file
The transcript and gene number tables contain the gene names as MSTRG only.
This scripts allows the recovery of the corresponding Bany names from the genome, which simplifies annotation and allows alignment in IGV (for sashimi plot building for example).
It imputs the transcriptome assembly file in gtf and outputs a two columns txt file with the MSTRG and Bany corresponding names. 
Output table can be opened in excel and used with VLOOKUP to add names to gene lists.
