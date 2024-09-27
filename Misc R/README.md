## This folder contains R scripts used for misc purposes

### Extract gene and genome ids from the gtf assembly file
The transcript and gene number tables contain the gene names as MSTRG only.
This scripts allows the recovery of the corresponding Bany names from the genome, which simplifies annotation and allows alignment in IGV (for sashimi plot building for example).
It imputs the transcriptome assembly file in gtf and outputs a two columns txt file with the MSTRG and Bany corresponding names. 
Output table can be opened in excel and used with VLOOKUP to add names to gene lists.

### Subsetting gene sequences from the assembly 
Uses the package seqinr.
This extracts sequences of genes of interest from the assembly and write them to a new file.
Inputs the assembly containing all genes of the transcriptome, or can input genome, inputs a list of gene names of interest, and output the sequences of the gene of interest extracted from the assembly.
Gene names in assembly and in the list of gene of interest must correspond. 
Useful for annotation of a subset of gene, like DEG and DSG.

### DEG&DSG number plot
Uses packaes ggplot2 to build the plot, and wesanderson to use the movies' color palettes.
This plots numbers of DEG and DSG next to each other as bar a plot. 

### Splicing proportions plot
Uses packaes ggplot2 to build the plot, and wesanderson to use the movies' color palettes.
This plots the percentage of each splicing type for each treatment comparison as a stacked bar chart to 100%. 
