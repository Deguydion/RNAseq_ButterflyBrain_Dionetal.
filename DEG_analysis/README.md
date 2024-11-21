# Purpose
This R script is used to identify differentially expressed genes (DEGs) between the exposure treatments and generate overall gene expression profiles (heatmap & PCA).
\
The analysis is performed in R, with the use of the package DEseq2 with details found at https://bioconductor.org/packages/release/bioc/html/DESeq2.html
\
Input are the abundance table (gene count matrix) obtained from stringtie and a txt frame containing sample name and corresponding treatments. 
