## Purpose
Gene Ontology defines concepts/classes used to describe gene function, and relationships between these concepts. 
\
It classifies functions along three aspects: MF: Molecular Function (molecular activities of gene products), CC: Cellular Component where gene products are active and  BP: Biological Process (pathways and larger processes made up of the activities of multiple gene products). 
\
Enrichment analysis: compare the nb of GO occurrences with the expected number of occurrences for each category
\
We use clusterProfiler package (https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html ; http://bioconductor.org/packages/release/bioc/html/clusterProfiler.html).

## Process
List of files necessary for the analysis:
GO2name.csv : a list of all GO ids + description of function
\
KEGG2name.csv : list of kegg id s+ description
\
go2gene.txt : list of all GO ids + genes from the assembly
\
kegg2gene.txt : list of all kegg ids + genes from the assembly
\
mylist.csv : table of gene names of interest
\
1. Create go2gene file with tidyR and dplyr. go2gene contains all the go terms of all the genes in the assembly.
\
2. Import GO2name and mylist.
\
3. Do GO enrichment analysis.
\
4. Build plots.
