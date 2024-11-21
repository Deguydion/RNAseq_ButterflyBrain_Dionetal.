#GENE ENRICHMENT ANALYSIS DONE IN R

#We use clusterProfiler package 
#Enrichment analysis: compare the nb of GO occurrences with the expected number of occurrences for each category

setwd("path/to/data/data")

#Creating the go2gene file with tidyR. go2gene contains all the go terms of all the genes in the assembly
library(tidyr)
library(dplyr)
data_frame<-read.csv('go2gene_ass_prep.csv',header=F,sep=",",fileEncoding = 'UTF-8-BOM') #file encoding removes the weird signs before the name of the first row
head(data_frame)
newframe<-pivot_longer(data_frame,cols=2:24,names_to="GOtype",values_to="GO") #pivot the frame
newframe<-subset(newframe,select=-GOtype) #removes the columns 
newframe<- newframe %>%  #fills the empty cells with NA, and remove the rows with NAs
mutate(across(everything(),~na_if(., ""))) 
newframe<-na.omit(newframe)
head(newframe) 
write.csv(as.data.frame(newframe),file="go2gene_assembly.csv") #Remove first column & line in excel, make GO first column, remove space before GO

# Import the GOannotations.
go2name_data=read.csv("GO2name.csv",header = FALSE,sep=",") #all go and corresponding names
GO2gene_data=read.table("go2gene_assembly.csv",header = FALSE, sep=",") #ref go in all genome or transcriptome
head(go2name_data)
head(GO2gene_data)

# Import gene list of interest
mylist_all <- read.csv("brain_list.csv",header = F,sep=",",fileEncoding = 'UTF-8-BOM')
names(mylist_all)=c("dsg","deg") #for comparing 2 samples
names(mylist_all)=c("naives versus exposed","males versus females") #for comparing males with females

mylist_all<-as.character(mylist_all[,1]) #converted as characters
head(mylist_all)

# GO enrichment analysis
background=GO2gene_data
#for comparison between 2 treatments
clusterenrich <- compareCluster(mylist_all, fun='enricher',TERM2GENE=background,TERM2NAME=go2name_data,pvalueCutoff = 0.1, pAdjustMethod = "BH") # pvalueCutoff can be lifted (changed to 1) if very few terms can be enriched. 
clusterenrich <- 0389
#For one set of genes
enricher(mylist_all,TERM2GENE=background,TERM2NAME=go2name_data,pvalueCutoff= 1, pAdjustMethod="BH",qvalueCutoff = 0.1) 
clusterenrich
write.csv(as.data.frame(clusterenrich),file="GO_ComparisonOfInterest.csv") # Export all GO enrichment results.


#Building the plots
dot=dotplot(clusterenrich,showCategory=8,includeAll=TRUE,color="pvalue")
dot
dot=dot+scale_color_continuous(low='orange', high='navy')
dot=dot+ggpubr::rotate_x_text()
dot=dot+scale_y_discrete(label = function(x) stringr::str_trunc(x, 50))
dot=dot+scale_x_discrete(label=c("naives versus exposed","males versus females"))
dot

# KEGG enrichment analysis
background_kegg=kegg2gene_data
clusterenrich_kegg <- compareCluster(mylist_all, fun='enricher',TERM2GENE=background_kegg,TERM2NAME=kegg2name_data,pvalueCutoff = 1, pAdjustMethod = "BH")
# Export all KEGG enrichment results
write.csv(as.data.frame(clusterenrich_kegg),file="KEGG_DSG.csv")
# Draw dot plot
dot=dotplot(clusterenrich_kegg, showCategory=100,includeAll=TRUE,color="pvalue")
dot=dot+scale_color_continuous(low='orange', high='navy')
dot=dot+ggpubr::rotate_x_text()
dot=dot+scale_yiscrete(label = function(x) stringr::str_trunc(x, 50))
dot=dot+scale_x_discrete(label=c("Wr60","PP50","P15")) #In this case one of the comparisons, P50, has no enriched KEGG terms
dot
# Save image
ggsave(dot,width=6.5,height=4.5,filename='KEGG_DSG.tiff')
ggsave(dot,width=6.5,height=4.5,filename='KEGG_DSG.eps')
