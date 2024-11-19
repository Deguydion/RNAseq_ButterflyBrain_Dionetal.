#DSG ANALYSES FOLLOWING RMATS OUTPUTS in R 
#Code adapted from Shen Tian (Tian S, Monteiro A (2022) A Transcriptomic Atlas Underlying Developmental Plasticity of Seasonal Forms of Bicyclus anynana Butterflies. Mol Biol Evol 39.)

library(reshape)
setwd("C:/path/to/out_MBvsNB")

#1. DETECT ALTERNATIVE SPLICING EVENTS WITHOUT COMPARISON IN RMATS
#Useful to build heatmaps

#Load packages
library(reshape)
library(grid)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(yarrr)
library(plyr)
library(goseq)
library(clusterProfiler)
library(tximport)
library(GenomicFeatures)
library(ggfortify)
library(rgl)

# Before importing the rMTAS outpout files (JC.txt), remove the first row (header) and rename as fixed.txt
# Here we detect each type of alternative splicing events/genes across all clean samples 

# Importing the fixed.JC.txt 
setwd("C:/Users/molen/Desktop/transcriptomics_analysis/Splicing_analysis/splicing_brains_april23/New folder")
list.files("C:/Users/molen/Desktop/transcriptomics_analysis/Splicing_analysis/splicing_brains_april23/New folder")

#A3SS
A3SS=read.table("A3SS_all_fixed.txt",header=FALSE)
colnames(A3SS)=c("ID","GeneID","geneSymbol","chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
A3SS=transform(A3SS,IJC=colsplit(IJC,split=",",names=c('MB01','MB02','MB03','MB04','MB05','NB01','NB02','NB03','NB04','NB05','NBB01','NBB02','NBB03','NBB04','NBB05','NCB01','NCB02','NCB03','NCB04','NCB05','WTB01','WTB02','WTB03','WTB04','WTB05')))
A3SS=transform(A3SS,SJC=colsplit(SJC,split=",",names=c('MB01','MB02','MB03','MB04','MB05','NB01','NB02','NB03','NB04','NB05','NBB01','NBB02','NBB03','NBB04','NBB05','NCB01','NCB02','NCB03','NCB04','NCB05','WTB01','WTB02','WTB03','WTB04','WTB05')))
A3SS=transform(A3SS,IncLevel=colsplit(IncLevel,split=",",names=c('MB01','MB02','MB03','MB04','MB05','NB01','NB02','NB03','NB04','NB05','NBB01','NBB02','NBB03','NBB04','NBB05','NCB01','NCB02','NCB03','NCB04','NCB05','WTB01','WTB02','WTB03','WTB04','WTB05')))

A3SS$IJC_sum=rowSums(A3SS$IJC)
A3SS$SJC_sum=rowSums(A3SS$SJC)
A3SS$IncLevel_mean=rowMeans(A3SS$IncLevel,na.rm = TRUE)
write.csv(A3SS,file="A3SS_all.csv")
A3SS_count_over_20=A3SS[A3SS$IJC_sum>20&A3SS$SJC_sum>20,]
write.csv(A3SS_count_over_20,file="A3SS_all_count_over_20.csv")
A3SS_PSI_10=A3SS_count_over_20[A3SS_count_over_20$IncLevel_mean>0.1&0.9>A3SS_count_over_20$IncLevel_mean,]
write.csv(A3SS_PSI_10,file="A3SS_all_PSI_10.csv")

#A5SS
A5SS=read.table("A5SS_all_fixed.txt",header=FALSE)
colnames(A5SS)=c("ID","GeneID","geneSymbol","chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
A5SS=transform(A5SS,IJC=colsplit(IJC,split=",",names=c('MB01','MB02','MB03','MB04','MB05','NB01','NB02','NB03','NB04','NB05','NBB01','NBB02','NBB03','NBB04','NBB05','NCB01','NCB02','NCB03','NCB04','NCB05','WTB01','WTB02','WTB03','WTB04','WTB05')))
A5SS=transform(A5SS,SJC=colsplit(SJC,split=",",names=c('MB01','MB02','MB03','MB04','MB05','NB01','NB02','NB03','NB04','NB05','NBB01','NBB02','NBB03','NBB04','NBB05','NCB01','NCB02','NCB03','NCB04','NCB05','WTB01','WTB02','WTB03','WTB04','WTB05')))
A5SS=transform(A5SS,IncLevel=colsplit(IncLevel,split=",",names=c('MB01','MB02','MB03','MB04','MB05','NB01','NB02','NB03','NB04','NB05','NBB01','NBB02','NBB03','NBB04','NBB05','NCB01','NCB02','NCB03','NCB04','NCB05','WTB01','WTB02','WTB03','WTB04','WTB05')))

A5SS$IJC_sum=rowSums(A5SS$IJC)
A5SS$SJC_sum=rowSums(A5SS$SJC)
A5SS$IncLevel_mean=rowMeans(A5SS$IncLevel,na.rm = TRUE)
write.csv(A5SS,file="A5SS_all.csv")
A5SS_count_over_20=A5SS[A5SS$IJC_sum>20&A5SS$SJC_sum>20,]
write.csv(A5SS_count_over_20,file="A5SS_all_count_over_20.csv")
A5SS_PSI_10=A5SS_count_over_20[A5SS_count_over_20$IncLevel_mean>0.1&0.9>A5SS_count_over_20$IncLevel_mean,]
write.csv(A5SS_PSI_10,file="A5SS_all_PSI_10.csv")

#SE
SE=read.table("SE_all_fixed.txt",header=FALSE)
colnames(SE)=c("ID","GeneID","geneSymbol","chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
SE=transform(SE,IJC=colsplit(IJC,split=",",names=c('MB01','MB02','MB03','MB04','MB05','NB01','NB02','NB03','NB04','NB05','NBB01','NBB02','NBB03','NBB04','NBB05','NCB01','NCB02','NCB03','NCB04','NCB05','WTB01','WTB02','WTB03','WTB04','WTB05')))
SE=transform(SE,SJC=colsplit(SJC,split=",",names=c('MB01','MB02','MB03','MB04','MB05','NB01','NB02','NB03','NB04','NB05','NBB01','NBB02','NBB03','NBB04','NBB05','NCB01','NCB02','NCB03','NCB04','NCB05','WTB01','WTB02','WTB03','WTB04','WTB05')))
SE=transform(SE,IncLevel=colsplit(IncLevel,split=",",names=c('MB01','MB02','MB03','MB04','MB05','NB01','NB02','NB03','NB04','NB05','NBB01','NBB02','NBB03','NBB04','NBB05','NCB01','NCB02','NCB03','NCB04','NCB05','WTB01','WTB02','WTB03','WTB04','WTB05')))

SE$IJC_sum=rowSums(SE$IJC)
SE$SJC_sum=rowSums(SE$SJC)
SE$IncLevel_mean=rowMeans(SE$IncLevel,na.rm = TRUE)
write.csv(SE,file="SE_all.csv")
SE_count_over_20=SE[SE$IJC_sum>20&SE$SJC_sum>20,]
write.csv(SE_count_over_20,file="SE_all_count_over_20.csv")
SE_PSI_10=SE_count_over_20[SE_count_over_20$IncLevel_mean>0.1&0.9>SE_count_over_20$IncLevel_mean,]
write.csv(SE_PSI_10,file="SE_all_PSI_10.csv")

#MXE
MXE=read.table("MXE_all_fixed.txt",header=FALSE)
colnames(MXE)=c("ID","GeneID","geneSymbol","chr","strand","1stExonStart_0base","1stExonEnd","2stExonStart_0base","2stExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
MXE=transform(MXE,IJC=colsplit(IJC,split=",",names=c('MB01','MB02','MB03','MB04','MB05','NB01','NB02','NB03','NB04','NB05','NBB01','NBB02','NBB03','NBB04','NBB05','NCB01','NCB02','NCB03','NCB04','NCB05','WTB01','WTB02','WTB03','WTB04','WTB05')))
MXE=transform(MXE,SJC=colsplit(SJC,split=",",names=c('MB01','MB02','MB03','MB04','MB05','NB01','NB02','NB03','NB04','NB05','NBB01','NBB02','NBB03','NBB04','NBB05','NCB01','NCB02','NCB03','NCB04','NCB05','WTB01','WTB02','WTB03','WTB04','WTB05')))
MXE=transform(MXE,IncLevel=colsplit(IncLevel,split=",",names=c('MB01','MB02','MB03','MB04','MB05','NB01','NB02','NB03','NB04','NB05','NBB01','NBB02','NBB03','NBB04','NBB05','NCB01','NCB02','NCB03','NCB04','NCB05','WTB01','WTB02','WTB03','WTB04','WTB05')))

MXE$IJC_sum=rowSums(MXE$IJC)
MXE$SJC_sum=rowSums(MXE$SJC)
MXE$IncLevel_mean=rowMeans(MXE$IncLevel,na.rm = TRUE)
write.csv(MXE,file="MXE_all.csv")
MXE_count_over_20=MXE[MXE$IJC_sum>20&MXE$SJC_sum>20,]
write.csv(MXE_count_over_20,file="MXE_all_count_over_20.csv")
MXE_PSI_10=MXE_count_over_20[MXE_count_over_20$IncLevel_mean>0.1&0.9>MXE_count_over_20$IncLevel_mean,]
write.csv(MXE_PSI_10,file="MXE_all_PSI_10.csv")

#RI
RI=read.table("RI_all_fixed.txt",header=FALSE)
colnames(RI)=c("ID","GeneID","geneSymbol","chr","strand","riExonStart_0base","riExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
RI=transform(RI,IJC=colsplit(IJC,split=",",names=c('MB01','MB02','MB03','MB04','MB05','NB01','NB02','NB03','NB04','NB05','NBB01','NBB02','NBB03','NBB04','NBB05','NCB01','NCB02','NCB03','NCB04','NCB05','WTB01','WTB02','WTB03','WTB04','WTB05')))
RI=transform(RI,SJC=colsplit(SJC,split=",",names=c('MB01','MB02','MB03','MB04','MB05','NB01','NB02','NB03','NB04','NB05','NBB01','NBB02','NBB03','NBB04','NBB05','NCB01','NCB02','NCB03','NCB04','NCB05','WTB01','WTB02','WTB03','WTB04','WTB05')))
RI=transform(RI,IncLevel=colsplit(IncLevel,split=",",names=c('MB01','MB02','MB03','MB04','MB05','NB01','NB02','NB03','NB04','NB05','NBB01','NBB02','NBB03','NBB04','NBB05','NCB01','NCB02','NCB03','NCB04','NCB05','WTB01','WTB02','WTB03','WTB04','WTB05')))

RI$IJC_sum=rowSums(RI$IJC)
RI$SJC_sum=rowSums(RI$SJC)
RI$IncLevel_mean=rowMeans(RI$IncLevel,na.rm = TRUE)
write.csv(RI,file="RI_all.csv")
RI_count_over_20=RI[RI$IJC_sum>20&RI$SJC_sum>20,]
write.csv(RI_count_over_20,file="RI_all_count_over_20.csv")
RI_PSI_10=RI_count_over_20[RI_count_over_20$IncLevel_mean>0.1&0.9>RI_count_over_20$IncLevel_mean,]
write.csv(RI_PSI_10,file="RI_all_PSI_10.csv")

# Number of genes associated with the splicing events
genes_A3SS=unique(as.data.frame(A3SS_PSI_10[,2])) #Select all unique gene names
colnames(genes_A3SS)="gene" #Add them to a data frame
genes_A5SS=unique(as.data.frame(A5SS_PSI_10[,2]))
colnames(genes_A5SS)="gene"
genes_MXE=unique(as.data.frame(MXE_PSI_10[,2]))
colnames(genes_MXE)="gene"
genes_SE=unique(as.data.frame(SE_PSI_10[,2]))
colnames(genes_SE)="gene"
genes_RI=unique(as.data.frame(RI_PSI_10[,2]))
colnames(genes_RI)="gene"
genes_all_clean=rbind(genes_A3SS,genes_A5SS) #Combine the data frames
genes_all_clean=rbind(genes_all_clean,genes_MXE)
genes_all_clean=rbind(genes_all_clean,genes_SE)                 
genes_all_clean=rbind(genes_all_clean,genes_RI) 
genes_all_clean=unique(genes_all_clean)

# Create column data to plot all filtered splicing events across all clean samples
A3SS_PSI10=read.csv(file="A3SS_all_PSI_10.csv")
A5SS_PSI10=read.csv(file="A5SS_all_PSI_10.csv")
MXE_PSI10=read.csv(file="MXE_all_PSI_10.csv")
RI_PSI10=read.csv(file="RI_all_PSI_10.csv")
SE_PSI10=read.csv(file="SE_all_PSI_10.csv")
A3SS_PSI10$GeneID <- paste("A3SS", substring(A3SS_PSI10$GeneID,6), sep="_")
A5SS_PSI10$GeneID <- paste("A5SS", substring(A5SS_PSI10$GeneID,6), sep="_")
MXE_PSI10$GeneID <- paste("MXE", substring(MXE_PSI10$GeneID,6), sep="_")
RI_PSI10$GeneID <- paste("RI", substring(RI_PSI10$GeneID,6), sep="_")
SE_PSI10$GeneID <- paste("SE", substring(SE_PSI10$GeneID,6), sep="_")

colnames(A3SS_PSI10) #Check column range to select
colnames(A5SS_PSI10)
colnames(MXE_PSI10)
colnames(SE_PSI10)
colnames(RI_PSI10)
all_PSI10=rbind(A3SS_PSI10[68:92],A5SS_PSI10[68:92],MXE_PSI10[70:94],SE_PSI10[68:92],RI_PSI10[68:92]) #Edit the column nrange accordingly
colnames(all_PSI10)=substring(colnames(all_PSI10),10)
all_PSI10_no_na=all_PSI10[complete.cases(all_PSI10),]
colnames(all_PSI10_no_na)

Treatment <- factor(rep(c("Males","Naives","NB1","NB2","Wt"),each = 5)) #Edit the treatment name and sample size
Sex <- c(rep("Males", 5), rep("Females", 20)) #Edit the second factor, here sex
col_data <- data.frame(row.names = substring(colnames(A3SS_PSI10[68:92]),10),Treatment,Sex) #Edit column numbers
col_data


# Plot heatmap
ann_colors = list(Sex = c(Males="#003300", Females="#660000"),Treatment=c(Males="lightgrey",Naives="black",NB1="#009988",NB2="#0177BB",Wt="#EE7733"))
temp_hm=pheatmap(all_PSI10_no_na,border=NA,scale = "row",annotation_colors = ann_colors,color = colorRampPalette(c("purple","purple", "black", "yellow", "yellow"))(100), annotation=col_data, show_rownames=FALSE,show_colnames = TRUE,annotation_names_col=FALSE,annotation_legend = TRUE)
ggsave(temp_hm,width=6,height=6,filename='all_clean.tiff')
ggsave(temp_hm,width=6,height=6,filename='all_clean.eps')


#2. BUILD THE DIFFERENTIALLY SPLICED GENE TABLES 
#Example shown is MvsN, same used for other comparisons

#We detect each type of differentially spliced genes between 2 treatments
#A3SS: alternative 3' splice sites
A3SS=read.table("A3SS.MATS.JC.txt",header=TRUE)
A3SS=transform(A3SS,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
A3SS=transform(A3SS,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
A3SS=transform(A3SS,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
A3SS=transform(A3SS,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))

A3SS=transform(A3SS,IncLevel1=colsplit(IncLevel1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
A3SS=transform(A3SS,IncLevel2=colsplit(IncLevel2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
nrow(A3SS) #Nb of total splice variants

A3SS=A3SS[A3SS$FDR<0.05,] #Select the significant ones
write.csv(A3SS,file="A3SS_malevsfemale_FDR.csv") #Create the file with significant ones
nrow(A3SS) #Nb of significant variants

A3SS$IJC_sum=rowSums(A3SS$IJC_SAMPLE_1)+rowSums(A3SS$IJC_SAMPLE_2)
A3SS$SJC_sum=rowSums(A3SS$SJC_SAMPLE_1)+rowSums(A3SS$SJC_SAMPLE_2)
A3SS_count_over_20=A3SS[A3SS$IJC_sum>20&A3SS$SJC_sum>20,]
write.csv(A3SS_count_over_20,file="A3SS_malevsfemale_FDR_count20.csv")
nrow(A3SS_count_over_20) #Nb of significant variants with SJC and IJC >20

#A5SS: alternative 5'splice sites
A5SS=read.table("A5SS.MATS.JC.txt",header=TRUE)
A5SS=transform(A5SS,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
A5SS=transform(A5SS,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
A5SS=transform(A5SS,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
A5SS=transform(A5SS,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))

A5SS=transform(A5SS,IncLevel1=colsplit(IncLevel1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
A5SS=transform(A5SS,IncLevel2=colsplit(IncLevel2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
nrow(A5SS) #Nb of total splice variants

A5SS=A5SS[A5SS$FDR<0.05,] #Select the significant ones
write.csv(A5SS,file="A5SS_malevsfemale_FDR.csv") #Create the file with significant ones
nrow(A5SS) #Nb of significant variants

A5SS$IJC_sum=rowSums(A5SS$IJC_SAMPLE_1)+rowSums(A5SS$IJC_SAMPLE_2)
A5SS$SJC_sum=rowSums(A5SS$SJC_SAMPLE_1)+rowSums(A5SS$SJC_SAMPLE_2)
A5SS_count_over_20=A5SS[A5SS$IJC_sum>20&A5SS$SJC_sum>20,]
write.csv(A5SS_count_over_20,file="A5SS_malevsfemale_FDR_count20.csv")
nrow(A5SS_count_over_20) #Nb of significant variants with SJC and IJC >20

#SE: skipped exons
SE=read.table("SE.MATS.JC.txt",header=TRUE)
SE=transform(SE,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
SE=transform(SE,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
SE=transform(SE,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
SE=transform(SE,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))

SE=transform(SE,IncLevel1=colsplit(IncLevel1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
SE=transform(SE,IncLevel2=colsplit(IncLevel2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
nrow(SE) #Nb of total splice variants

SE=SE[SE$FDR<0.05,] #Select the significant ones
write.csv(SE,file="SE_malevsfemale_FDR.csv") #Create the file with significant ones
nrow(SE) #Nb of significant variants

SE$IJC_sum=rowSums(SE$IJC_SAMPLE_1)+rowSums(SE$IJC_SAMPLE_2)
SE$SJC_sum=rowSums(SE$SJC_SAMPLE_1)+rowSums(SE$SJC_SAMPLE_2)
SE_count_over_20=SE[SE$IJC_sum>20&SE$SJC_sum>20,]
write.csv(SE_count_over_20,file="SE_malevsfemale_FDR_count20.csv")
nrow(SE_count_over_20) #Nb of significant variants with SJC and IJC >20

#MXE: mutually exclusive exons
MXE=read.table("MXE.MATS.JC.txt",header=TRUE)
MXE=transform(MXE,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
MXE=transform(MXE,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
MXE=transform(MXE,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
MXE=transform(MXE,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))

MXE=transform(MXE,IncLevel1=colsplit(IncLevel1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
MXE=transform(MXE,IncLevel2=colsplit(IncLevel2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
nrow(MXE) #Nb of total splice variants

MXE=MXE[MXE$FDR<0.05,] #Select the significant ones
write.csv(MXE,file="MXE_malevsfemale_FDR.csv") #Create the file with significant ones
nrow(MXE) #Nb of significant variants

MXE$IJC_sum=rowSums(MXE$IJC_SAMPLE_1)+rowSums(MXE$IJC_SAMPLE_2)
MXE$SJC_sum=rowSums(MXE$SJC_SAMPLE_1)+rowSums(MXE$SJC_SAMPLE_2)
MXE_count_over_20=MXE[MXE$IJC_sum>20&MXE$SJC_sum>20,]
write.csv(MXE_count_over_20,file="MXE_malevsfemale_FDR_count20.csv")
nrow(MXE_count_over_20) #Nb of significant variants with SJC and IJC >20

#RI: retained intron
RI=read.table("RI.MATS.JC.txt",header=TRUE)
RI=transform(RI,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
RI=transform(RI,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
RI=transform(RI,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
RI=transform(RI,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))

RI=transform(RI,IncLevel1=colsplit(IncLevel1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
RI=transform(RI,IncLevel2=colsplit(IncLevel2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
nrow(RI) #Nb of total splice variants

RI=RI[RI$FDR<0.05,] #Select the significant ones
write.csv(RI,file="MXE_malevsfemale_FDR.csv") #Create the file with significant ones
nrow(RI) #Nb of significant variants

RI$IJC_sum=rowSums(RI$IJC_SAMPLE_1)+rowSums(RI$IJC_SAMPLE_2)
RI$SJC_sum=rowSums(RI$SJC_SAMPLE_1)+rowSums(RI$SJC_SAMPLE_2)
RI_count_over_20=RI[RI$IJC_sum>20&RI$SJC_sum>20,]
write.csv(RI_count_over_20,file="RI_malevsfemale_FDR_count20.csv")
nrow(RI_count_over_20) #Nb of significant variants with SJC and IJC >20

# Generate all DSG data together
A3SS_count_over_20$Type=rep("A3SS",nrow(A3SS_count_over_20))
A5SS_count_over_20$Type=rep("A5SS",nrow(A5SS_count_over_20))
MXE_count_over_20$Type=rep("MXE",nrow(MXE_count_over_20))
SE_count_over_20$Type=rep("SE",nrow(SE_count_over_20))
RI_count_over_20$Type=rep("RI",nrow(RI_count_over_20))

DSG=rbind(rbind(rbind(rbind(A3SS_count_over_20[,c(2,20,23,26)]
                            ,A5SS_count_over_20[,c(2,20,23,26)])
                      ,MXE_count_over_20[,c(2,22,25,28)])
                ,SE_count_over_20[,c(2,20,23,26)])
          ,RI_count_over_20[,c(2,20,23,26)])
DSG$absIncleveldifference=abs(DSG$IncLevelDifference)
write.csv(DSG,file="DS_events_malevsfemale.csv") # Total DS events between treatments

# The DS event with the largest absIncleveldifference of a DSG represent the differential splicing level of the gene.
DSG_unique_gene=DSG[order(DSG[,'GeneID'],-DSG[,'absIncleveldifference']),]
DSG_unique_gene=DSG_unique_gene[!duplicated(DSG_unique_gene$GeneID),]
write.csv(DSG_unique_gene,file="DSG_malevsfemale_unique.csv") # Total unique DS genes between treatments

#Generate all gene splicing data, including those not differentialy spliced. 
#This is used to build the scatter plot for DEG vs DSG comparison.
A3SS=read.table("A3SS.MATS.JC.txt",header=TRUE)
A3SS=transform(A3SS,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
A3SS=transform(A3SS,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
A3SS=transform(A3SS,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
A3SS=transform(A3SS,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
A3SS=transform(A3SS,IncLevel1=colsplit(IncLevel1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
A3SS=transform(A3SS,IncLevel2=colsplit(IncLevel2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
A3SS$IJC_sum=rowSums(A3SS$IJC_SAMPLE_1)+rowSums(A3SS$IJC_SAMPLE_2)
A3SS$SJC_sum=rowSums(A3SS$SJC_SAMPLE_1)+rowSums(A3SS$SJC_SAMPLE_2)
A3SS_all=A3SS[A3SS$IJC_sum>20&A3SS$SJC_sum>20,]

A5SS=read.table("A5SS.MATS.JC.txt",header=TRUE)
A5SS=transform(A5SS,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
A5SS=transform(A5SS,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
A5SS=transform(A5SS,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
A5SS=transform(A5SS,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
A5SS=transform(A5SS,IncLevel1=colsplit(IncLevel1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
A5SS=transform(A5SS,IncLevel2=colsplit(IncLevel2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
A5SS$IJC_sum=rowSums(A5SS$IJC_SAMPLE_1)+rowSums(A5SS$IJC_SAMPLE_2)
A5SS$SJC_sum=rowSums(A5SS$SJC_SAMPLE_1)+rowSums(A5SS$SJC_SAMPLE_2)
A5SS_all=A5SS[A5SS$IJC_sum>20&A5SS$SJC_sum>20,]

SE=read.table("SE.MATS.JC.txt",header=TRUE)
SE=transform(SE,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
SE=transform(SE,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
SE=transform(SE,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
SE=transform(SE,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
SE=transform(SE,IncLevel1=colsplit(IncLevel1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
SE=transform(SE,IncLevel2=colsplit(IncLevel2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
SE$IJC_sum=rowSums(SE$IJC_SAMPLE_1)+rowSums(SE$IJC_SAMPLE_2)
SE$SJC_sum=rowSums(SE$SJC_SAMPLE_1)+rowSums(SE$SJC_SAMPLE_2)
SE_all=SE[SE$IJC_sum>20&SE$SJC_sum>20,]

MXE=read.table("MXE.MATS.JC.txt",header=TRUE)
MXE=transform(MXE,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
MXE=transform(MXE,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
MXE=transform(MXE,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
MXE=transform(MXE,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
MXE=transform(MXE,IncLevel1=colsplit(IncLevel1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
MXE=transform(MXE,IncLevel2=colsplit(IncLevel2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
MXE$IJC_sum=rowSums(MXE$IJC_SAMPLE_1)+rowSums(MXE$IJC_SAMPLE_2)
MXE$SJC_sum=rowSums(MXE$SJC_SAMPLE_1)+rowSums(MXE$SJC_SAMPLE_2)
MXE_all=MXE[MXE$IJC_sum>20&MXE$SJC_sum>20,]

RI=read.table("RI.MATS.JC.txt",header=TRUE)
RI=transform(RI,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
RI=transform(RI,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
RI=transform(RI,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
RI=transform(RI,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
RI=transform(RI,IncLevel1=colsplit(IncLevel1,split=",",names=c('MB01','MB02','MB03','MB04','MB05')))
RI=transform(RI,IncLevel2=colsplit(IncLevel2,split=",",names=c('NB01','NB02','NB03','NB04','NB05')))
RI$IJC_sum=rowSums(RI$IJC_SAMPLE_1)+rowSums(RI$IJC_SAMPLE_2)
RI$SJC_sum=rowSums(RI$SJC_SAMPLE_1)+rowSums(RI$SJC_SAMPLE_2)
RI_all=RI[RI$IJC_sum>20&RI$SJC_sum>20,]

# Combine
A3SS_all$Type=rep("A3SS",nrow(A3SS_all))
A5SS_all$Type=rep("A5SS",nrow(A5SS_all))
MXE_all$Type=rep("MXE",nrow(MXE_all))
SE_all$Type=rep("SE",nrow(SE_all))
RI_all$Type=rep("RI",nrow(RI_all))

SG_all=rbind(rbind(rbind(rbind(A3SS_all[,c(2,20,23,26)]
                               ,A5SS_all[,c(2,20,23,26)])
                         ,MXE_all[,c(2,22,25,28)])
                   ,SE_all[,c(2,20,23,26)])
             ,RI_all[,c(2,20,23,26)])
SG_all$absIncleveldifference=abs(SG_all$IncLevelDifference)

# Only keep the splicing event with the max asbIncleveldiff for a unique gene list.

SG_unique_gene=SG_all[order(SG_all[,'GeneID'],-SG_all[,'absIncleveldifference']),]
SG_unique_gene=SG_unique_gene[!duplicated(SG_unique_gene$GeneID),]
write.csv(SG_unique_gene,file="SG_malevsfemale_unique.csv") #This is used for the DEG vs DSG scatter plot.
