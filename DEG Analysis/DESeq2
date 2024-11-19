#DEG ANALYSIS WITH R WITHIN RSTUDIO USING DESEQ2 PACKAGE
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(RColorBrewer)
library(vsn)

#DATA INPUT AND FILTERING
#Data import
setwd("C:/Users/folders/DEG")
brain_countData <-as.matrix(read.csv("brain_gene_count_matrix.csv",row.names="gene_id"))
brain_colData<-read.csv("brain_pheno_data.csv",row.names=1)

#Check that all sample IDs in colData are also in CountData and match their orders
all(rownames(brain_colData) %in% colnames(brain_countData))
brain_countData <- brain_countData[, rownames(brain_colData)]
all(rownames(brain_colData) == colnames(brain_countData))

#Create a DESeqDataSet from count matrix and labels
brain_dds<-DESeqDataSetFromMatrix(countData=brain_countData,colData=brain_colData,design=~treatment)

#Pre-filtering of low count genes. Rows with less than 10 reads in total are removed
brain_keep<-rowSums(counts(brain_dds)) >= 10
brain_dds<-brain_dds[brain_keep,]

#Generate normalized counts (this is for extra analysis, during DE normalization is done in the glm)
brain_dds<-estimateSizeFactors(brain_dds)
sizeFactors(brain_dds) #normalization factors of each sample 
normalized_brain_counts<-counts(brain_dds,normalized=TRUE) #retrieve the normalized counts matrix
write.table(normalized_brain_counts,file="brain_normalized_counts.txt",sep="\t",quote=F,col.names=NA)

colSums(counts(brain_dds)) # Total number of raw counts per sample
colSums(counts(brain_dds,normalized=T))# Total number of normalized counts per sample

#Data transformation
brain_rld<-rlog(brain_dds,blind=TRUE)
head(assay(brain_rld),3)

#Principal component analysis PCA 
plotPCA(brain_rld, intgroup="treatment")
plotPCA(brain_rld, intgroup="treatment") + geom_text(aes(label=name))

#Heatmap of the count matrix
pcaData<-plotPCA(brain_rld,intgroup="treatment",ntop=500,returnData=TRUE)
percentVar<-round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData,aes(PC1, PC2, color=treatment)) +
  geom_point(size=5) +
  scale_color_manual(values = c("#BBBBBB","#000000","#009988","#0177BB","#EE7733")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  theme(text=element_text(size=20))
ggsave(pcaData,filename='PCA.eps')

#Plot of individual gene counts (here we use MSTRG.18812 as an example)
genecount<-plotCounts(brain_dds,gene="MSTRG.18812",intgroup="treatment",returnData=TRUE)
ggplot(genecount,aes(x=treatment, y=count)) +
  geom_point(aes(color=treatment),position=position_jitter(w=0.1,h=0),size=2.5) +
  scale_x_discrete(limits=c("male", "naive_fem", "wt_exp_fem","nb1_exp_fem","nb2_exp_fem")) +
  scale_color_manual(values=c("darkgrey","black","#08B59C","#0796DD","#F49867"))+
  theme(panel.background= element_rect(fill = "white", colour = "black"),panel.grid.major=element_line(size=0.5,linetype ='dashed',colour="grey"))+
  theme(text=element_text(size=14,face="bold"))+
  labs(subtitle="MSTRG.18812")
ggsave(genecount,filename='MSTRG.18812_counts.eps')

#DESEQ ANALYSIS ANS SHRINKAGE FOR VIZUALIZATION AND RANKING OF GENES

#Males versus naive females
brain_dds <- DESeq(brain_dds)
resultsNames(brain_dds) # to see what names to use
brain_res_mvsnaive <- results(brain_dds, contrast=c("treatment","naive_fem","male"),alpha=0.05) #Males=base, results are changes in females
summary(brain_res_mvsnaive) #Print a summary of the contrast analyse
brain_res_mvsnaive <- brain_res_mvsnaive[order(brain_res_mvsnaive$pvalue),]#Re-ordering the result datasets by the smallest p-values
write.csv(as.data.frame(brain_res_mvsnaive),file="brain_mvsnaive_DEG.csv")#Exporting the results to csv files
resLFC_brain_NvsM<-lfcShrink(brain_dds,coef = "treatment_naive_fem_vs_male",type="apeglm")
resLFC_brain_NvsM<-resLFC_brain_NvsM[order(resLFC_brain_NvsM$pvalue),]
write.csv(as.data.frame(resLFC_brain_NvsM),file="resLFC_brain_mvsnaive_DEG.csv")

#Naive females versus Wt-exp
brain_dds$treatment<-relevel(brain_dds$treatment,ref="naive_fem") #Adjust factor level for shrinkage purposes, here naive vs wt
brain_dds <- DESeq(brain_dds)
resultsNames(brain_dds)
brain_res_naivevswt <- results(brain_dds,contrast=c("treatment","wt_exp_fem","naive_fem"),alpha=0.05) #changes in wt females
summary(brain_res_naivevswt)
brain_res_naivevswt <- brain_res_naivevswt[order(brain_res_naivevswt$pvalue),]
write.csv(as.data.frame(brain_res_naivevswt),file="brain_naivevswt_DEG.csv")
resLFC_brain_NvsWt<-lfcShrink(brain_dds,coef = "treatment_wt_exp_fem_vs_naive_fem",type="apeglm")
resLFC_brain_NvsWt<-resLFC_brain_NvsWt[order(resLFC_brain_NvsWt$pvalue),]
write.csv(as.data.frame(resLFC_brain_NvsWt),file="resLFC_brain_nvswt_DEG.csv")

#Wt-exp versus nb1- and nb2-exposed 
brain_dds$treatment<-relevel(brain_dds$treatment,ref="wt_exp_fem")
brain_dds <- DESeq(brain_dds)
resultsNames(brain_dds)
brain_res_wtvsnb1 <- results(brain_dds, contrast=c("treatment","nb1_exp_fem","wt_exp_fem"),alpha=0.05) #changes in NB1
brain_res_wtvsnb2 <- results(brain_dds, contrast=c("treatment","nb2_exp_fem","wt_exp_fem"),alpha=0.05) #changes in NB2
summary(brain_res_wtvsnb1)
summary(brain_res_wtvsnb2)
brain_res_wtvsnb1 <- brain_res_wtvsnb1[order(brain_res_wtvsnb1$pvalue),]
brain_res_wtvsnb2 <- brain_res_wtvsnb2[order(brain_res_wtvsnb2$pvalue),]
write.csv(as.data.frame(brain_res_wtvsnb2),file="brain_wtvsnb2_DEG.csv")
write.csv(as.data.frame(brain_res_wtvsnb1),file="brain_wtvsnb1_DEG.csv")
resLFC_brain_wtvsnb1<-lfcShrink(brain_dds,coef = "treatment_nb1_exp_fem_vs_wt_exp_fem",type="apeglm")
resLFC_brain_wtvsnb2<-lfcShrink(brain_dds,coef = "treatment_nb2_exp_fem_vs_wt_exp_fem",type="apeglm")
resLFC_brain_wtvsnb1 <- resLFC_brain_wtvsnb1[order(resLFC_brain_wtvsnb1$pvalue),]
resLFC_brain_wtvsnb2 <- resLFC_brain_wtvsnb2[order(resLFC_brain_wtvsnb2$pvalue),]
write.csv(as.data.frame(resLFC_brain_wtvsnb1),file="resLFC_brain_wtvsnb1_DEG.csv")
write.csv(as.data.frame(resLFC_brain_wtvsnb2),file="resLFC_brain_wtvsnb2_DEG.csv")

#nb1-exp versus nb2-exposed females
brain_dds$treatment<-relevel(brain_dds$treatment,ref="nb1_exp_fem")
brain_dds <- DESeq(brain_dds)
resultsNames(brain_dds)
brain_res_nb1vsnb2 <- results(brain_dds, contrast=c("treatment","nb1_exp_fem","nb2_exp_fem"),alpha=0.05) #Changes in NB1
summary(brain_res_nb1vsnb2)
brain_res_nb1vsnb2 <- brain_res_nb1vsnb2[order(brain_res_nb1vsnb2$pvalue),]
write.csv(as.data.frame(brain_res_nb1vsnb2),file="brain_nb1vsnb2_DEG.csv")
resLFC_brain_nb1vsnb2<-lfcShrink(brain_dds,coef = "treatment_nb2_exp_fem_vs_nb1_exp_fem",type="apeglm")
write.csv(as.data.frame(resLFC_brain_nb1vsnb2),file="resLFC_brain_nb1vsnb2_DEG.csv")


#BUILDING DEG PLOTS

#MA plots
plotMA(resLFC_brain_NvsM,ylim=c(-4,3),cex=0.75,cex.lab=1.5,cex.axis=1.5)
plotMA(resLFC_brain_NvsWt,ylim=c(-2.5,5),cex=0.75,cex.lab=1.5,cex.axis=1.5)
plotMA(resLFC_brain_wtvsnb1,ylim=c(-2,5),cex=0.75,cex.lab=1.5,cex.axis=1.5)
plotMA(resLFC_brain_wtvsnb2,ylim=c(-1,4),cex=0.75,cex.lab=1.5,cex.axis=1.5)
plotMA(resLFC_brain_nb1vsnb2,ylim=c(-1,6),cex=0.75,cex.lab=1.5,cex.axis=1.5)


#Volcanoe plots  
#To remove all the labels, use NA for select lab
#Example shown for MvsN, same is done for other comparisons.
EnhancedVolcano(res_brain_NvsM,
                lab = rownames(resLFC_brain_NvsM),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = c(NA),
                xlab = bquote(~Log[2]~ 'fold change'),
                xlim = c(-4, 4),
                ylim = c(0, -log10(10e-7)),
                pCutoff = 50e-3,
                FCcutoff = 0.0,
                pointSize = 2,
                labSize = 2.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = FALSE,
                colAlpha = 4/5,
                legendPosition = 'none',
                legendLabSize = 14,
                legendIconSize = 4.0,
                typeConnectors ='closed',
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                col=c('black','grey','grey','black'),
                subtitle = NULL,
                title=NULL)

