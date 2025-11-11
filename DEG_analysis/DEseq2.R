#DEG ANALYSIS WITH R WITHIN RSTUDIO USING DESEQ2 PACKAGE
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(RColorBrewer)
library(vsn)
library(SmartSVA)
library(ggrepel)

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

#Generate variance-stabilize counts and generate mnatrixes for smartSVA
vsd<-vst(brain_dds,blind=TRUE)
exprs<-assay(vsd)   
class(exprs)
mod<-model.matrix(~treatment,colData(vsd)) #full model
mod0<-model.matrix(~1,colData(vsd)) #null model

#Estimate surrogate variables. We find none, and proceed to DEG analysis without extra variable
smart_sv<-smartsva.cpp(dat=exprs,mod=mod,mod0=mod0,n.sv=NULL)
n_sv<-ncol(smart_sv$sv)
cat("Number of surrogate variables detected:",n_sv,"\n")
head(smart_sv$sv)

#Unbiased principal component analysis PCA on 500 random genes
plotPCA(brain_rld, intgroup="treatment")
plotPCA(brain_rld, intgroup="treatment") + geom_text(aes(label=name))

#Unbiased heatmap of the count matrix
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


#PERFORMING LRT ANALYSIS ON EXPOSED FEMALES (AS PER REVIEWER SUGGESTION, JUNE 2025) 

brain_countData <- as.matrix(read.csv("brain_gene_count_matrix.csv", row.names = "gene_id"))
brain_colData <- read.csv("brain_pheno_data.csv", row.names = 1)

#Subset to only samples with the three exposure treatments
keep_samples <- brain_colData$treatment %in% c("wt_exp_fem", "nb1_exp_fem", "nb2_exp_fem")
counts_sub <- brain_countData[, keep_samples]
colData_sub <- brain_colData[keep_samples, , drop = FALSE]
colData_sub$treatment <- factor(colData_sub$treatment, levels = c("wt_exp_fem", "nb1_exp_fem", "nb2_exp_fem")) #Ensure 'treatment' is a factor 

#DESeq2 analysis
dds_lrt <- DESeqDataSetFromMatrix(countData = counts_sub, colData = colData_sub, design = ~ treatment)
dds_lrt <- DESeq(dds_lrt, test = "LRT", reduced = ~ 1)
res_lrt <- results(dds_lrt)
res_lrt_ordered <- res_lrt[order(res_lrt$padj), ] #Orders by adjusted p-value
head(res_lrt_ordered)
write.csv(as.data.frame(sig_res_lrt), file = "LRT_DEGs_wt_nb1_nb2.csv")

#Do the pairwise with wald on the DEG from the LRT
sig_genes <- rownames(sig_res_lrt)
dds_wald <- DESeqDataSetFromMatrix(countData = counts_sub[sig_genes, ], colData = colData_sub, design = ~ treatment)
dds_wald <- DESeq(dds_wald)
res_nb1_vs_wt <- results(dds_wald, contrast = c("treatment", "nb1_exp_fem", "wt_exp_fem"))
res_nb2_vs_wt <- results(dds_wald, contrast = c("treatment", "nb2_exp_fem", "wt_exp_fem"))
res_nb2_vs_nb1 <- results(dds_wald, contrast = c("treatment", "nb2_exp_fem", "nb1_exp_fem"))

sig_nb1_vs_wt <- res_nb1_vs_wt[which(res_nb1_vs_wt$padj < 0.05), ]
sig_nb2_vs_wt <- res_nb2_vs_wt[which(res_nb2_vs_wt$padj < 0.05), ]
sig_nb2_vs_nb1 <- res_nb2_vs_nb1[which(res_nb2_vs_nb1$padj < 0.05), ]

write.csv(as.data.frame(sig_nb1_vs_wt), file = "DEGs_nb1_vs_wt.csv")
write.csv(as.data.frame(sig_nb2_vs_wt), file = "DEGs_nb2_vs_wt.csv")
write.csv(as.data.frame(sig_nb2_vs_nb1), file = "DEGs_nb2_vs_nb1.csv")


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


#BUILD THE HEATMAPS FROM THE COUNT MATRIX: UNBIASED MALE VS FEMALE; UNBIASED NAIVE VERSUS WT-EXP FEMALE AMD EXPOSED FEMALES
#AS PER REVIEWER'S SUGGESTION, JULY 2025

library(pheatmap)
library(RColorBrewer)

# Define custom colors for annotation
treatment_colors <- c(male = "lightgrey", naive_fem = "black", wt_exp_fem = "#EE7733", nb1_exp_fem = "#009988", nb2_exp_fem = "#0177BB")
heat_colors <- colorRampPalette(c("purple", "purple", "black", "yellow", "yellow"))(100)

#Naive individuals
subset1 <- brain_colData[brain_colData$treatment %in% c("male", "naive_fem"), ]
dds1 <- DESeqDataSetFromMatrix(countData = brain_countData[, rownames(subset1)],colData = subset1, design = ~ treatment)
brain_keep<-rowSums(counts(dds1)) >= 10
dds1<-dds1[brain_keep,]

rld1 <- rlog(dds1, blind = TRUE)
mat1 <- assay(rld1)
mat1_scaled <- t(scale(t(mat1)))
annotation1 <- data.frame(Treatment = subset1$treatment)
rownames(annotation1) <- rownames(subset1)
treatment_colors_subset1 <- treatment_colors[levels(droplevels(subset1$treatment))]

pheatmap(
  mat1_scaled[sample(rownames(mat1_scaled), 500), ],
  color = heat_colors,
  annotation_col = annotation1,
  annotation_colors = list(Treatment = treatment_colors_subset1),
  show_rownames = FALSE,
  main = "Unbiased Heatmap: Male vs Naive",
  filename = "unbiased_male_vs_naive_heatmap.png"
)

#Exposed females
subset2 <- brain_colData[brain_colData$treatment %in% c("wt_exp_fem","nb1_exp_fem","nb2_exp_fem"), ]
dds2 <- DESeqDataSetFromMatrix(countData = brain_countData[, rownames(subset2)], colData = subset2, design = ~ treatment)
brain_keep<-rowSums(counts(dds2)) >= 10
dds2<-dds2[brain_keep,]

rld2 <- rlog(dds2, blind = TRUE)
mat2 <- assay(rld2)
mat2_scaled <- t(scale(t(mat2)))
annotation2 <- data.frame(Treatment = subset2$treatment)
rownames(annotation2) <- rownames(subset2)
treatment_colors_subset2 <- treatment_colors[levels(droplevels(subset2$treatment))]

pheatmap(
  mat2_scaled[sample(rownames(mat2_scaled), 500), ],
  color = heat_colors,
  annotation_col = annotation2,
  annotation_colors = list(Treatment = treatment_colors_subset2),
  show_rownames = FALSE,
  main = "Unbiased Heatmap: Exposed females",
  filename = "unbiased_expfemales_heatmap.png"
  )

#Naives versus Wt females
subset3 <- brain_colData[brain_colData$treatment %in% c("naive_fem", "wt_exp_fem"), ]
dds3 <- DESeqDataSetFromMatrix(countData = brain_countData[, rownames(subset3)], colData = subset3, design = ~ treatment)

brain_keep<-rowSums(counts(dds3)) >= 10
dds3<-dds3[brain_keep,]

rld3 <- rlog(dds3, blind = TRUE)
mat3 <- assay(rld3)
mat3_scaled <- t(scale(t(mat3)))
annotation3 <- data.frame(Treatment = subset3$treatment)
rownames(annotation3) <- rownames(subset3)
treatment_colors_subset3 <- treatment_colors[levels(droplevels(subset3$treatment))]

pheatmap(
  mat3_scaled[sample(rownames(mat3_scaled), 500), ],
  color = heat_colors,
  annotation_col = annotation3,
  annotation_colors = list(Treatment = treatment_colors_subset3),
  show_rownames = FALSE,
  main = "Unbiased Heatmap: Naives and wt-exposed females",
  filename = "unbiased_naive&wt-expfem_heatmap.png"
)

#BUILDING THE PCA AND HEATMAPS FROM the DEGS. WE CHOSE TO USE DEGS AT P<0.1 BECAUSE WE FIND VERY FEW AT P<0.05
#AS PER REVIEWER'S REQUIREMENT NOV2025

# Define custom colors for annotation
treatment_colors<-c(
  male="lightgrey",
  naive_fem="black",
  wt_exp_fem="#EE7733",
  nb1_exp_fem="#009988",
  nb2_exp_fem="#0177BB"
)

# Define heatmap color palette
heat_colors<-colorRampPalette(c("purple","purple","black","yellow","yellow"))(100)

#Run DESeq2, subset to DEG at p<0.1 and keep first 500
dds<-DESeq(brain_dds)
res<-results(dds,alpha=0.1)
deg<-subset(res,padj<0.1 & !is.na(padj))
deg<-deg[order(deg$padj), ]
top500_genes<-rownames(head(deg,500))

#Transform data for stabilizing variance
vsd<-vst(dds,blind=FALSE)
vsd_top500<-assay(vsd)[top500_genes, ]
vsd_top500_scaled<-t(scale(t(vsd_top500))) #scale genes to zscore for better visual contrast for heatmap

#Preparing PCA on DEGs
pca<-prcomp(t(vsd_top500_scaled),scale.=TRUE)
nrow(vsd_top500)  #provides total number of 0.1DEGs genes used in PCA

pcaData<-as.data.frame(pca$x[, 1:2])  #keeps only PC1 and PC2
pcaData$name<-colnames(vsd_top500_scaled) 
pcaData$treatment<-colData(vsd)$treatment  
percentVar<-round(100 *(pca$sdev^2/sum(pca$sdev^2))[1:2], 1) #Calculates variance explained by PCs

#Building PCA plots 
xlim_val<-max(abs(pcaData$PC1)) * 1.1
ylim_val<-max(abs(pcaData$PC2)) * 1.1

ggplot(pcaData,aes(x=PC1,y=PC2,color=treatment)) +
  geom_point(size=4) +
  scale_color_manual(values=treatment_colors) +
  scale_x_continuous(limits=c(-xlim_val, xlim_val),expand=c(0, 0)) +
  scale_y_continuous(limits=c(-ylim_val, ylim_val),expand=c(0, 0)) +
  labs(
    title="All data: Top 139 DEGs at adj.p<0.1",   
    x=paste0("PC1:",percentVar[1],"% variance"),
    y=paste0("PC2:",percentVar[2],"% variance")
  ) +
  coord_fixed() +
  theme_bw(base_size=14) +
  theme(
    axis.text.x=element_text(face="bold"),  
    axis.text.y=element_text(face="bold"),
    panel.background=element_rect(fill ="grey97",color=NA),  #design for same colors as the unbiased PCA
    panel.grid.major=element_line(color="white"),              
    panel.grid.minor=element_line(color="white"),              
    panel.border=element_blank(),                                 
    axis.line=element_line(color="white"),                      
    legend.title=element_blank(),
    plot.title=element_text(hjust=0.5,size=16,face="bold")
  )

#Subsetting data by treatment comparisons for heatmaps
#Males versus naive females. Same procedure is used for the other comparisons. 
subset1<-brain_colData[brain_colData$treatment %in% c("male","naive_fem"), ]
count_subset1<-brain_countData[,rownames(subset1)]
dds1<-DESeqDataSetFromMatrix(countData=count_subset1,colData=subset1,design=~treatment)
dds1 <- dds1[rowSums(counts(dds1)) >= 10, ]
dds1 <- DESeq(dds1)

res1<-results(dds1 alpha=0.1) 
deg1<-subset(res1,padj<0.1 & !is.na(padj))
nrow(deg1) #how much genes passed the padj
deg1<-deg1[order(deg1$padj), ] #Select top 500 DEGs
topN<-min(500,nrow(deg1))
top500_genes1<-rownames(head(deg1,topN))
rld1<-rlog(dds1,blind=TRUE) #Normalize data
mat1<-assay(rld1)[top500_genes1, ]
mat1_scaled<-t(scale(t(mat1)))

annotation1<-data.frame(Treatment=subset1$treatment) #Create annotation
rownames(annotation1)<-rownames(subset1)
treatment_colors_subset1<-treatment_colors[levels(droplevels(subset1$treatment))] #Restrict treatment_colors to only this subset

#Plot heatmap
pheatmap(
  mat1_scaled,
  color = heat_colors,
  annotation_col = annotation1,
  annotation_colors = list(Treatment = treatment_colors_subset1),
  show_rownames = FALSE,
  main = "Top 31 DEGs at adj. p <0.1: Males vs Naive Females"
)

