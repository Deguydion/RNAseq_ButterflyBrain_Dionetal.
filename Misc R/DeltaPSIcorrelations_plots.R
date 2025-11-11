### DELTAPSI CORRELATIONS ACROSS TREATMENT COMPARISONS ###
### FROM REVIEWER COMMENT, NOV25 ###

#For this we create a data table with deltaPSI of genes with at least 1 significant deltaPSI across treatment comparisons.
#For easier vizualisation, deltaPSI values are absolute, and MSP2 changes are illustrated with different colors
#Frame has 5 columns: GeneID, DeltaPSI_WtvsNB1, DeltaPSI_WtvsNB2, DeltaPSI_NB1vsNB2 and direction of MSP2 change (higher or lower). 
#Data are in Additional File 3 Table 7

#Load libraries and data
library(ggplot2)
library(ggpubr)
cor_data<-read.table("path/to/correlation_data.txt",header=TRUE,sep ="\t",stringsAsFactors=FALSE)
head(cor_data)
cor_data$Direction<-factor(cor_data$Direction,levels=c("higher","lower"))

#Scatterplot WTvsNB2 vs NB1vsNB2
#Same method used for the two other comparisons.
ggscatter(cor_data,
          x="PSI_Wt_NB2",
          y="PSI_NB1_NB2",
          color="Direction",                  
          palette=c("blue","orange"),        
          size=4,                             
          add=NULL,                           
          conf.int=FALSE,                     
          cor.coef=TRUE,                      
          cor.method="pearson",               
          cor.coef.size=6,                    
          xlab="ΔPSI WT vs NB2",
          ylab="ΔPSI NB1 vs NB2",
          title="ΔPSI correlation") +
  geom_smooth(method="lm",color="black",se=FALSE) +  #linear regression fitted by ordinary least squares
  labs(color="Changes in MSP2 amounts") + 
  theme_bw(base_size=14) +
  theme(
    axis.title=element_text(face="bold",size=16),
    axis.text=element_text(face="bold",size=14),
    plot.title=element_text(face="bold",size=18,hjust=0.5),
    legend.title=element_text(face="bold",size=14),
    legend.text=element_text(face="bold",size=13),
    legend.position=c(0.80,0.2),     
    legend.background=element_blank(), 
    legend.key=element_blank(),        
    panel.border=element_rect(color="black",size = 1)
  )

