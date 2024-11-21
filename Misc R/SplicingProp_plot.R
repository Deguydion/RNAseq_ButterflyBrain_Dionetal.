#BUILD PROPORTIONS OF SPLICING TYPES PLOT IN R

library(ggplot2) #to build the plot
library(wesanderson) #for the colors

dsgprop<-read.table("DSG_prop.txt",header=T)
head(dsgprop)

plot<-ggplot(dsgprop, aes(x=comparison,y=number,fill=type)) +
  geom_col(colour="black", position="fill") + 
  scale_y_continuous(labels=scales::percent) + 
  scale_fill_manual(values=wes_palette("Darjeeling1", n = 5)) +
  ylab("% Splicing type") +
  xlab("") + 
  labs(fill="Splicing type") + 
  theme(axis.title.y = element_text(face="bold",size=12)) + 
  theme(legend.title=element_text(face="bold",size=12)) +
  theme(legend.text = element_text(color="black",size=12)) +
  theme(axis.line=element_line(color="black", linewidth=1)) +
  theme(axis.text=element_text(color="black",size=12)) +
  theme(axis.ticks=element_line(linewidth = 1)) +
  theme(axis.ticks.length =unit(0.2, "cm")) +
  theme(panel.background=element_rect(fill="white"))
plot
ggsave("DSGprop.png",plot=plot())
