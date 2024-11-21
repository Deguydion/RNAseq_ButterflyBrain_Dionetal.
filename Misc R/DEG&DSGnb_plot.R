#BUILD NB OF DSG AND DEG PLOT IN R

library(ggplot2) #plot building
library(wesanderson) #plot colours

plot<-ggplot(genenb,aes(x=comparison,y=number,fill=genetype)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values=wes_palette("Darjeeling1", n = 2)) +
  ylab("Gene number") +
  xlab("") + 
  labs(fill="") + 
  scale_x_discrete(limits = c("MvsN", "NvsWt", "WtvsNB1" , "WtvsNB2", "NB1vsNB2")) +
  theme(axis.title.y = element_text(face="bold",size=12)) + 
  theme(legend.title=element_text(face="bold",size=12)) +
  theme(legend.text = element_text(color="black",size=12)) +
  theme(axis.line=element_line(color="black", linewidth=1)) +
  theme(axis.text=element_text(color="black",size=12)) +
  theme(axis.ticks=element_line(linewidth = 1)) +
  theme(axis.ticks.length =unit(0.2, "cm")) +
  theme(panel.background=element_rect(fill="white"))
plot
ggsave("dsgxdeg_nb.png",plot=plot())
