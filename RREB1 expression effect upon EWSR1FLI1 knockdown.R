
## Effect of EWSR1-ETS (EWSR1-FLI1, EWSR1-ERG) knockdown on RREB1 expression


library(ggplot2)
library(ggpubr)


RES=read.table(file="101222 rreb1 all cell lines EWSFFLI1 ERG.txt",header=TRUE,sep="\t")


RES$Condition=as.character(RES$Condition)
RES$Fusion=as.factor(RES$Fusion)
RES$Fusion<-relevel(RES$Fusion, "EWSR1-FLI1")

RES[RES$Condition=="Control",]$Condition="Dox- (Ctrl)"
RES[RES$Condition=="Doxycycline induced EWSR1-ETS knockdown",]$Condition="Dox+ (EWSR1-ETS KD)"



p<-ggplot(RES)+
  facet_grid(. ~ Fusion, scales = "free_x", space = "free_x")+
  geom_boxplot(lwd=0.3, aes(x=cell_line,y=TC0600006855.hg.1.Value, fill=Condition))+
  theme_bw()+
  ylab("Normalized expression (log2)")+
  ggtitle("RREB1 mRNA expression with EWSR1-ETS KD")+
  theme(text=element_text(size=10), axis.title.y = element_text(size=9),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), 
        plot.title = element_text(hjust=0.5, face = "bold", size=10), 
        legend.title = element_blank(), 
        legend.text = element_text(size=9), legend.position = "bottom")+
  xlab("")+
  ylim(10.9, 14.2)+
  scale_fill_manual(values = c("brown2","cyan4"))


p



# RREB1 fold change vs EF1 KD efficiency



RREB1.table1<-read.table(file = "EF1kd_vs_RREB1_fold_change_051222.txt", sep="\t", header = TRUE)


RREB1.table1$Sig<-RREB1.table1$T_test_pval<0.05


p<-ggplot(data = RREB1.table1, mapping = aes(x=EF1_KD_efficiencyperc, y= RREB1_FC))+
  geom_point(aes(colour=factor(Sig)))+
  geom_smooth(method="lm", se=FALSE, color="dodgerblue3", size=0.5, formula= y~x)+
  ggtitle("RREB1 fold change v/s \nEWSR1-FLI1 KD efficiency")+
  xlab("EWSR1-FL1 KD efficiency %")+
  ylab("RREB1 fold change")+
  theme_bw(base_size = 9)+
  ylim(c(0.45, 1))+
  theme(plot.title = element_text(face = "bold", size=10), legend.text = element_text(size=8))

p+scale_color_brewer(palette="Dark2")+
  labs(color = "p-value<0.05:")