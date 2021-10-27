
libraries <- c('RColorBrewer', 'pheatmap', 'tidyverse','scales','readxl','ggthemes', 'ggrepel')
lapply(libraries,library, character.only = TRUE)

setwd('C:/*')

EMT_4<-read_excel("Resultados para ggplot.xlsx", sheet=4)
EMT_COR <- EMT_4[,2:3]
rownames(EMT_COR)<- EMT_4$Genes

####Correlation proteome-qPCR######

cor(EMT_COR$qPCR, EMT_COR$Proteome, method = "pearson")

ggplot(EMT_COR, aes(x= qPCR, y=Proteome)) +
  geom_point(colour="#1F77B4", size=3)+
  xlab("Log2FC qPCR")+
  ylab("Log2FC Proteome")+
  ylim(-5,5)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_vline(xintercept =0, linetype="dashed", color = "red")+
  geom_label_repel(data = EMT_COR, aes(label=rownames(EMT_COR)),colour ="#1F77B4", size= 5)+
  geom_smooth(method=lm, colour = "#1F77B4")+
  stat_cor(method = "pearson")+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(linetype="dashed"),
        panel.grid.minor = element_line(linetype="dashed"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title=element_text(size=20,face="bold"),
        panel.border = element_blank())
