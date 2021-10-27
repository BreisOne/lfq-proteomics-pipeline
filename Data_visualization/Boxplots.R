
libraries <- c('RColorBrewer', 'pheatmap', 'tidyverse','scales','readxl','ggthemes', 'ggrepel')
lapply(libraries,library, character.only = TRUE)

setwd('C:/Users/Brise/OneDrive - Universidade de Vigo/Tesis - Universidad de Vigo/KO BJ-5TA/KO ALMS1/EMT biomarkers')

EMT_24H <-read_excel("Resultados para ggplot.xlsx", sheet=1)
EMT_48H <-read_excel("Resultados para ggplot.xlsx", sheet=2)

EMT_24H <- data.frame(EMT_24H)
EMT_48H <- data.frame(EMT_48H)

#####Boxplot with Rgpraph######


##Paletta de colores de ggplot para boxplot##

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Then I make the boxplot, asking to use the 2 factors : variety (in the good order) AND genotype:

colors = gg_color_hue(2) 

myplot <- boxplot(FC ~ Genotype+Gene , data = EMT_48H,
          boxwex=0.4 , ylab="Fold Change",
          main="TGF-B 48H BJ-5TA",
          ylim = c(0, 3),
          col= c(colors[2],colors[1]),  
          xaxt="n",
          frame.plot = FALSE)

# To add the label of x axis
my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
my_names <- my_names[seq(1 , length(my_names) , 2)]
axis(side = 1, 
     at = seq(1.5 , 14 , 2), 
     labels = my_names, cex=0.3)

# Add the grey vertical lines
for(i in seq(0.5 , 20 , 2)){ 
  abline(v=i,lty=15, col="grey")
}

# Add a legend
legend("topright",                   # Create legend outside of plot
       legend = c("ALMS1 +/+", "ALMS1 -/-"),
       col=c(colors[2],colors[1]),
       pch = 15,
       pt.cex = 2,
       inset = c(0.035, 0))

# legend("bottomright", legend = c("ALMS +/+", "ALMS -/-"), 
#        col=c(colors[2],colors[1]),
#        pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = c(0, 0))


#####Boxplot with ggplot#####

ggplot(EMT_24H, aes(x=Gene, y= FC, fill= Genotype))+
  geom_boxplot()+
  ylim(0,2)+
  labs(title = "TGF-B 24H" , x = NULL, y = "Fold Change" , fill = "Genotype")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(EMT_48H, aes(x=Gene, y= FC, fill=Genotype))+
  geom_boxplot()+
  ylim(0,2.5)+
  labs(title = "TGF-B 48H" , x = NULL, y = "Fold Change" , fill = "Genotype")

ggplot(EMT_72H, aes(x=Gene, y= FC, fill=Genotype))+
  geom_boxplot()+
  ylim(0,2.5)+
  labs(title = "TGF-B 72H" , x = NULL, y = "Fold Change" , fill = "Genotype")

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

######Heatmap BJ and Hela######

EMT_24H_Heatmap <-read_excel("Resultados para ggplot.xlsx", sheet=3)
EMT_24H_Heatmap_2 <- EMT_24H_Heatmap[1:7,]
EMT_24H_Heatmap_3 <- EMT_24H_Heatmap_2[,2:7]
row.names(EMT_24H_Heatmap_3) <- EMT_24H_Heatmap_2$...1

heat_colors <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)

pheatmap(EMT_24H_Heatmap_3,
         color = heat_colors,
         cluster_rows = FALSE,
         #cluster_cols = mat_cluster_cols,
         show_rownames = TRUE,
         #legend_breaks = c(-1.9,-1,0,1,1.9),
         #annotation_col = select(BJ_Metadata, replicates, condition),
         #annotation_colors = colours_heatmap_annotation,
         scale = "column",
         main = "qPCR BJ-5TA and Hela TGF-B 24h"
)
