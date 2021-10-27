
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
