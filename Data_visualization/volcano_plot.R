#Join Gene names

data.joined <- left_join(Hela_data, Gene_names, c("ID"="entry"))
View(data.joined)

Hela_df <- data.joined[,c(10,1:5,7)]
colnames(Hela_df) <- c("gene.names","gene.names_ID","ID","pval","padj","significant","log2FC")

Hela_df <- make_unique(Hela_df, "gene.names", "ID", delim = ";")
Hela_df <- Hela_df[,c(8,3:7)]
colnames(Hela_df)[1]<- "gene.names"

##### VOLCANO PLOT #########
# Generate logical column 
Hela_df <- Hela_df %>% mutate(threshold_padj = padj < 0.05,threshold_fc = abs(log2FC) >= 1.5)

#Order by padj
Hela_df <- Hela_df[order(Hela_df$padj),]

# Create volcano plot
ggplot(Hela_df, aes(x = log2FC, y = -log10(padj), color = interaction(threshold_padj, threshold_fc), label = gene.names)) + 
  geom_point() +
  geom_label_repel(data = head(Hela_df,40), aes(label=gene.names)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")+
  geom_vline(xintercept = 1.5, linetype="dashed", color = "red")+
  geom_vline(xintercept = -1.5, linetype="dashed", color = "red")+
  scale_x_continuous(limits = c(-7, 7), breaks = seq(from = -7.5, to = 7.5, by = 1.5))+
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  theme_bw()+
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25)))
