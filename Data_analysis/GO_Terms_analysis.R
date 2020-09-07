####ClusterPorfile####

BiocManager::install("clusterProfiler")
install.packages('GOplot')
library(org.Hs.eg.db)
library(clusterProfiler)
library(GOplot)

#First dataframe with the list of genes, logFC and p-value adj.
Hela_df_clusPro1 <- Hela_df[,c(1,3:4,7)]
colnames(Hela_df_clusPro1) <- c("ID","P.Value","adj.P.Val","logFC")

#Second dataframe with GO enrichment analysis.
genes <- Hela_df[,1]
genes_logFC <- Hela_df[,c(1,7)]
colnames(genes_logFC)<- c("ID", "logFC")
genes_logFC$logFC <- gsub(",",".",genes_logFC$logFC)
genes_logFC$logFC <- as.numeric(genes_logFC$logFC)
genes.df <- bitr(genes, fromType = "SYMBOL",
                 toType = c("ENSEMBL", "ENTREZID"),
                 OrgDb = org.Hs.eg.db)

ego_CC <- enrichGO(gene          = genes.df$ENTREZID,
                keyType       = "ENTREZID",
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

ego_BP <- enrichGO(gene          = genes.df$ENTREZID,
                   keyType       = "ENTREZID",
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

ego_MF <- enrichGO (gene          = genes.df$ENTREZID,
                   keyType       = "ENTREZID",
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

ego_BP@result <- mutate(ego_BP@result, Category = "BP")
ego_MF@result <- mutate(ego_MF@result, Category = "MF")
ego_CC@result <- mutate(ego_CC@result, Category = "CC")

ego <- rbind(ego_BP@result, ego_CC@result) %>%
  rbind(ego_MF@result)

FA.results <- as.data.frame(ego)
rownames(FA.results) <- NULL
FA.results <- FA.results[,c(10,1:9)]
colnames(FA.results)[c(3,7,9)]<- c("Terms","adj_pval", "genes")
FA.results2 <- FA.results[,c(1:4,7,9)]
FA.results2$genes <- gsub('/',',',FA.results2$genes)
