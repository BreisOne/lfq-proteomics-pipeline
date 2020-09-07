##########DEP: Differential Enrichment analysis of Proteomics data########
#generate the data frame of the experimental design

install.packages("BiocManager")
BiocManager::install("DEP")
library(DEP)

# Check if there are duplicates in the gene names
data.joined$gene.names %>% duplicated() %>% any()

# Evaluate exactly which gene names are duplicated
data.joined %>% group_by(gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)

#Generate unique identifiers for each Gene / Protein
data_unique <- make_unique(data.joined, "gene.names", "protein.IDs", delim = ";")

#Validate that all names are unique
data_unique$name %>% duplicated() %>% any()

# Generate the dataframe of the experimental design
l <- colnames(data_unique)[6:11]
c <- c("Control","Control","Control", "ALMS1KO","ALMS1KO","ALMS1KO")
r <- c("1","2","3","1","2","3")
colnames <- c("label", "condition","replicate")
exp_design <- data.frame(l,c,r, stringsAsFactors = FALSE)
colnames(exp_design) <- colnames

#Generate a SummarizedExperiment object from the experimental design
exp_Cols <- c(6:11)
data_se <- make_se(data_unique,exp_Cols, exp_design)

####PLOTS FOR DATA EVALUATION####

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr= 1)

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

# Normalize the data
data_norm <- normalize_vsn(data_filt)

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)


#######IMPUTE MISSING DATA USING RANDOM DRAWS FROM A GAUSSIAN DISTRIBUTION CENTERED AROUND A MINIMAL VALUE (FOR MNAR)########
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)

########DIFFERENTIAL ENRICHMENT ANALYSIS  BASED ON LINEAR MODELS AND EMPHERICAL BAYES STATISTICS########

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "Control")

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = 1.5)

# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = TRUE,
             indicate = c("condition", "replicate"))

# Plot a volcano plot for the contrast "ALMS1KO_vs_Control""
plot_volcano(dep, contrast = "ALMS1KO_vs_Control", label_size = 2.5, add_names = TRUE)

# Plot a barplot for USP15 and IKBKG
plot_single(dep, proteins = c("ALDH16A1","S100A6", "CRABP2", "EDIL3", "DCK", "DSP", "VIM","EEF1A2", "FKBP10", "IDH2",
                              "ITGB5", "LMNB1", "LMNB2", "NAPRT","O95379-4", "P0DMM9-3", "P98175-5", "PHAX", "PSIP1", 
                              "Q32P28-3", "Q86X76-3", "Q9UKX5-2")) #Asociados con EMT
plot_single(dep, proteins = c("EEF1A2","MICU2", "Q4KMQ2-2")) #relacionados con la homeostasis del calcio

plot_single(dep, proteins = "TGFBI", type = "centered")

# Plot a frequency plot of significant proteins for the different conditions
plot_cond(dep)

# Generate a results table
data_results <- get_results(dep)
view(data_results)
write.csv(data_results, file="data_total_HelaKO.csv",row.names = FALSE)

# Number of significant proteins
data_significant <- data_results %>% filter(significant)
view(data_significant)
write.csv(data_significant, file="data_significant_HelaKO.csv",row.names = FALSE)

#######GUARDAR SUMMARIZED EXPERIMENT PARA FUTUROS ANÁLISIS############
# Generate a wide data.frame
df_wide <- get_df_wide(dep)
write.csv(df_wide, file="df_wide_total_HelaKO.csv",row.names = FALSE)

# Generate a long data.frame
df_long <- get_df_long(dep)
write.csv(df_long, file="df_long_total_HelaKO.csv",row.names = FALSE)

# Save analyzed data
save(data_se, data_norm, data_imp, data_diff, dep, file = "data_HelaALMS1KO.RData")
