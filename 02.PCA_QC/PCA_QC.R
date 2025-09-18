############################################################
##### Differential Expression & PCA using all samples   ####
############################################################

# Author: Yedomon Ange Bovys Zoclanclounon, PhD
# Affiliation: Rothamsted Research, Plant Sciences and the Bioeconomy
# Contact: X (Twitter): @angeomics
############################################################

# Clean environment and set random seed
rm(list = ls())   # Clear environment
gc()              # Garbage collection
set.seed(1000)    # Ensure reproducibility

# Set working directory
work_dir <- "D:/Ange_DELL/My_folder/Smita_Kurup/Sent/Code_availability/PCA_QC"
setwd(work_dir)

# Load required packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(DESeq2, tidyverse, ggplot2)

############################################################
# Load Input Data
############################################################

# Count matrix
counts_data <- read.csv("Zhongshuang11_vs_Ragged_Jack_est_count_data.csv", 
                        header = TRUE, row.names = 1)

# Sample design metadata
colData <- read.csv("design.csv", header = TRUE, row.names = 1)

# Verify sample metadata alignment with count data
stopifnot(all(colnames(counts_data) %in% rownames(colData)))  # Should return TRUE
stopifnot(all(colnames(counts_data) == rownames(colData)))    # Should return TRUE

############################################################
# Construct DESeq Object
############################################################

dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData   = colData,
  design    = ~ tissue + stage + species
)

dds

############################################################
# PCA Analysis
############################################################

# Regularized log transformation
rld <- rlogTransformation(dds)

# Extract PCA data
data <- plotPCA(rld, intgroup = c("tissue", "stage", "species"), returnData = TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

# PCA Plot
p <- ggplot(data, aes(PC1, PC2, color = group, shape = species)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_polygon(aes(fill = group), alpha = 0.1, position = position_jitter(width = 0.3, height = 0.3)) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 18, family = "sans"),
    legend.position = "bottom",
    legend.text = element_text(size = 18, family = "sans"),
    axis.text.x = element_text(color = "black", size = 16, family = "sans"),
    axis.text.y = element_text(color = "black", size = 16, family = "sans"),
    axis.title.x = element_text(color = "black", size = 18, family = "sans"),
    legend.title = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(size = 10, hjust = 0.5),
    legend.key = element_rect(fill = "white", color = "white")
  )

p

############################################################
# Custom Color Palette
############################################################

cbp1 <- c(
  "#413344", "#614c65", "#806485", "#936397", "#a662a8", "#664972",
  "#463c57", "#6e8da9", "#91bcdd", "#567d99", "#395e77", "#305662",
  "#264d4d", "#315c45", "#8a9a65", "#b6b975", "#b65d54", "#b60033",
  "#98062d", "#800022", "#008585", "#d9a7b0"
)

p <- p + 
  scale_color_manual(values = cbp1) +
  scale_fill_manual(values = cbp1)

############################################################
# Save PCA Plot
############################################################

ggsave("PCA_plot_QC.pdf", plot = p, width = 12, height = 8, dpi = 300)

############################################################
# End of Script
############################################################
