################################################################################
### Differential Expression Analysis in Endosperm Across Developmental Stages ###
################################################################################

# Author: Yedomon Ange Bovys Zoclanclounon, PhD
# Affiliation: Rothamsted Research, Plant Sciences and the Bioeconomy
# Contact: X: @angeomics
################################################################################

# Clean environment and set random seed
rm(list = ls())
gc()
set.seed(1000)

# Set working directory
work_dir <- "D:/Ange_DELL/My_folder/Smita_Kurup/Sent/Code_availability/DEG"
setwd(work_dir)

# Load required packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(DESeq2, tidyverse, pheatmap, ggplot2)

# ------------------------------------------------------------------------------ #
# Load input data
# ------------------------------------------------------------------------------ #
counts_data <- read.csv("endosperm_data_only.csv", header = TRUE, row.names = 1)
colData     <- read.csv("design_endosperm_only.csv", header = TRUE, row.names = 1)

# Ensure design file matches the counts data
stopifnot(all(colnames(counts_data) == rownames(colData)))

# Construct DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData   = colData,
  design    = ~ Species + Stage
)
saveRDS(dds, "ExpressionSet_endosperm_only.RDS")

# ------------------------------------------------------------------------------ #
# PCA analysis
# ------------------------------------------------------------------------------ #
rld <- rlog(dds)  
data <- plotPCA(rld, intgroup = c("Species", "Stage"), returnData = TRUE)

# Extract variance explained
percentVar <- round(100 * attr(data, "percentVar"))

# PCA plot
p <- ggplot(data, aes(PC1, PC2, color = group, shape = Species)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_polygon(aes(fill = group), alpha = 0.1, position = position_jitter(width = 0.3, height = 0.3)) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.key   = element_rect(fill = "white", color = "white"),
    plot.title   = element_text(hjust = 0.5)
  )

# Custom colors
cbp1 <- c("#008585", "#62380F", "#AB4817", "#F07575", "#b6b975", "#315c45")
p + scale_color_manual(values = cbp1) + scale_fill_manual(values = cbp1)

# Save plot
ggsave("PCA_plot_endosperm.pdf", plot = p, width = 10, height = 8, dpi = 300)

# ------------------------------------------------------------------------------ #
# DEG analysis (Species × Stage interaction)
# ------------------------------------------------------------------------------ #
dds_endo <- DESeqDataSetFromMatrix(countData = counts_data, colData = colData, design = ~ Species_Stage)
dds_endo <- dds_endo[rowSums(counts(dds_endo)) >= 10, ]  # pre-filtering
dds_endo <- DESeq(dds_endo)

saveRDS(dds_endo, "dds_DESeq_object_endo_only.RDS")

# Results
res_endo <- results(dds_endo)
res_endo_df <- as.data.frame(res_endo)

# Top 1000 DEGs
top_1000_genes <- res_endo_df[order(res_endo_df$padj, na.last = NA), ][1:1000, ]
normalized_counts <- counts(dds_endo, normalized = TRUE)
top_1000_counts <- normalized_counts[rownames(top_1000_genes), ]

top_1000_combined <- cbind(GeneID = rownames(top_1000_genes), top_1000_genes, top_1000_counts)
write.csv(top_1000_combined, "Top1000_DEGs_Endosperm.csv", row.names = FALSE)

# ------------------------------------------------------------------------------ #
# Stage-wise DEG analysis
# ------------------------------------------------------------------------------ #
stages <- c("Heart", "Torpedo", "Green")
results_list <- list()

for (stage in stages) {
  contrast <- c("Species_Stage", paste0("ZS11", stage), paste0("RJ", stage))
  res <- results(dds_endo, contrast = contrast)
  
  saveRDS(res, paste0("res_DESeq_endosperm_ZS11_", stage, "_vs_RJ_", stage, ".RDS"))
  write.csv(as.data.frame(res), paste0("DEG_results_", stage, ".csv"), row.names = TRUE)
  
  # Count DEGs
  upregulated   <- sum(res$padj < 0.05 & res$log2FoldChange > 2, na.rm = TRUE)
  downregulated <- sum(res$padj < 0.05 & res$log2FoldChange < -2, na.rm = TRUE)
  
  results_list[[stage]] <- list(upregulated = upregulated, downregulated = downregulated)
}

# Save DEG summary
sink("DEG_summary.txt")
for (stage in names(results_list)) {
  cat("\nStage:", stage, "\n")
  cat("Upregulated genes:", results_list[[stage]]$upregulated, "\n")
  cat("Downregulated genes:", results_list[[stage]]$downregulated, "\n")
}
sink()
cat("DEG summary has been saved to DEG_summary.txt\n")

# ------------------------------------------------------------------------------ #
# Heatmap of selected genes
# ------------------------------------------------------------------------------ #
expredata <- read.csv("gene_of_interest_endo_only_v0.1.csv", header = TRUE, row.names = 1)
data_matrix <- as.matrix(expredata[1:18])

rowannot_data      <- read.csv("rowannotation_file.csv", header = TRUE, row.names = 1)
annotation_col_data <- read.csv("annotation_col_endo.csv", row.names = 1)

ann_colors <- list(
  species = c(ZS11 = "#b6b975", RJ = "firebrick"),
  stage   = c(Heart = "#1B9E77", Torpedo = "#D95F02", Green = "#8a9a65", Mature = "#F1C40F"),
  Type    = c(
    Seed_cell_expansion                    = "#c7522a",
    Seed_cell_elongation                   = "#e5c185",
    Cell_wall_development_and_lignification = "#fbf2c4",
    Maternal_control_of_seed_size          = "#74a892",
    Transcriptional_regulator              = "#008585",
    Seed_development_and_size              = "#ff6361",
    Phytohormone_signaling_and_homeostasis = "#194a7a",
    Seed_size_control_by_HAIKU_pathway_regulation = "#183e4b"
  )
)

set.seed(12345)
pheatmap(
  data_matrix,
  cellwidth  = 20,
  cellheight = 20,
  scale = "row",
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = annotation_col_data,
  annotation_row = rowannot_data,
  annotation_colors = ann_colors,
  cutree_cols = 2,
  border_color = "white",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  fontsize = 20
)

# ------------------------------------------------------------------------------ #
# Mirror barplot of DEG counts
# ------------------------------------------------------------------------------ #
endo_data <- read.csv("mirror_barplot_endo_only.csv", header = TRUE)

bar_endo <- ggplot(endo_data, aes(x = Stage, y = Count, fill = Type)) + 
  geom_bar(stat = "identity", position = "identity") +
  geom_label(aes(label = Count), hjust = 0.5, vjust = 1.5, colour = "white", size = 6)

mycolor <- c("Up-regulated" = "#2980B9", "Down-regulated" = "#C0392B")

bar_endo + 
  scale_fill_manual(values = mycolor) +
  coord_flip()

# ------------------------------------------------------------------------------ #
# Venn diagrams
# ------------------------------------------------------------------------------ #
# For venn diagram use: https://bioinformatics.psb.ugent.be/webtools/Venn/
# Use the DEG lists generated per stage above
# ------------------------------------------------------------------------------ #

