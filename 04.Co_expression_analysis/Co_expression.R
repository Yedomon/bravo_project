####################################################################################################
# Multivariate Co-expression Analysis with Coseq
####################################################################################################
# Author      : Yedomon Ange Bovys Zoclanclounon, PhD
# Affiliation : Rothamsted Research, Plant Sciences and the Bioeconomy
# Contact     : X (Twitter) @angeomics
####################################################################################################

# --- Environment Setup ---------------------------------------------------------------------------

# Clean workspace and set seed for reproducibility
rm(list = ls())   # Remove all objects from environment
gc()              # Perform garbage collection
set.seed(1000)    # Global seed for reproducibility

# Define working directory
work.dir <- "D:/Ange_DELL/My_folder/Smita_Kurup/Sent/Code_availability/Co_seq"
setwd(work.dir)

# Load required packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(DIANE, tidyverse)

# --- Data Import ---------------------------------------------------------------------------------

# Load normalized gene expression dataset for co-expression analysis
data_normalized_1000 <- read.csv("data_endo_reformatted_for_coseq_analysis.csv",
                                 header = TRUE, sep = ",", row.names = 1)

# --- Co-expression Clustering --------------------------------------------------------------------

# Run co-expression clustering (Coseq via DIANE)
clustering <- DIANE::run_coseq(
  conds   = c("A1", "A2", "B1", "B2", "C1", "C2"),  # Experimental conditions
  data    = data_normalized_1000,                   # Normalized expression data
  genes   = rownames(data_normalized_1000),         # Gene IDs
  K       = 6:9,                                    # Range of clusters to explore
  transfo = "arcsin",                               # Transformation method
  model   = "Normal",                               # Clustering model
  seed    = 123                                     # Internal seed for reproducibility
)

# Save cluster membership to file
write.csv(as.data.frame(clustering$membership),
          "endospermonly_top_1000_coseq_clustering_membership.csv")

# --- Visualization of Clustering Results ---------------------------------------------------------

# 1. Model probability bar plots
pdf("coseq_model_probability_barplots.pdf", width = 10, height = 8)
DIANE::draw_coseq_run(clustering$model, plot = "barplots")
dev.off()

# 2. Elbow plot (Integrated Completed Likelihood) for optimal K selection
pdf("coseq_elbow_plot.pdf", width = 10, height = 8)
DIANE::draw_coseq_run(clustering$model, plot = "ICL")
dev.off()

# 3. Cluster profile visualization
profiles_plot <- DIANE::draw_profiles(
  data  = data_normalized_1000,
  clustering$membership,
  conds = c("A1", "A2", "B1", "B2", "C1", "C2")
) + theme_bw()

ggsave("gene_cluster_profiles.pdf", plot = profiles_plot,
       width = 10, height = 8, device = "pdf")

####################################################################################################
# End of Script
####################################################################################################
