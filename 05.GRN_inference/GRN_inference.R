####################################################################################################
# Regulatory Network Inference in the Top 1000 Differentially Expressed Genes in Endorsperm Tissues
####################################################################################################
# Author      : Yedomon Ange Bovys Zoclanclounon, PhD
# Affiliation : Rothamsted Research, Plant Sciences and the Bioeconomy
# Contact     : X (Twitter) @angeomics
####################################################################################################

# --- Environment Setup ---------------------------------------------------------------------------

# Clean workspace and set seed for reproducibility
rm(list = ls())
gc()
set.seed(1000)

# Define working directory
work.dir <- "D:/Ange_DELL/My_folder/Smita_Kurup/Sent/Code_availability/GRN"
setwd(work.dir)

# Load required packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(DIANE, visNetwork, tidyverse)

# --- Data Import ---------------------------------------------------------------------------------

# Load expression counts and transcription factor list
data_set_top1000_endo <- read.csv("expression_count_data.csv", header = TRUE, sep = ",", row.names = 1)
brassica_regulators   <- read.csv("Transcription_factors_data.csv", header = TRUE)

# Extract gene and regulator IDs
genes <- rownames(data_set_top1000_endo)
brassica_regulators_list <- as.character(brassica_regulators$TF)
regressors <- intersect(genes, brassica_regulators_list)

# --- Regulator Grouping --------------------------------------------------------------------------

# Group regressors to identify putative intermediary regulators
grouping <- DIANE::group_regressors(data_set_top1000_endo, genes, regressors)

# Define grouped data objects
grouped_counts    <- grouping$counts
grouped_targets   <- grouping$grouped_genes
grouped_regressors <- grouping$grouped_regressors

# --- Network Inference ---------------------------------------------------------------------------

set.seed(123) # ensure reproducibility
conds <- c("ZS11heart", "ZS11torpedo", "ZS11green", "RJheart", "RJtorpedo", "RJgreen")

mat <- DIANE::network_inference(grouped_counts, conds, grouped_targets, grouped_regressors,
                                nCores = 1, verbose = TRUE)

# Threshold the network
network <- DIANE::network_thresholding(mat, n_edges = length(genes))

# --- Network Evaluation --------------------------------------------------------------------------

# Edge testing based on density
nGenes      <- length(grouped_targets)
nRegulators <- length(grouped_regressors)

res <- DIANE::test_edges(mat,
                         normalized_counts = grouped_counts,
                         density = 0.03,
                         nGenes = nGenes,
                         nRegulators = nRegulators,
                         nTrees = 1000,
                         verbose = TRUE)

# Save results
write.csv(res$links, "endosperm_only_top_1000_diane_res_links_calculations.csv")

# --- Edge Density Curve --------------------------------------------------------------------------

density_values <- seq(0.001, 0.1, length.out = 20)
res_density <- data.frame(density = density_values,
                          nEdges   = sapply(density_values, get_nEdges, nGenes, nRegulators))

ggplot(res_density, aes(x = density, y = nEdges)) +
  geom_line(size = 1) +
  ggtitle("Number of Network Edges as a Function of Density") +
  geom_vline(xintercept = 0.03, color = "darkgreen", linetype = 2, size = 1.2)

# --- Network Annotation --------------------------------------------------------------------------

genes_description <- read.csv("gene_description_data.csv", header = TRUE, sep = ",", row.names = 1)

data <- network_data(network, brassica_regulators_list, genes_description)

# Community detection (Louvain algorithm)
data$nodes$group <- data$nodes$community
louvain_membership <- data$nodes$community
names(louvain_membership) <- data$nodes$id
write.csv(louvain_membership, "endospermonly_louvain_membership.csv")

# --- Network Visualization -----------------------------------------------------------------------

network_plot <- DIANE::draw_network(data$nodes, data$edges) %>%
  visLegend() %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visEdges(arrows = "from")

# Save network visualization to HTML
network_plot %>% visSave(file = "network_top1000_endo_selection_id_legend_arrow.html",
                         background = "white")

# Export nodes and edges
write.csv(data$nodes, "top1000_endo_nodes.csv")
write.csv(data$edges, "top1000_endo_edges.csv")

# --- Network Visualization using Gephi--------------------------------------------------------------

#Download Gephi at https://gephi.org/users/download/ add nodes and edges data to visualise the network.

####################################################################################################
# End of Script
####################################################################################################
