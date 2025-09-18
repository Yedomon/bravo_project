############################################################
##### Subgenome Dominance Analysis #########################
############################################################

# Author: Yedomon Ange Bovys Zoclanclounon, PhD
# Affiliation: Rothamsted Research, Plant Sciences and the Bioeconomy
# Contact: X: angeomics
############################################################

# --- Prepare Environment ---
rm(list = ls())   # Clear workspace
gc()              # Garbage collection
set.seed(1000)    # Reproducibility

# --- Working Directory ---
work.dir <- "D:/laura"
setwd(work.dir)

# --- Load Required Packages ---
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  dplyr,
  tidyverse,
  ggstatsplot,
  patchwork
)

############################################################
## Tissue: Endosperm
############################################################

# Load data
data_zs11 <- read.csv("sub_genome_dominance_endo_zs11.csv", header = TRUE)

# Generate plots for each stage
set.seed(123)

endo_heart <- data_zs11 %>% filter(Stage == "Heart")
plot_heart <- ggbetweenstats(
  data  = endo_heart,
  x     = Subgenome,
  y     = log2Normalized_count,
  xlab  = "Subgenome",
  ylab  = "log2(Normalized count + 1)",
  p.adjust.method = "bonferroni",
  title = "Heart"
)

endo_torpedo <- data_zs11 %>% filter(Stage == "Torpedo")
plot_torpedo <- ggbetweenstats(
  data  = endo_torpedo,
  x     = Subgenome,
  y     = log2Normalized_count,
  xlab  = "Subgenome",
  ylab  = "log2(Normalized count + 1)",
  p.adjust.method = "bonferroni",
  title = "Torpedo"
)

endo_green <- data_zs11 %>% filter(Stage == "Green")
plot_green <- ggbetweenstats(
  data  = endo_green,
  x     = Subgenome,
  y     = log2Normalized_count,
  xlab  = "Subgenome",
  ylab  = "log2(Normalized count + 1)",
  p.adjust.method = "bonferroni",
  title = "Green"
)

# Combined analysis (Top 1000, all stages)
plot_top1000 <- ggbetweenstats(
  data  = data_zs11,
  x     = Subgenome,
  y     = log2Normalized_count,
  xlab  = "Subgenome",
  ylab  = "log2(Normalized count + 1)",
  p.adjust.method = "bonferroni",
  title = "Top 1000 (All Stages)"
)

# Combine endosperm plots
plot_heart + plot_torpedo + plot_green + plot_top1000 + plot_layout(nrow = 1)

############################################################
## Tissue: Seed Coat
############################################################

# Load and reshape data
seedcoat_data <- read.csv("data_subgenome_top1000_seed_coat_heart_ZS11.csv", header = TRUE, sep = ",") %>%
  pivot_longer(cols = -c(Gene, Subgenome), names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = factor(Variable, levels = paste0(c("Heart", "Torpedo", "Green", "Mature"), rep(1:4, each = 4))))

# Generate stage-specific plots
set.seed(123)

seedcoat_heart <- seedcoat_data %>%
  filter(grepl("Heart", Variable)) %>%
  mutate(Value = log2(Value + 1))
plot_sc_heart <- ggbetweenstats(
  data  = seedcoat_heart,
  x     = Subgenome,
  y     = Value,
  xlab  = "Subgenome",
  ylab  = "log2(Normalized count + 1)",
  p.adjust.method = "bonferroni",
  title = "Top1000: Seed Coat (ZS11, Heart)"
)

seedcoat_torpedo <- seedcoat_data %>%
  filter(grepl("Torpedo", Variable)) %>%
  mutate(Value = log2(Value + 1))
plot_sc_torpedo <- ggbetweenstats(
  data  = seedcoat_torpedo,
  x     = Subgenome,
  y     = Value,
  xlab  = "Subgenome",
  ylab  = "log2(Normalized count + 1)",
  p.adjust.method = "bonferroni",
  title = "Top1000: Seed Coat (ZS11, Torpedo)"
)

seedcoat_green <- seedcoat_data %>%
  filter(grepl("Green", Variable)) %>%
  mutate(Value = log2(Value + 1))
plot_sc_green <- ggbetweenstats(
  data  = seedcoat_green,
  x     = Subgenome,
  y     = Value,
  xlab  = "Subgenome",
  ylab  = "log2(Normalized count + 1)",
  p.adjust.method = "bonferroni",
  title = "Top1000: Seed Coat (ZS11, Green)"
)

seedcoat_mature <- seedcoat_data %>%
  filter(grepl("Mature", Variable)) %>%
  mutate(Value = log2(Value + 1))
plot_sc_mature <- ggbetweenstats(
  data  = seedcoat_mature,
  x     = Subgenome,
  y     = Value,
  xlab  = "Subgenome",
  ylab  = "log2(Normalized count + 1)",
  p.adjust.method = "bonferroni",
  title = "Top1000: Seed Coat (ZS11, Mature)"
)

# Combine seed coat plots
plot_sc_heart + plot_sc_torpedo + plot_sc_green + plot_sc_mature + plot_layout(nrow = 1)

############################################################
## Tissue: Embryo
############################################################

# Load and reshape data
embryo_data <- read.csv("data_subgenome_top1000_embryo_heart_ZS11.csv", header = TRUE, sep = ",") %>%
  pivot_longer(cols = -c(Gene, Subgenome), names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = factor(Variable, levels = paste0(c("Heart", "Torpedo", "Green", "Mature"), rep(1:4, each = 4))))

# Generate stage-specific plots
set.seed(123)

embryo_heart <- embryo_data %>%
  filter(grepl("Heart", Variable)) %>%
  mutate(Value = log2(Value + 1))
plot_em_heart <- ggbetweenstats(
  data  = embryo_heart,
  x     = Subgenome,
  y     = Value,
  xlab  = "Subgenome",
  ylab  = "log2(Normalized count + 1)",
  p.adjust.method = "bonferroni",
  title = "Top1000: Embryo (ZS11, Heart)"
)

embryo_torpedo <- embryo_data %>%
  filter(grepl("Torpedo", Variable)) %>%
  mutate(Value = log2(Value + 1))
plot_em_torpedo <- ggbetweenstats(
  data  = embryo_torpedo,
  x     = Subgenome,
  y     = Value,
  xlab  = "Subgenome",
  ylab  = "log2(Normalized count + 1)",
  p.adjust.method = "bonferroni",
  title = "Top1000: Embryo (ZS11, Torpedo)"
)

embryo_green <- embryo_data %>%
  filter(grepl("Green", Variable)) %>%
  mutate(Value = log2(Value + 1))
plot_em_green <- ggbetweenstats(
  data  = embryo_green,
  x     = Subgenome,
  y     = Value,
  xlab  = "Subgenome",
  ylab  = "log2(Normalized count + 1)",
  p.adjust.method = "bonferroni",
  title = "Top1000: Embryo (ZS11, Green)"
)

embryo_mature <- embryo_data %>%
  filter(grepl("Mature", Variable)) %>%
  mutate(Value = log2(Value + 1))
plot_em_mature <- ggbetweenstats(
  data  = embryo_mature,
  x     = Subgenome,
  y     = Value,
  xlab  = "Subgenome",
  ylab  = "log2(Normalized count + 1)",
  p.adjust.method = "bonferroni",
  title = "Top1000: Embryo (ZS11, Mature)"
)

# Combine embryo plots
plot_em_heart + plot_em_torpedo + plot_em_green + plot_em_mature + plot_layout(nrow = 1)

############################################################
# End of Script
############################################################
