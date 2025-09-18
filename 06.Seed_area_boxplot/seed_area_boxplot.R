####################################################################################################
# Seed area distribution analysis with boxplot construction
####################################################################################################

# Author: Yedomon Ange Bovys Zoclanclounon, PhD
# Affiliation: Rothamsted Research, Plant Sciences and the Bioeconomy
# Contact: X: @angeomics
####################################################################################################

# Clean environment and set random seed
rm(list = ls())
gc()
set.seed(1000)

# Set working directory
work.dir <- "D:/laura/Figure2_plot"
setwd(work.dir)

# Load required packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, ggplot2, ggstatsplot)

# Load data
data_sample <- read.csv("data_withn_300.csv")

# Boxplot with stats
ggbetweenstats(
  data = data_sample,  # corrected from 'filter_data_sample'
  x = genotype,
  y = Area,
  title = "Comparison of seed area between genotypes (equal number of observations)",
  xlab = "Genotype",
  ylab = "Seed area (cm²)"  # used superscript ² for clarity
)

######################################################################
# Mean and standard deviation calculation
######################################################################

summary_stats <- data_sample %>%
  group_by(genotype) %>%
  summarise(
    Mean_Area = mean(Area, na.rm = TRUE),
    SD_Area   = sd(Area, na.rm = TRUE)
  )

# Print the results
print(summary_stats)

######################################################################





