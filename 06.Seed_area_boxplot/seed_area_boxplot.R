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
pacman::p_load(tidyverse, ggplot2, ggstatsplot,car)

# Load data
data_sample <- read.csv("data_withn_300.csv")

# Shapiro–Wilk test per genotype

data_sample %>%
  group_by(genotype) %>%
  summarise(
    n = n(),
    p_value = shapiro.test(Area)$p.value
  )


##########################################################
# A tibble: 4 × 3
#  genotype     n  p_value
#  <chr>    <int>    <dbl>
#1 Express    300 2.32e-11
#2 RJ         300 2.45e- 2
#3 Stellar    300 5.85e- 3
#4 ZS11       300 1.37e- 1
##########################################################

# Levene’s Test (data are not following a normal distribution)

leveneTest(Area ~ genotype, data = data_sample)

#########################################################
#Levene's Test for Homogeneity of Variance (center = median)
#        Df F value    Pr(>F)    
#group    3  36.735 < 2.2e-16 ***
#      1196                      

#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#########################################################

# Boxplot with stats (Welch's one-way ANOVA)
ggbetweenstats(
  data = data_sample,  
  x = genotype,
  y = Area,
  title = "Comparison of seed area between genotypes (equal number of observations)",
  xlab = "Genotype",
  ylab = "Seed area (cm²)"  
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






