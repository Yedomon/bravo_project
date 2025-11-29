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
  cowplot
)

tpm_data = read.csv("zho_tpm_data.csv", header = TRUE)

df_tidy <- tpm_data %>%
  pivot_longer(
    cols = -GeneID,
    names_to = c("Stage", "Organ", "Replicate"),
    names_pattern = "(.*)_(.*)_rep(\\d+)",
    values_to = "TPM"
  ) %>%
  mutate(
    Stage     = str_to_title(Stage),
    Organ     = str_to_title(Organ),
    Replicate = paste0("Rep", Replicate)
  ) %>%
  mutate(
    Subgenome = case_when(
      str_sub(GeneID, 4, 4) %in% c("A", "C", "U") ~ str_sub(GeneID, 4, 4),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(Subgenome %in% c("A", "C")) %>%      #  KEEP ONLY A and C
  filter(TPM >= 1) %>%
  mutate(
    log2TPM = log2(TPM + 1)
  )

df_tidy

 

############################################################
## Tissue: Seed Coat
############################################################


# A function to generate the plots

plot_seedcoat_stage <- function(stage_name) {
  df_tidy %>%
    filter(Organ == "Seedcoat", Stage == stage_name) %>%
    ggbetweenstats(
      x     = Subgenome,
      y     = log2TPM,
      xlab  = "Subgenome",
      ylab  = "log2(TPM + 1)",
      title = paste("Seed Coat (ZS11", stage_name, ")"),
      ggtheme = ggplot2::theme_classic(),
      ggplot.component = list(
        scale_fill_manual(
          values = c(
            "A" = "#5e0c3f66",   # color for A
            "C" = "#5c389966"    # color for C
          )
        ),
        scale_color_manual(
          values = c(
            "A" = "#5e0c3f66",
            "C" = "#5c389966"
          )
        )
      )
    )
}


# Generate the plots

plot_sc_heart   <- plot_seedcoat_stage("Heart")
plot_sc_torpedo <- plot_seedcoat_stage("Torpedo")
plot_sc_green   <- plot_seedcoat_stage("Green")
plot_sc_mature  <- plot_seedcoat_stage("Mature")


# Combine seed coat plots

library(cowplot)


combined_seedcoat_plot <- cowplot::plot_grid(
  plot_sc_heart,
  plot_sc_torpedo,
  plot_sc_green,
  plot_sc_mature,
  nrow = 1,
  labels = NULL,
  align = "h",
  rel_widths = c(1,1,1,1)
)

ggsave(
  filename = "seedcoat_subgenome_boxplots.png",
  plot = combined_seedcoat_plot,
  width = 16, height = 5, dpi = 300
)

ggsave(
  filename = "seedcoat_subgenome_boxplots.pdf",
  plot = combined_seedcoat_plot,
  width = 16, height = 5
)

ggsave(
  filename = "seedcoat_subgenome_boxplots.svg",
  plot = combined_seedcoat_plot,
  width = 16, height = 5
)




############################################################
## Tissue: Embryo
############################################################

plot_embryo_stage <- function(stage_name) {
  df_tidy %>%
    filter(Organ == "Embryo", Stage == stage_name) %>%
    ggbetweenstats(
      x     = Subgenome,
      y     = log2TPM,
      xlab  = "Subgenome",
      ylab  = "log2(TPM + 1)",
      title = paste("Embryo (ZS11", stage_name, ")"),
      ggtheme = ggplot2::theme_classic(),
      ggplot.component = list(
        scale_fill_manual(
          values = c("A" = "#5e0c3f66", "C" = "#5c389966")
        ),
        scale_color_manual(
          values = c("A" = "#5e0c3f66", "C" = "#5c389966")
        )
      )
    )
}

# Generate plots
plot_em_heart   <- plot_embryo_stage("Heart")
plot_em_torpedo <- plot_embryo_stage("Torpedo")
plot_em_green   <- plot_embryo_stage("Green")
plot_em_mature  <- plot_embryo_stage("Mature")

# Combine
combined_embryo_plot <- cowplot::plot_grid(
  plot_em_heart,
  plot_em_torpedo,
  plot_em_green,
  plot_em_mature,
  nrow = 1,
  labels = NULL,
  align = "h"
)

# Save
ggsave("embryo_subgenome_boxplots.png", combined_embryo_plot, width = 16, height = 5, dpi = 300)
ggsave("embryo_subgenome_boxplots.pdf", combined_embryo_plot, width = 16, height = 5)
ggsave("embryo_subgenome_boxplots.svg", combined_embryo_plot, width = 16, height = 5)



############################################################
## Tissue: Endosperm
############################################################

plot_endosperm_stage <- function(stage_name) {
  df_tidy %>%
    filter(Organ == "Endosperm", Stage == stage_name) %>%
    ggbetweenstats(
      x     = Subgenome,
      y     = log2TPM,
      xlab  = "Subgenome",
      ylab  = "log2(TPM + 1)",
      title = paste("Endosperm (ZS11", stage_name, ")"),
      ggtheme = ggplot2::theme_classic(),
      ggplot.component = list(
        scale_fill_manual(
          values = c("A" = "#5e0c3f66", "C" = "#5c389966")
        ),
        scale_color_manual(
          values = c("A" = "#5e0c3f66", "C" = "#5c389966")
        )
      )
    )
}

# Generate plots
plot_en_heart   <- plot_endosperm_stage("Heart")
plot_en_torpedo <- plot_endosperm_stage("Torpedo")
plot_en_green   <- plot_endosperm_stage("Green")
#plot_en_mature  <- plot_endosperm_stage("Mature") # Not available so I will replace with all combined TPM

## Combined All TPM

combined_all_plot <- ggbetweenstats(
  data  = df_tidy,
  x     = Subgenome,
  y     = log2TPM,
  xlab  = "Subgenome",
  ylab  = "log2(TPM + 1)",
  title = "All Organs and All Stages Combined",
  ggtheme = ggplot2::theme_classic(),
  ggplot.component = list(
    scale_fill_manual(
      values = c(
        "A" = "#5e0c3f66",   # color for A
        "C" = "#5c389966"    # color for C
      )
    ),
    scale_color_manual(
      values = c(
        "A" = "#5e0c3f66",
        "C" = "#5c389966"
      )
    )
  )
)




# Combine
combined_endosperm_plot <- cowplot::plot_grid(
  plot_en_heart,
  plot_en_torpedo,
  plot_en_green,
  combined_all_plot,
  nrow = 1,
  labels = NULL,
  align = "h"
)

# Save
ggsave("endosperm_subgenome_boxplots.png", combined_endosperm_plot, width = 16, height = 5, dpi = 300)
ggsave("endosperm_subgenome_boxplots.pdf", combined_endosperm_plot, width = 16, height = 5)
ggsave("endosperm_subgenome_boxplots.svg", combined_endosperm_plot, width = 16, height = 5)


#Final Combined 3-Panel Figure (A–C)

final_threepanel_plot <- cowplot::plot_grid(
  combined_endosperm_plot,
  combined_embryo_plot,
  combined_seedcoat_plot,
  nrow = 3,
  labels = c("(A)", "(B)", "(C)"),   # A = Endosperm, B = Embryo, C = Seed coat
  label_size = 22,
  align = "v",
  rel_heights = c(1, 1, 1)     # equal height
)

# SAVE in PNG
ggsave(
  "subgenome_endosperm_embryo_seedcoat_3panel.png",
  final_threepanel_plot,
  width = 16, height = 15, dpi = 300
)

# SAVE in PDF
ggsave(
  "subgenome_endosperm_embryo_seedcoat_3panel.pdf",
  final_threepanel_plot,
  width = 16, height = 15
)

# SAVE in SVG
ggsave(
  "subgenome_endosperm_embryo_seedcoat_3panel.svg",
  final_threepanel_plot,
  width = 16, height = 15
)
############################################################
# End of Script
############################################################

