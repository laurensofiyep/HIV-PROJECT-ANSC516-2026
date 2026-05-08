#Install the packages, IF YOU NEED TO :)
install.packages("tidyverse")
install.packages("vegan")
install.packages("devtools")

devtools::install_github("jbisanz/qiime2R")

#Load the packages. Everyone needs to do this.
library(tidyverse)
library(vegan)
library(qiime2R)
library(devtools)

###Set your working directory path/to/ANSC516/ANSC-repo/data/moving-pictures
setwd("/Users/lauren/Downloads/core-metrics-results-hiv") #you can paste your directory to move there

list.files()

if(!dir.exists("output"))
  dir.create("output")

metadata<-read_q2metadata("metadatagrouped.txt")
str(metadata)

bc <- read_qza("bray_curtis_pcoa_results.qza")
jaccard <- read_qza("jaccard_pcoa_results.qza")
uwu <- read_qza("unweighted_unifrac_pcoa_results.qza")
wu <- read_qza("weighted_unifrac_pcoa_results.qza")

bc_meta <- bc$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

jaccard_meta <- jaccard$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

uwu_meta <- uwu$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

wu_meta <- wu$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

# ---------------------------
# COLOR + LABELS + SHAPES SETUP
# ---------------------------

treatment_colors <- c(
  "blank" = "coral",
  "FOS" = "darkolivegreen2",
  "FM" = "turquoise2",
  "blank36" = "purple"
)

treatment_labels <- c(
  "blank" = "Blank",
  "FOS" = "FOS",
  "FM" = "FM",
  "blank36" = "Blank 36h"
)

treatment_shapes <- c(
  "blank" = 16,   # circle
  "FOS" = 17,     # triangle
  "FM" = 15,      # square
  "blank36" = 18  # diamond
)

# reusable scales
color_scale <- scale_color_manual(
  values = treatment_colors,
  labels = treatment_labels,
  name = "Treatment"
)

shape_scale <- scale_shape_manual(
  values = treatment_shapes,
  labels = treatment_labels,
  name = "Treatment"
)

bcplot <- ggplot(bc_meta, aes(x = PC1, y = PC2, color = treatment, shape = treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  theme_q2r() +
  color_scale +
  shape_scale +
  labs(
    title = "Bray-Curtis (HIV)",
    x = paste0("PC1 (", round(100 * bc$data$ProportionExplained[1], 2), "%)"),
    y = paste0("PC2 (", round(100 * bc$data$ProportionExplained[2], 2), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
bcplot
ggsave("output/brayc_hiv.png", bcplot, height = 3, width = 4)


jaccplot <- ggplot(jaccard_meta, aes(x = PC1, y = PC2, color = treatment, shape = treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  theme_q2r() +
  color_scale +
  shape_scale +
  labs(
    title = "Jaccard (HIV)",
    x = paste0("PC1 (", round(100 * jaccard$data$ProportionExplained[1], 2), "%)"),
    y = paste0("PC2 (", round(100 * jaccard$data$ProportionExplained[2], 2), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
jaccplot
ggsave("output/jacc_hiv.png", jaccplot, height = 3, width = 4)

uwuplot <- ggplot(uwu_meta, aes(x = PC1, y = PC2, color = treatment, shape = treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  theme_q2r() +
  color_scale +
  shape_scale +
  labs(
    title = "Unweighted Unifrac (HIV)",
    x = paste0("PC1 (", round(100 * uwu$data$ProportionExplained[1], 2), "%)"),
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
uwuplot
ggsave("output/unweighunifrac_hiv.png", uwuplot, height = 3, width = 4)

wuplot <- ggplot(wu_meta, aes(x = PC1, y = PC2, color = treatment, shape = treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  theme_q2r() +
  color_scale +
  shape_scale +
  labs(
    title = "Weighted Unifrac (HIV)",
    x = paste0("PC1 (", round(100 * wu$data$ProportionExplained[1], 2), "%)"),
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

wuplot <- ggplot(wu_meta, aes(x = PC1, y = PC2, color = treatment, shape = treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "euclid", level = 0.95) +
  theme_q2r() +
  color_scale +
  shape_scale +
  labs(
    title = "Weighted Unifrac (HIV)",
    x = paste0("PC1 (", round(100 * wu$data$ProportionExplained[1], 2), "%)"),
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

wuplot
ggsave("output/weighunifrac_hiv.png", wuplot, height = 3, width = 4)

