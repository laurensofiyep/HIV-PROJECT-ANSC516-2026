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
setwd("/Users/lauren/Downloads/core-metrics-results-healthy") #you can paste your directory to move there

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
##FILTERING FOR ONLY 36HR REMOVING BLANK 0HR

bc_meta_36 <- bc_meta %>%
  filter(treatment %in% c("blank36", "FOS", "FM"))

jaccard_meta_36 <- jaccard_meta %>%
  filter(treatment %in% c("blank36", "FOS", "FM"))

uwu_meta_36 <- uwu_meta %>%
  filter(treatment %in% c("blank36", "FOS", "FM"))

wu_meta_36 <- wu_meta %>%
  filter(treatment %in% c("blank36", "FOS", "FM"))

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

#Bray-Curtis PLOTS
# ---------------------------
bcplot <- ggplot(bc_meta_36, aes(x = PC1, y = PC2, color = treatment, shape = treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  theme_q2r() +
  color_scale +
  shape_scale +
  labs(
    title = "Bray-Curtis (Healthy)",
    x = paste0("PC1 (", round(100 * bc$data$ProportionExplained[1], 2), "%)"),
    y = paste0("PC2 (", round(100 * bc$data$ProportionExplained[2], 2), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
bcplot
ggsave("output/brayc_healthy36hr.png", bcplot, height = 3, width = 4)

#Jaccard PLOTS
# ---------------------------
jaccplot <- ggplot(jaccard_meta_36, aes(x = PC1, y = PC2, color = treatment, shape = treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  theme_q2r() +
  color_scale +
  shape_scale +
  labs(
    title = "Jaccard (Healthy)",
    x = paste0("PC1 (", round(100 * jaccard$data$ProportionExplained[1], 2), "%)"),
    y = paste0("PC2 (", round(100 * jaccard$data$ProportionExplained[2], 2), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
jaccplot
ggsave("output/jacc_healthy36.png", jaccplot, height = 3, width = 4)

# Unweighted Unifrac PLOTS
# ---------------------------
uwuplot <- ggplot(uwu_meta_36, aes(x = PC1, y = PC2, color = treatment, shape = treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  theme_q2r() +
  color_scale +
  shape_scale +
  labs(
    title = "Unweighted Unifrac (Healthy)",
    x = paste0("PC1 (", round(100 * uwu$data$ProportionExplained[1], 2), "%)"),
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
uwuplot
ggsave("output/unweighunifrac_healthy36.png", uwuplot, height = 3, width = 4)

# Weighted Unifrac PLOTS
# ---------------------------
wuplot <- ggplot(wu_meta_36, aes(x = PC1, y = PC2, color = treatment, shape = treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  theme_q2r() +
  color_scale +
  shape_scale +
  labs(
    title = "Weighted Unifrac (Healthy)",
    x = paste0("PC1 (", round(100 * wu$data$ProportionExplained[1], 2), "%)"),
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
wuplot
ggsave("output/weighunifrac_healthy36.png", wuplot, height = 3, width = 4)

# ---------------------------
#PERMANOVA
# ---------------------------
# Import distance matrices 
bc_dist <- read_qza("bray_curtis_distance_matrix.qza")
wu_dist <- read_qza("weighted_unifrac_distance_matrix.qza")
#Match metadata correctlya and check
bc_dm <- as.matrix(bc_dist$data)
wu_dm <- as.matrix(wu_dist$data)
#Filter ONLY 36h samples 
samples_36 <- metadata_sub$SampleID[metadata_sub$treatment %in% c("blank36", "FOS", "FM")]

bc_dm_36 <- bc_dm[samples_36, samples_36]
wu_dm_36 <- wu_dm[samples_36, samples_36]

meta_36 <- metadata_sub[metadata_sub$SampleID %in% samples_36, ]

#Run PERMANOVA

# Bray-Curtis
adonis_bc <- adonis2(bc_dm_36 ~ treatment, data = meta_36, permutations = 999)

# Weighted UniFrac
adonis_wu <- adonis2(wu_dm_36 ~ treatment, data = meta_36, permutations = 999)

adonis_bc
adonis_wu

