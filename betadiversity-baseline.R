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
setwd("/Users/lauren/Downloads/core-metrics-results-baseline") #you can paste your directory to move there

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

group_colors <- c(
  "hiv" = "coral",
  "healthy" = "darkolivegreen2"
)
group_labels <- c(
  "hiv" = "HIV patients",
  "healthy" = "Healthy"
)
group_shapes <- c(
  "hiv" = 16,   # circle
  "healthy" = 17     # triangle
)
# reusable scales
color_scale <- scale_color_manual(
  values = group_colors,
  labels = group_labels,
  name = "group"
)
shape_scale <- scale_shape_manual(
  values = group_shapes,
  labels = group_labels,
  name = "group"
)
# ---------------------------
#PLOTS
# ---------------------------
#Bray curtis PLOTS
# ---------------------------

bcplot <- ggplot(bc_meta, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  theme_q2r() +
  color_scale +
  shape_scale +
  labs(
    title = "Bray-Curtis Baseline Healthy vs. HIV",
    x = paste0("PC1 (", round(100 * bc$data$ProportionExplained[1], 2), "%)"),
    y = paste0("PC2 (", round(100 * bc$data$ProportionExplained[2], 2), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
bcplot
ggsave("output/brayc_baseline.png", bcplot, height = 3, width = 4)

#Jaccard PLOTS
# ---------------------------
jaccplot <- ggplot(jaccard_meta, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  theme_q2r() +
  color_scale +
  shape_scale +
  labs(
    title = "Jaccard Baseline Healthy vs. HIV",
    x = paste0("PC1 (", round(100 * jaccard$data$ProportionExplained[1], 2), "%)"),
    y = paste0("PC2 (", round(100 * jaccard$data$ProportionExplained[2], 2), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
jaccplot
ggsave("output/jacc_baseline.png", jaccplot, height = 3, width = 4)

#Unweighted unifrac PLOTS
# ---------------------------
uwuplot <- ggplot(uwu_meta, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  theme_q2r() +
  color_scale +
  shape_scale +
  labs(
    title = "Unweighted Unifrac Baseline Healthy vs. HIV",
    x = paste0("PC1 (", round(100 * uwu$data$ProportionExplained[1], 2), "%)"),
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
uwuplot
ggsave("output/unweighunifrac_baseline.png", uwuplot, height = 3, width = 4)

#Weighted unifrac PLOTS
# ---------------------------
wuplot <- ggplot(wu_meta, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  theme_q2r() +
  color_scale +
  shape_scale +
  labs(
    title = "Weighted Unifrac Baseline Healthy vs. HIV",
    x = paste0("PC1 (", round(100 * wu$data$ProportionExplained[1], 2), "%)"),
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
wuplot
ggsave("output/weighunifrac_baseline.png", wuplot, height = 3, width = 4)

# ---------------------------
#PERMANOVA
# ---------------------------
bc_dm <- read_qza("bray_curtis_distance_matrix.qza")$data
bc_dm <- as.matrix(bc_dm)

bc_meta <- bc$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = "SampleID")

adonis2(bc_dm ~ group, data = metadata)

bc_dist_mat<-read_qza("bray_curtis_distance_matrix.qza")
bc_dm <- as.matrix(bc_dist_mat$data) 
rownames(bc_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(bc_dm),metadata$SampleID),]
rownames(bc_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(bc_dm ~ group, data = metadata_sub)

write.table(PERMANOVA_out,"output/Body.site_Adonis_overall.csv",sep=",", row.names = TRUE) 


adonis2(bc_dm ~ group, data = metadata)
