# Alpha diversity in R - qiime2 outpu
# We will use the following files created using the qiime2 moving pictures tutorial.

# core-metrics-results/evenness_vector.qza (alpha diversity)
# core-metrics-results/faith_pd_vector.qza (alpha diversity)
# core-metrics-results/observed_features_vector.qza (alpha diversity)
# core-metrics-results/shannon_vector.qza (alpha diversity)
# sample-metadata.tsv (metadata)

setwd("/Users/lauren/Downloads")
list.files()

library(tidyverse)
library(qiime2R)
library(ggpubr)

##---DIVIDED INTO 2 GROUPS
#---HEALTHY

meta<-read_q2metadata("metadatagrouped.txt")
str(meta)
colnames(meta)[3] <- "body.site"
colnames(meta)[8] <- "reported.antibiotic.usage"
colnames(meta)[9] <- "days.since.experiment.start"
str(meta)

#the here below we are reading this file and we save it as evenness
evenness = read_qza("evenness_vector.qza")
evenness<-evenness$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

observed_features = read_qza("observed_features_vector.qza")
observed_features<-observed_features$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

shannon = read_qza("shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

faith_pd = read_qza("faith_pd_vector.qza")
faith_pd<-faith_pd$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged\

#fixing faith
faith_pd <- faith_pd %>%
  select(V1, V2) %>%
  rename(
    SampleID = V1,
    faith_pd = V2
  )

head(faith_pd)
str(faith_pd)

##ALPHA DIVERSITY TABLES
#Step 1: Combine all alpha diversity metrics with "merge"

alpha_diversity = merge(x=faith_pd, y=evenness, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, observed_features, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, shannon, by.x = "SampleID", by.y = "SampleID")
meta = merge(meta, alpha_diversity, by.x = "SampleID", by.y = "SampleID")
row.names(meta) <- meta$SampleID
#meta = meta[,-1]
str(meta)

#Plots
#hist command makes histogram
hist(meta$shannon_entropy, main="Shannon diversity", xlab="", breaks=10)
hist(meta$faith_pd, main="Faith phylogenetic diversity", xlab="", breaks=10)
hist(meta$pielou_e, main="Evenness", xlab="", breaks=10)
hist(as.numeric(meta$observed_features), main="Observed Features", xlab="", breaks=10)

#Plots the qq-plot for residuals
ggqqplot(meta$shannon_entropy, title = "Shannon")
ggqqplot(meta$faith_pd, title = "Faith PD")
ggqqplot(meta$pielou_e, title = "Evenness")
ggqqplot(meta$observed_features, title = "Observed Features")

#after looking at them some turn out to be normal others not
install.packages("ggpubr")
library("ggpubr")

# To test for normalcy statistically, we can run the Shapiro-Wilk 
# test of normality.

shapiro.test(meta$shannon)
shapiro.test(meta$faith_pd)
shapiro.test(meta$pielou_e)
shapiro.test(meta$observed_features)

aov.evenness.body_site = aov(pielou_evenness ~ treatment, data=meta)





#---HIV