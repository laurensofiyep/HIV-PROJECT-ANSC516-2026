## This script goes through making a taxa bar plot and then running
# to install phyloseq:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

library(qiime2R)
library(phyloseq)
#library(zoo)
library(tidyverse)
library(DESeq2)

##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
##############################################
getwd
setwd('/Users/lauren/Downloads/core-metrics-results-baseline')
list.files()

if(!dir.exists("output-taxa"))
  dir.create("output-taxa")

##Qiime2r method of reading in the taxonomy, metadata and table files individually. 
##Run these lines to troubleshoot if you have trouble on line 55.
taxonomy<-read_qza("taxonomy.qza")
head(taxonomy$data)
taxonomy_table<-parse_taxonomy(taxonomy$data)

metadata<-read_q2metadata("metadatagrouped.txt")
str(metadata)

rare_table <- read_qza("rarefied_table.qza")
feature_table <- rare_table$data

##Qiime2R method of creating a phyloseq object
physeq <- qza_to_phyloseq(
  features="rarefied_table.qza",
  tree="rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadatagrouped.txt"
)

asv_table <- data.frame(otu_table(physeq), check.names = F)
metadata <- data.frame(sample_data(physeq), check.names = F)
taxonomy <- data.frame(tax_table(physeq), check.names = F)

#Clean up metadata, slightly
metadata$group.ord = factor(
  metadata$group,
  c("healthy", "hiv")
)

#Clean up taxonomy
head(taxonomy)
tax.clean <- taxonomy

#All this is OK except that in future use of the taxonomy table, 
#these ASVs designated as NA will be ignored because they are not classified. 
#Why are ASVs not classified? Its because there is not a close enough 
#match in the database. Just because there is not a good match in 
#the database does not mean they don’t exist, so I wanted to make 
#sure this data was not lost. So in my new code, from lines 200 – 224 
#I make it so that ASVs that are unclassified at any level are 
#classified as the lowest taxonomic level for which there is a 
#classification.


tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("uncl_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("uncl_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("uncl_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("uncl_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("uncl_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("uncl_",tax.clean$Genus[i], sep = "_")
  }
}



#################################################################
##Taxa barplot
#################################################################


#Assign our edited and formatted tables as variables to be feed into phyloseq
OTU.physeq = otu_table(as.matrix(asv_table), taxa_are_rows=TRUE)
tax.physeq = tax_table(as.matrix(tax.clean))    
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_bar_plot = phyloseq(OTU.physeq, tax.physeq, meta.physeq)


# Set colors for plotting
my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)

#If you want different taxonomic level, find and replace the taxonomic level listed here
my_level <- c("Phylum", "Family", "Genus")
my_column <- "group"  #this is the metadata column that we will use in the taxa barplot
my_column_ordered <- c("healthy", "hiv")

abund_filter <- 0.01  # Our abundance threshold
#ml ="Family"

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plot %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  physeq.taxa.mean <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.mean=mean(Abundance.average))
  
  physeq.taxa.mean <- as.data.frame(physeq.taxa.mean)
  colnames(physeq.taxa.mean)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  physeq_meta <- merge (physeq_meta, physeq.taxa.mean)
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  
  physeq_meta_filtered$my_column_ordered = factor(physeq_meta_filtered[[my_column]], my_column_ordered)
  
  physeq_meta_filtered[[ml]] <- factor(physeq_meta_filtered[[ml]])
  y = tapply(physeq_meta_filtered$overall.mean, physeq_meta_filtered[[ml]], function(y) max(y))
  y = sort(y, TRUE)
  physeq_meta_filtered[[ml]] = factor(as.character(physeq_meta_filtered[[ml]]), levels=names(y))
  levels(physeq_meta_filtered[[ml]])
  
  # Plot 
  ggplot(physeq_meta_filtered, aes(x = my_column_ordered, y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~.) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = my_colors) +
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = TRUE, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    ylab("Relative Abundance") +
    xlab(my_column) +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%) in at least 1 sample")) 
  ggsave(paste0("output-taxa/", ml, "BarPlot_", my_column, ".png"), height = 5, width = 4)
}


