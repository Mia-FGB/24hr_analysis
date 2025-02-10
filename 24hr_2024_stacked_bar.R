#Analysis of the 2024 24hr data, using a similar script to the 2023 data

#Library -----------------------------------------------------------------------
library(phyloseq)   # Facilitate the import, storage, analysis, and graphical display of microbiome census data.
library(vegan)      # Analysis of variance using distance matrices using the adonis2 function
library(Maaslin2)   # Determine multivariable association between metadata and microbial meta-omics features; differential abundance analysis
library(ggplot2)    # Generate visualization plots 
library(ggsignif)   # Visualize comparisons and the significance between two groups
library(dplyr)
library(RColorBrewer)
library(BiocManager)
library(Biostrings)
library(knitr)
library(scater)
library(patchwork)
library(stringr)
library(fantaxtic)
library(data.table)
library(tidyr)

# Set plotting theme
theme_set(theme_bw())

# Load in data ------------------------------------------------------------

#OTU table 
otu <- read.delim("taxa_counts/2024_data/10_feb_assigned.tsv", check.names = FALSE)
otu <- otu[,-c(1,3)]          # remove Name & Rank
rownames(otu) <- otu[,1]      # Making the rownames the NCBI ID
otu <- otu[,-1]             # Then drop the repeated col 

#taxa table - need to generate
taxa <- read.csv("taxa_counts/2024_data/marti_assignments_noML_lca_0.1_all_levels_2025-FEB-10_10-36-18_taxaID_lineage.csv")
taxa<- taxa %>% 
  mutate_if(is.character, as.factor)
#Remove first row (NCBI ID which is empty)
taxa <- taxa[-c(1), ]

rownames(taxa) <- rownames(otu) #same rownames as the otu table
taxa <- taxa[,-1]               

#sample metadata
meta <- read.delim("metadata/All metadata - 24 hour collections (2023 & 2024).tsv")

meta <- meta %>%
  filter(Year == 2024) %>%  # Keep only rows where Year is 2024
  select(-Year)  # Drop the Year column

rownames(meta) <- meta[,1]       #Make the rownames the Sample ID

#tables need to be matrices
otu_mat<- as.matrix(otu)
tax_mat<- as.matrix(taxa)

# Phyloseq ----------------------------------------------------------------

#Make phyloseq object
phylo_OTU<- otu_table(otu_mat, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax_mat)
phylo_samples <- sample_data(meta)
phylo_object<- phyloseq(phylo_OTU, phylo_TAX, phylo_samples) #Bring them together


#Community composition ----
samp_phylum <- phylo_object %>%
  tax_glom(taxrank = "phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance - compared to pass reads
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa <2%
  arrange(phylum)                                      # Sort data frame alphabetically by phylum

phy_col <- c('#593225','#89411E','#B74F15','#D96521','#E48A56',
             '#F0B18B','#F9CEB5','#F9E0D1','#212E22' ) #Need more colours to use this (11)

phylum_stacked_bar <- ggplot(
  samp_phylum, aes(x = Sample, y = Abundance, fill = phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set3") +
  # scale_fill_manual(values = phy_col) +
  # facet_grid(rows = 'Month') + #Not working as expected - may need to plot sep
  theme(axis.title.x = element_blank()) +   # Remove x axis title
  ylab("Relative Abundance (Phyla > 0.01%) \n") +
  ggtitle("Phylum Composition from different collection times") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
