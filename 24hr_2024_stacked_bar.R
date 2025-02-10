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

#OTU table - need to generate
otu <- read.csv("", check.names = FALSE)
rownames(otu) <- otu[,1] #Making the rownames the NCBI ID
otu <- otu[,-c(1,2)]          #remove 1st & 2nd col

#taxa table - need to generate
taxa <- read.csv("")
taxa<- taxa %>% 
  mutate_if(is.character, as.factor)

rownames(taxa) <- rownames(otu) #same rownames as the otu table
taxa <- taxa[,-1]               

#sample meta data - need to change this to take a tsv file
meta <- read.csv("../metadata/2024_metadata/24h 2024 metadata - Sheet1.tsv")
rownames(meta) <- meta[,1]       #1st col is rownames

#tables need to be matrices
otu_mat<- as.matrix(otu)
tax_mat<- as.matrix(taxa)