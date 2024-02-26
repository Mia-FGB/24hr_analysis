#Phyloseq analysis of 24hr Cub test using data from MARTi

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
otu <- read.csv("marti_assigned_read_count_24hrCub.csv", check.names = FALSE)
rownames(otu) <- otu[,1] #Making the rownames the NCBI ID
otu <- otu[,-c(1,2)]          #remove 1st & 2nd col

#taxa table - missing species :( 
taxa <- read.csv("24hrCub_taxaIDs_lineage_sep.csv")
taxa<- taxa %>% 
  mutate_if(is.character, as.factor)

rownames(taxa) <- rownames(otu) #same rownames as the otu table
taxa <- taxa[,-1]               

#sample meta data - needs more info
meta <- read.csv("24hr_Cub_meta.csv")
rownames(meta) <- meta[,1]       #1st col is rownames

#tables need to be matrices
otu_mat<- as.matrix(otu)
tax_mat<- as.matrix(taxa)

# Phyloseq ----------------------------------------------------------------


#Make phyloseq object
phylo_OTU<- otu_table(otu_mat, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax_mat)
phylo_samples <- sample_data(meta)
phylo_object<- phyloseq(phylo_OTU, phylo_TAX, phylo_samples) #Bring them together

#Histogram - distribution of read data
sample_sum_df <- data.frame(sum = sample_sums(phylo_object)) #sample_sums = no assigned reads per sample

ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "darkred", binwidth = 100000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") + ylab("Count") +
  theme(axis.title.y = element_blank())

#Community composition ----
samp_phylum <- phylo_object %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance - compared to pass reads
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa <2%
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# Plot 
#phylum_colours <- c('#9e0142', '#d53e4f', '#f46d43', '#fdae61', '#fee08b', '#e6f598', '#abdda4', '#66c2a5', 
#                    '#3288bd', '#5e4fa2', '#b5b3bd')

phylum_colours <- c('#9e0142', '#d53e4f', '#f46d43', '#fee08b', '#abdda4',
                    '#66c2a5', '#3288bd', '#5e4fa2', '#b5b3bd')

org_phy_col <- c('#593225','#89411E','#B74F15','#D96521','#E48A56',
                   '#F0B18B','#F9CEB5','#F9E0D1','#212E22' )

phylum_stacked_bar <- ggplot(
  samp_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colours) +
  theme(axis.title.x = element_blank()) +   # Remove x axis title
  ylab("Relative Abundance (Phyla > 0.01%) \n") +
  ggtitle("Phylum Composition from different collection times") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))

#Filtered plot for poster --------
# Filter out Negative, 12 & 24hr
samp_phylum_filter <- samp_phylum[samp_phylum$Time_hrs < 11 & samp_phylum$Sample != "Negative", ]

#Create time column
#Extract the time range string
samp_phylum_filter$time_range <- sub("^\\d+h_(\\d+:\\d+)_(\\d+:\\d+)$", "\\1 - \\2",samp_phylum_filter$Sample)

# Extract start and end times from the 'time_range' column
start_end_times <- strsplit(samp_phylum_filter$time_range, " - ", fixed = TRUE)
start_times <- sapply(start_end_times, function(x) x[1])

# Convert start and end times to POSIXct format
samp_phylum_filter$start_datetime <- as.POSIXct(paste(samp_phylum_filter$Date_collected, start_times),
                                                format = "%d/%m/%Y %H:%M", tz = "GMT")

# Calculate end datetime using 'time_hrs'
samp_phylum_filter$end_datetime <- samp_phylum_filter$start_datetime + 
  as.difftime(samp_phylum_filter$Time_hrs, units = "hours")

# Clean up
cols_to_keep <- c("OTU", "Sample", "Abundance", "Time_hrs", "Yield",
                  "Basecalled_reads", "Kingdom", "Phylum", "time_range",
                  "start_datetime", "end_datetime")

samp_phylum_filter <- samp_phylum_filter[, cols_to_keep, drop = FALSE]

#Renaming
phylum_kingdom_mapping <- c(
  "Actinomycetota" = "Bacteria",
  "Ascomycota" = "Fungi",
  "Bacillota" = "Bacteria",
  "Bacteroidota" = "Bacteria",
  "Basidiomycota" = "Fungi",
  "Oomycota" = "Oomycete",
  "Pseudomonadota" = "Bacteria",
  "Streptophyta" = "Plant"
)

samp_phylum_filter$Phylum <- paste(samp_phylum_filter$Phylum, " (", phylum_kingdom_mapping[samp_phylum_filter$Phylum], ")", sep = "")

# Order levels of 'Phylum' based on abundance
phylum_order <- samp_phylum_filter %>%
  group_by(Phylum) %>%
  summarize(mean_abundance = mean(Abundance)) %>%
  arrange(mean_abundance) %>%
  pull(Phylum)

# Apply the order to the 'Phylum' factor
samp_phylum_filter$Phylum <- factor(samp_phylum_filter$Phylum, levels = phylum_order)

#Plot
phylum_filter_stacked_bar <- ggplot(
  samp_phylum_filter,
  aes(x = time_range, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = org_phy_col) +
  theme(axis.title.x = element_blank()) +   # Remove x axis title
  ylab("Relative Abundance (Phyla > 0.01%)") +
  #theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  facet_wrap(~Time_hrs, scales = "free_x", ncol = 1,
  labeller = labeller(Time_hrs = c("4" = "4hr Samples", "6" = "6hr Samples")))

phyla <- sort(unique(samp_phylum_filter$Phylum))

#DNA conc with collection time ---------
df_phylo <- psmelt(phylo_object)

#Maybe plot as curve
avg_DNA_conc_time <-  meta %>% 
  group_by(Sampler, Time_hrs) %>% 
  summarise(avg_DNA_conc = mean(DNA_conc))

#Diversity plots------------
#Not sure iI like the look of these / not sure whow they have been calculated

#This method involves subsampling the libraries with replacement to estimate the species abundance of the real population while standardising sampling effort 
min_lib <- min(sample_sums(phylo_object))

#We will subsample to the minimum number of reads.
#We will repeat this 100 times and average the diversity estimates from each trial.

# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(phylo_object)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(phylo_object)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(phylo_object)

# It is always important to set a seed when you subsample so your result is replicable 
set.seed(3)

for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(phylo_object, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
}

#Calculate the Mean and sd per sample for observed richness and inverse simpsons index to store in a df
# Create a new dataframe to hold the means and standard deviations of richness estimates
Sample <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(Sample, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
Sample <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(Sample, mean, sd, measure)

#combine richness and eveness into one df
alpha <- rbind(rich_stats, even_stats)

#add sample metadata using merge 
s <- data.frame(sample_data(phylo_object))
alphadiv <- merge(alpha, s, by = "Sample") 

#Plot alpha diversity measures for different samplers in a facet
ggplot(alphadiv, aes(x = Sample, y = mean, color = Time_hrs, group = Time_hrs)) +
  geom_point(size = 3) + 
  facet_wrap(~measure, ncol = 1, scales = "free") +
  #scale_color_manual(values = c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231')) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

