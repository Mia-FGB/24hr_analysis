#Analysis of the 2024 24hr data, using a similar script to the 2023 data

#Library -----------------------------------------------------------------------
library(phyloseq)   # Facilitate the import, storage, analysis, and graphical display of microbiome census data.
library(vegan)      # Analysis of variance using distance matrices using the adonis2 function
library(Maaslin2)   # Determine multivariable association between metadata and microbial meta-omics features; differential abundance analysis
library(ggplot2)    # Generate visualization plots 
library(ggsignif)   # Visualize comparisons and the significance between two groups
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(forcats)

# Set plotting theme
theme_set(theme_bw())

# Load in data ------------------------------------------------------------

#OTU table --
#generated the assigned only read tsv from MARTi output using Scripts/split_marti_taxa.py 
otu <- read.delim("taxa_counts/2024_data/11_feb_assigned.tsv", check.names = FALSE)
otu <- otu[,-c(1,3)]          # remove Name & Rank
rownames(otu) <- otu[,1]      # Making the rownames the NCBI ID
otu <- otu[,-1]             # Then drop the repeated col 

#Taxa table --
#generated the lineages using - /Scripts/get_lineage_from_marti.sh
taxa <- read.csv("taxa_counts/2024_data/marti_assignments_lca_0.1_all_levels_2025-FEB-11_9-58-14_taxaID_lineage.csv")
taxa<- taxa %>% 
  mutate_if(is.character, as.factor)
#Remove first row (NCBI ID which is empty)
taxa <- taxa[-c(1), ]

rownames(taxa) <- rownames(otu) #same rownames as the otu table
taxa <- taxa[,-1]               

#METADATA --
#Downloaded from google sheet 24hr tab - https://docs.google.com/spreadsheets/d/1Oh7zeWlQewzo9bDmnu5cenVM5a9zddVSlzS5OcUnhaE/edit?gid=196309316#gid=196309316
meta <- read.delim("metadata/All metadata - 24 hour collections (2023 & 2024).tsv")

meta <- meta %>%
  filter(Year == 2024)  # Keep only rows where Year is 2024 for this analysis

meta$Duration_min <- as.integer(meta$Duration_min)
meta$Duration_Hrs <- as.integer(meta$Duration_Hrs)

rownames(meta) <- meta[,2]       #Make the rownames the Sample ID

#tables need to be matrices
otu_mat<- as.matrix(otu)
tax_mat<- as.matrix(taxa)

# Phyloseq ----------------------------------------------------------------

#Make phyloseq object
phylo_OTU<- otu_table(otu_mat, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax_mat)
phylo_samples <- sample_data(meta)
phylo_object_24 <- phyloseq(phylo_OTU, phylo_TAX, phylo_samples) #Bring them together


# Plotting -----------------------------------------------------------------


#Stacked phyla  ----
samp_phylum <- phylo_object_24 %>%
  tax_glom(taxrank = "phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance - compared to pass reads
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa <0.01%
  arrange(phylum)                                      # Sort data frame alphabetically by phylum

#Get unique phylum in this dataset 
uniq_phyla <- unique(samp_phylum$phylum)
#Aplly consistent colours from a colour blind friendly palette
phy_colours <- colorRampPalette(brewer.pal(8, "Set3"))(length(uniq_phyla))
phy_colours <- setNames(phy_colours, uniq_phyla)

#Including N/A's which are negative & Lambda
samp_phylum <- samp_phylum %>% mutate(Duration_Hrs = ifelse(is.na(Duration_Hrs), 0, Duration_Hrs))


#Function to plot stacked phylum for duration, month, year
plot_phylum_by_duration <- function(data, month, year, save_path = "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/24_hour/Graphs/") {
  
  # Ensure Start_Time is in a proper date format (handle NAs if needed)
  data$Start_Time <- as.POSIXct(data$Start_Time, format = "%Y-%m-%d %H:%M", tz = "GMT")
  
  # Filter data by the selected month
  filtered_data <- data %>% filter(Month == month, Year == year)
  
  # Order data by Start_Time (ascending or descending)
  filtered_data <- filtered_data %>% arrange(Start_Time)
  # Ensure Samples are uniquely ordered by Start_Time
  filtered_data$Sample <- factor(filtered_data$Sample, levels = unique(filtered_data$Sample))
  
  # Get unique durations
  durations <- unique(filtered_data$Duration_Hrs)
  
  # Create folder if it doesn't exist
  if (!dir.exists(save_path)) {
    dir.create(save_path)
  }
  
  # Loop through each duration and create a plot
  for (dur in durations) {
    plot_data <- filtered_data %>% filter(Duration_Hrs == dur)
    
    p <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = phylum)) +  
      geom_bar(stat = "identity") +
      geom_text(aes(
        label = ifelse(is.na(Start_Time), "N/A", format(Start_Time, "%d")),
        y = 1.05),  # Place labels slightly above the bars
        angle = 0,  # Keep horizontal
        size = 3,
        color = "black") +
      scale_fill_manual(values = phy_colours) +
      theme(axis.title.x = element_blank()) +
      labs(fill = "Phylum") +
      ylab("Relative Abundance (Phyla > 0.01%) \n") +
      ggtitle(paste("Phylum Composition -", month, year, dur, "(Hrs)")) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
    
    # Save the plot
    file_name <- paste0(save_path, "Phylum_", month, year, "_Duration_", dur, ".svg")
    ggsave(file_name, p, width = 10, height = 6)
    
    print(p)  # Display the plot
  }
}

#Run the function - Will overweite saved graphs
# plot_phylum_by_duration(samp_phylum, "May", 2024)
# plot_phylum_by_duration(samp_phylum, "June", 2024)


# Fungi specific -----

#Filter the dataset to just be fungi
#First want to filter to genus level 
samp_genus <- phylo_object_24 %>%
  tax_glom(taxrank = "genus") %>%                      # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance - compared to pass reads
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa diff threshold to before
  arrange(genus)                                      


fungal_phyla <- c("Opisthosporidia", "Chytridiomycota", "Neocallimastigomycota", 
                  "Blastocladiomycota", "Zoopagomycota", "Mucoromycota", 
                  "Glomeromycota", "Basidiomycota", "Ascomycota")

fungi_data <- samp_genus %>% 
  filter(phylum %in% fungal_phyla)

#Get unique phylum in this dataset 
uniq_genera <- unique(fungi_data$genus)
#Aplly consistent colours from a colour blind friendly palette
fung_colours <- colorRampPalette(brewer.pal(8, "Set3"))(length(uniq_genera))
fung_colours <- setNames(fung_colours, uniq_genera)

#Including N/A's which are negative & Lambda
fungi_data <- fungi_data %>% mutate(Duration_Hrs = ifelse(is.na(Duration_Hrs), 0, Duration_Hrs))

#Function to plot stacked phylum for duration, month, year
plot_fungi_by_duration <- function(data, month, year, save_path = "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/24_hour/Graphs/") {
  
  # Ensure Start_Time is in a proper date format (handle NAs if needed)
  data$Start_Time <- as.POSIXct(data$Start_Time, format = "%Y-%m-%d %H:%M", tz = "GMT")
  
  # Filter data by the selected month
  filtered_data <- data %>% filter(Month == month, Year == year)
  
  # Order data by Start_Time (ascending or descending)
  filtered_data <- filtered_data %>% arrange(Start_Time)
  # Ensure Samples are uniquely ordered by Start_Time
  filtered_data$Sample <- factor(filtered_data$Sample, levels = unique(filtered_data$Sample))
  
  # Get unique durations
  durations <- unique(filtered_data$Duration_Hrs)
  
  # Create folder if it doesn't exist
  if (!dir.exists(save_path)) {
    dir.create(save_path)
  }
  
  # Loop through each duration and create a plot
  for (dur in durations) {
    plot_data <- filtered_data %>% filter(Duration_Hrs == dur)
    
    p <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = genus)) +  
      geom_bar(stat = "identity") +
      geom_text(aes(
        label = ifelse(is.na(Start_Time), "N/A", format(Start_Time, "%d")),
        y = 1.05),  # Place labels slightly above the bars
        angle = 0,  # Keep horizontal
        size = 3,
        color = "black") +
      scale_fill_manual(values = fung_colours) +
      theme(axis.title.x = element_blank()) +
      labs(fill = "Genus") +
      ylab("Relative Abundance (Fungi > 0.02%) \n") +
      ggtitle(paste("Fungal & Oomycete -", month, year, dur, "(Hrs)")) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
    
    # Save the plot
    file_name <- paste0(save_path, "Fungi_Genus_", month, year, "_Duration_", dur, ".svg")
    ggsave(file_name, p, width = 10, height = 6)
    
    print(p)  # Display the plot
  }
}

#Run the function
plot_fungi_by_duration(fungi_data, "May", 2024)
plot_fungi_by_duration(fungi_data, "June", 2024)

