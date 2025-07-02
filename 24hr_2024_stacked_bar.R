# Analysis of the 2023 and 2024 24hr data
# To plot phylum level stacked bar plots
# Built using a similar script originally for just the 2023 data

setwd("/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/24_hour")


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
library(patchwork)
library(grid)

# Set plotting theme
custom_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
  )

theme_set(custom_theme)

# Load in data ------------------------------------------------------------

#OTU table --
#generated the assigned only read tsv from MARTi output using Scripts/split_marti_taxa.py 
# Use assigned reads for phyloseq as builds the objevt using lineage data
otu <- read.delim("taxa_counts/all_24hr_data_assigned.tsv", check.names = FALSE) #2023 & 2024 samples
otu <- otu[,-c(1,3)]          # remove Name & Rank
rownames(otu) <- otu[,1]      # Making the rownames the NCBI ID
otu <- otu[,-1]             # Then drop the repeated col 

#Taxa table --
#generated the lineages using - /Scripts/get_lineage_from_marti.sh from the marti output
taxa <- read.csv("taxa_counts/all_24hr_data_taxaID_lineage.csv")
taxa<- taxa %>% 
  mutate_if(is.character, as.factor)
#Remove first row (NCBI ID which is empty)
taxa <- taxa[-c(1), ]

rownames(taxa) <- rownames(otu) #same rownames as the otu table
taxa <- taxa[,-1]               

#METADATA --
#Downloaded from google sheet 24hr tab
# https://docs.google.com/spreadsheets/d/1Oh7zeWlQewzo9bDmnu5cenVM5a9zddVSlzS5OcUnhaE/edit?gid=196309316#gid=196309316
meta <- read.delim("metadata/All metadata - 24 hour collections (2023 & 2024).tsv")

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
phylo_object <- phyloseq(phylo_OTU, phylo_TAX, phylo_samples) #Bring them together


# Phylum data -----------------------------------------------------------------


#Collapse to phylum level and filter on abundance  ----
samp_phylum <- phylo_object %>%
  tax_glom(taxrank = "phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance - compared to pass reads
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa >0.01%
  arrange(phylum)                                      # Sort data frame alphabetically by phylum

#Get unique phylum in this dataset 
uniq_phyla <- unique(samp_phylum$phylum)
# Get unique phylum names, ensuring 'Higher_Taxa' is last
ordered_phyla <- sort(setdiff(uniq_phyla, "Higher_Taxa"))  # Alphabetical order
ordered_phyla <- c(ordered_phyla, "Higher_Taxa")  # Append 'Higher_Taxa' at the end

# Assign colors, keeping the existing ones but setting 'Higher_Taxa' to grey
phy_colours <- colorRampPalette(brewer.pal(8, "Set3"))(length(ordered_phyla) - 1)  # Exclude 'Higher_Taxa'
phy_colours <- c(phy_colours, "darkgrey")  # Add grey for 'Higher_Taxa'
phy_colours <- setNames(phy_colours, ordered_phyla)

# Ensure phylum column is a factor with ordered levels
samp_phylum$phylum <- factor(samp_phylum$phylum, levels = ordered_phyla)

#Including N/A's which are negative & Lambda
samp_phylum <- samp_phylum %>% mutate(Duration_Hrs = ifelse(is.na(Duration_Hrs), 0, Duration_Hrs))



# Fungi data  ---------

#Filter the data to just be fungi
#First want to filter to genus level - takes a little while to run
samp_genus <- phylo_object %>%
  tax_glom(taxrank = "genus") %>%                      # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance - compared to pass reads
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa diff threshold to before
  arrange(genus)                                      

#Fungi & Oomycetes being selected
fungal_phyla <- c("Opisthosporidia", "Chytridiomycota", "Neocallimastigomycota", 
                  "Blastocladiomycota", "Zoopagomycota", "Mucoromycota", 
                  "Glomeromycota", "Basidiomycota", "Ascomycota", 'Oomycota')

fungi_data <- samp_genus %>% 
  filter(phylum %in% fungal_phyla)

#Get unique phylum in this dataset 
uniq_genera <- unique(fungi_data$genus)
# Get unique phylum names, ensuring 'Higher_Taxa' is last
ordered_genera <- sort(setdiff(uniq_genera, "Higher_Taxa"))  # Alphabetical order
ordered_genera <- c(ordered_genera, "Higher_Taxa")  # Append 'Higher_Taxa' at the end

# Assign colors, keeping the existing ones but setting 'Higher_Taxa' to grey
fung_colours <- colorRampPalette(brewer.pal(8, "Set3"))(length(ordered_genera) - 1)  # Exclude 'Higher_Taxa'
fung_colours <- c(fung_colours, "darkgrey")  # Add grey for 'Higher_Taxa'
fung_colours <- setNames(fung_colours, ordered_genera)

# Ensure genus is a factor with correct ordering
fungi_data$genus <- factor(fungi_data$genus, levels = ordered_genera)

#Including N/A's which are negative & Lambda
fungi_data <- fungi_data %>% mutate(Duration_Hrs = ifelse(is.na(Duration_Hrs), 0, Duration_Hrs))



# Plotting Functions --------

# Updated function that can be used for phylum or fungi genus 
plot_taxa_by_duration <- function(data, month, year,
                                  taxonomic_var,
                                  fill_colours,
                                  fill_label,
                                  ylab_text = "Relative Abundance \n",
                                  save_prefix = "Taxa",
                                  save_path = "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/24_hour/Graphs/") {
  
  plot_list <- list()
  data$Start_Time <- as.POSIXct(data$Start_Time, format = "%Y-%m-%d %H:%M", tz = "GMT")
  
  filtered_data <- data %>%
    filter(Month == month, Year == year) %>%
    arrange(Start_Time, Repeat)
  
  # Reorder 'Name' by Start_Time and Repeat
  filtered_data$Name <- factor(filtered_data$Name, levels = unique(filtered_data$Name))
  
  durations <- unique(filtered_data$Duration_Hrs)
  
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  for (dur in durations) {
    plot_data <- filtered_data %>% filter(Duration_Hrs == dur)
    
    # Dynamically assign taxonomic column to 'Taxon'
    plot_data$Taxon <- plot_data[[taxonomic_var]]
    
    # Optional ordering
    if (taxonomic_var == "phylum" && exists("ordered_phyla")) {
      plot_data$Taxon <- factor(plot_data$Taxon, levels = ordered_phyla)
    }
    
    p <- ggplot(plot_data, aes(x = Name, y = Abundance, fill = Taxon)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = fill_colours, drop = FALSE) +
      labs(fill = fill_label) +
      ylab(ylab_text) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
            axis.title.x = element_blank())
    
    file_name <- paste0(save_path, save_prefix, "_", month, year, "_Duration_", dur, ".svg")
    ggsave(file_name, p, width = 10, height = 6)
    
    plot_list[[paste0(month, "_", year, "_", dur, "h")]] <- p
  }
  
  return(plot_list)
  View(plot_list)
}

# Call the function and create panel plot for phylum --------

month_years <- c("May" = 2024, "June" = 2024, "August" = 2023)
plot_results <- list()

for (m in names(month_years)) {
  y <- month_years[[m]]
  plot_results[[m]] <- plot_taxa_by_duration(
    data = samp_phylum,
    month = m,
    year = y,
    taxonomic_var = "phylum",
    fill_colours = phy_colours,
    fill_label = "Phylum",
    ylab_text = "Relative Abundance (Phyla > 0.01%) \n",
    save_prefix = "Phylum"
  )
}

plot_list <- list(
  p_august_4h_2023 = plot_results[["August"]][["August_2023_4h"]],
  p_august_6h_2023 = plot_results[["August"]][["August_2023_6h"]],
  p_may_2h_2024    = plot_results[["May"]][["May_2024_2h"]],
  p_may_6h_2024    = plot_results[["May"]][["May_2024_6h"]],
  p_june_2h_2024   = plot_results[["June"]][["June_2024_2h"]],
  p_june_6h_2024   = plot_results[["June"]][["June_2024_6h"]]
)

# Create dummy panels with month labels using wrap_elements
aug_label <- wrap_elements(grid::textGrob("August 2023", gp = gpar(fontsize = 12)))
may_label <- wrap_elements(grid::textGrob("May 2024", gp = gpar(fontsize = 12)))
june_label <- wrap_elements(grid::textGrob("June 2024", gp = gpar(fontsize = 12)))

# Row for top labels
top_row <- aug_label + may_label + june_label + plot_layout(ncol = 3)

# First row of plots
plot_row1 <- plot_list$p_august_4h_2023 + 
  plot_list$p_may_2h_2024 + 
  plot_list$p_june_2h_2024 + 
  plot_layout(ncol = 3)

# Second row of plots
plot_row2 <- plot_list$p_august_6h_2023 + 
  plot_list$p_may_6h_2024 + 
  plot_list$p_june_6h_2024 + 
  plot_layout(ncol = 3)

# Combine all
final_plot <- top_row / plot_row1 / plot_row2 +
  plot_layout(heights = c(0.1, 1, 1), guides = "collect") & 
  theme_minimal(base_size = 12) &
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "plain"),
    plot.margin = margin(t = 4, r = 6, b = 10, l = 10)  # Increase margins
  )

print(final_plot)

# To save with these axis need to export manually as a PDF 12 x 8 inches phylum_panel
grid.text("Collection Month and Time", x = 0.5, y = 0.02, gp = gpar(fontsize = 12))
grid.text("Relative Abundance of Phyla (>0.01%)", x = 0.01, y = 0.5, rot = 90, gp = gpar(fontsize = 12))


# Calling the function on the fungi data and creating panel -----
# Same months years as created earlier
fungi_plot_results <- list()

for (m in names(month_years)) {
  y <- month_years[[m]]
  fungi_plot_results[[m]] <- plot_taxa_by_duration(
    data = fungi_data,
    month = m,
    year = y,
    taxonomic_var = "genus",
    fill_colours = fung_colours,
    fill_label = "Genus",
    ylab_text = "Relative Abundance (Fungi & Oomycete > 0.02%) \n",
    save_prefix = "Fungi_Genus"
  )
}

fungi_plot_list <- list(
  f_august_4h_2023 = fungi_plot_results[["August"]][["August_2023_4h"]],
  f_august_6h_2023 = fungi_plot_results[["August"]][["August_2023_6h"]],
  f_may_2h_2024    = fungi_plot_results[["May"]][["May_2024_2h"]],
  f_may_6h_2024    = fungi_plot_results[["May"]][["May_2024_6h"]],
  f_june_2h_2024   = fungi_plot_results[["June"]][["June_2024_2h"]],
  f_june_6h_2024   = fungi_plot_results[["June"]][["June_2024_6h"]]
)

# Use same labels as with the phylum plot
top_f_row <- aug_label + may_label + june_label + plot_layout(ncol = 3)

# Plot rows
plot_f_row1 <- fungi_plot_list$f_august_4h_2023 + 
  fungi_plot_list$f_may_2h_2024 + 
  fungi_plot_list$f_june_2h_2024 + 
  plot_layout(ncol = 3)

plot_f_row2 <- fungi_plot_list$f_august_6h_2023 + 
  fungi_plot_list$f_may_6h_2024 + 
  fungi_plot_list$f_june_6h_2024 + 
  plot_layout(ncol = 3)

# Combine all
final_fungi_plot <- top_f_row / plot_f_row1 / plot_f_row2 +
  plot_layout(heights = c(0.1, 1, 1), guides = "collect") & 
  theme_minimal(base_size = 12) &
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "plain"),
    plot.margin = margin(t = 4, r = 6, b = 10, l = 10)
  )

print(final_fungi_plot)

# To save with these axis need to export manually 12 x 8 inches phylum_panel
grid.text("Collection Month and Time", x = 0.5, y = 0.02, gp = gpar(fontsize = 12))
grid.text("Relative Abundance of Fungi & Oomycetes (>0.02%)", x = 0.01, y = 0.5, rot = 90, gp = gpar(fontsize = 12))


# Just the 24 hour samples ----------
# Filtering the data to only have 200 LPM samples for now 

plot_single_month_24hr <- function(data, month, year,
                                   taxonomic_var = "phylum",
                                   fill_colours,
                                   fill_label = "Taxon",
                                   ylab_text = "Relative Abundance",
                                   save_prefix = "Taxon_24hr",
                                   save_path = "Graphs/") {
  
  # Ensure Start_Time is a proper datetime
  data$Start_Time <- as.POSIXct(data$Start_Time, format = "%Y-%m-%d %H:%M", tz = "GMT")
  
  # Filter for 24-hour samples from the selected month/year
  filtered_data <- data %>%
    filter(Year == year, Month == month, Duration_Hrs == 24, LPM == 200) %>%
    arrange(Start_Time, Repeat)
  
  # Create informative sample labels
  filtered_data <- filtered_data %>%
    mutate(
      Sample_Label = paste0(format(Start_Time, "%b %d\n%H:%M"), " - 12:00\nRep ", Repeat)
    )
  
  # Order Sample_Label by Start_Time and Repeat
  sample_order <- filtered_data %>%
    distinct(Sample, .keep_all = TRUE) %>%
    arrange(Start_Time, Repeat) %>%
    pull(Sample_Label)
  
  filtered_data$Sample_Label <- factor(filtered_data$Sample_Label, levels = sample_order)
  
  # Assign the taxonomic variable dynamically
  filtered_data$Taxon <- filtered_data[[taxonomic_var]]
  
  # Optional: enforce order (e.g. phylum level)
  if (taxonomic_var == "phylum" && exists("ordered_phyla")) {
    filtered_data$Taxon <- factor(filtered_data$Taxon, levels = ordered_phyla)
  }
  
  # Plot
  p <- ggplot(filtered_data, aes(x = Sample_Label, y = Abundance, fill = Taxon)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = fill_colours, drop = FALSE) +
    labs(
      x = "Sample Collection Time and Replicate",
      y = ylab_text,
      fill = fill_label
    ) +
    ggtitle(paste(month, year)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.line = element_line(color = "black", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Save
  file_name <- paste0(save_path, save_prefix, "_", month, "_", year, ".pdf")
  ggsave(file_name, p, width = 12, height = 4)
  
  return(p)
}


# Run the function for each month phylum 
may_24hr_plot <- plot_single_month_24hr(
  data = samp_phylum,
  month = "May",
  year = 2024,
  taxonomic_var = "phylum",
  fill_colours = phy_colours,
  fill_label = "Phylum",
  ylab_text = "Relative Abundance of Phyla (>0.01%)",
  save_prefix = "Phylum_24hr"
)

june_24hr_plot <- plot_single_month_24hr(
  data = samp_phylum,
  month = "May",
  year = 2024,
  taxonomic_var = "phylum",
  fill_colours = phy_colours,
  fill_label = "Phylum",
  ylab_text = "Relative Abundance of Phyla (>0.01%)",
  save_prefix = "Phylum_24hr"
)



# Remove axis labels
may_24hr_plot <- may_24hr_plot +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
june_24hr_plot <- june_24hr_plot +
  theme(axis.title.y = element_blank())

# Combine the months 
combined_24hr_plot <- may_24hr_plot / june_24hr_plot +
  plot_layout(ncol = 1, heights = c(1, 1), guides = "collect") &
  theme(
    legend.position = "right",
    plot.margin = margin(t = 4, r = 6, b = 10, l = 10)
  )

# Add the y axis centred - then manually export as a PDF 12 x 8
grid.text("Relative Abundance of Phyla (>0.01%)", x = 0.01, y = 0.5, rot = 90, gp = gpar(fontsize = 12))


# Optional code ----

# To add labels of the date to the bars:
# geom_text(aes(
#   label = ifelse(is.na(Start_Time), "N/A", format(Start_Time, "%d")),
#   y = 1.05),  # Place labels slightly above the bars
#   angle = 0,  # Keep horizontal
#   size = 3,
#   color = "black") 