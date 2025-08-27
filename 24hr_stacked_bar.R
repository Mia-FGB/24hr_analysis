# Analysis of the 2023 and 2024 24hr data
# To plot phylum level stacked bar plots
# Built using a similar script originally for just the 2023 data

setwd("/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/24_hour")

# Libraries --------------------------------------------------------------------
library(phyloseq)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(forcats)
library(patchwork)
library(grid)

# Theme ------------------------------------------------------------------------
custom_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )
theme_set(custom_theme)

# Colours ----------------------------------------------------------------------
phy_colours <- c(
  Actinomycetota = "#BEE6BE",
  Arthropoda = "#FBF3A1",
  Ascomycota = "#E49599",
  Bacillota = "#CECBD0",
  Bacteroidota = "#B0A7DA",
  Basidiomycota = "#FFD8C4",
  Chordata = "#D96C4E",
  Oomycota = "#CAD066",
  Pseudomonadota = "#DCD4AE",
  Streptophyta = "#FCCDE5",
  Higher_Taxa = "darkgrey"
)

phy_colours_24 <- c(
  Actinomycetota = "#BEE6BE",
  Arthropoda = "#FBF3A1",
  Ascomycota = "#E49599",
  Bacillota = "#CECBD0",
  Bacteroidota = "#B0A7DA",
  Basidiomycota = "#FFD8C4",
  Oomycota = "#CAD066",
  Pseudomonadota = "#DCD4AE",
  Streptophyta = "#FCCDE5",
  Higher_Taxa = "darkgrey"
)

fung_colours <- c(
  Alternaria ="#8DD3C7",
  Epichloe = "#F0F0BB",
  Diaporthe = "#C5B2CC",
  Dioszegia = "#31A465",
  Parastagonospora = "#C6B293",
  Peronospora = "#4B6587",
  Puccinia ="#CFD799",
  Pyrenophora = "#D5B030",
  Ramularia = "#F5A469",
  Septoria =  "#8FCF3C",
  Ustilago = "#A9746E",
  Higher_Taxa = "darkgrey"
)

fung_colours_24 <- c(
  Alternaria ="#8DD3C7",
  Dioszegia = "#31A465",
  Parastagonospora = "#C6B293",
  Peronospora = "#4B6587",
  Puccinia ="#CFD799",
  Ramularia = "#F5A469",
  Ustilago = "#A9746E",
  Higher_Taxa = "darkgrey"
)


# Load in data -----------------------------------------------------------------
otu <- read.delim("taxa_counts/all_24hr_data_assigned.tsv", check.names = FALSE)
otu <- otu[,-c(1,3)]
rownames(otu) <- otu[,1]
otu <- otu[,-1]

taxa <- read.csv("taxa_counts/all_24hr_data_taxaID_lineage.csv")
taxa <- taxa %>% mutate_if(is.character, as.factor)
taxa <- taxa[-c(1), ]
rownames(taxa) <- rownames(otu)
taxa <- taxa[,-1]

meta <- read.delim("metadata/All metadata - 24 hour collections (2023 & 2024).tsv")
meta$Duration_min <- as.integer(meta$Duration_min)
meta$Duration_Hrs <- as.integer(meta$Duration_Hrs)
rownames(meta) <- meta[,2]

otu_mat<- as.matrix(otu)
tax_mat<- as.matrix(taxa)

# Phyloseq ---------------------------------------------------------------------
phylo_OTU<- otu_table(otu_mat, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax_mat)
phylo_samples <- sample_data(meta)
phylo_object <- phyloseq(phylo_OTU, phylo_TAX, phylo_samples)

# Phylum data ------------------------------------------------------------------
samp_phylum <- phylo_object %>%
  tax_glom(taxrank = "phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() %>%
  filter(Abundance > 0.01) %>%
  arrange(phylum)

uniq_phyla <- unique(samp_phylum$phylum)
ordered_phyla <- sort(setdiff(uniq_phyla, "Higher_Taxa"))
ordered_phyla <- c(ordered_phyla, "Higher_Taxa")
samp_phylum$phylum <- factor(samp_phylum$phylum, levels = ordered_phyla)

samp_phylum <- samp_phylum %>% mutate(Duration_Hrs = ifelse(is.na(Duration_Hrs), 0, Duration_Hrs))

# Fungi data -------------------------------------------------------------------
samp_genus <- phylo_object %>%
  tax_glom(taxrank = "genus") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() %>%
  filter(Abundance > 0.02) %>%
  arrange(genus)

fungal_phyla <- c("Opisthosporidia","Chytridiomycota","Neocallimastigomycota",
                  "Blastocladiomycota","Zoopagomycota","Mucoromycota",
                  "Glomeromycota","Basidiomycota","Ascomycota","Oomycota")

fungi_data <- samp_genus %>% filter(phylum %in% fungal_phyla)

uniq_genera <- unique(fungi_data$genus)
ordered_genera <- sort(setdiff(uniq_genera, "Higher_Taxa"))
ordered_genera <- c(ordered_genera, "Higher_Taxa")
fungi_data$genus <- factor(fungi_data$genus, levels = ordered_genera)

fungi_data <- fungi_data %>% mutate(Duration_Hrs = ifelse(is.na(Duration_Hrs), 0, Duration_Hrs))

# Dummy row so May_12_12_3 appears as empty bar in 24h graph -------------------
meta_missing <- meta %>% filter(Name == "May_12_12_3")
dummy_row <- meta_missing %>%
  mutate(
    OTU = "dummy_OTU",
    Abundance = 0,
    Sample = "May_12_12_3",
    kingdom = "Eukaryota",
    phylum = "Ascomycota",
    class = "Dothideomycetes",
    order = "Pleosporales",
    family = "Phaeosphaeriaceae",
    genus = "Parastagonospora"
  )
fungi_data <- bind_rows(fungi_data, dummy_row)

# Functions --------------------------------------------------------------------
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
  
  filtered_data$Name <- factor(filtered_data$Name, levels = unique(filtered_data$Name))
  durations <- unique(filtered_data$Duration_Hrs)
  
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  for (dur in durations) {
    plot_data <- filtered_data %>% filter(Duration_Hrs == dur)
    plot_data$Taxon <- plot_data[[taxonomic_var]]
    
    if (taxonomic_var == "phylum" && exists("ordered_phyla")) {
      plot_data$Taxon <- factor(plot_data$Taxon, levels = ordered_phyla)
    } else if (taxonomic_var == "genus" && exists("ordered_genera")) {
      plot_data$Taxon <- factor(plot_data$Taxon, levels = ordered_genera)
    }
    
    p <- ggplot(plot_data, aes(x = Name, y = Abundance, fill = Taxon)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = fill_colours, drop = FALSE) +          # <- keep shared legend
      labs(fill = fill_label) +
      ylab(ylab_text) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
        axis.title.x = element_blank()
      )
    
    file_name <- paste0(save_path, save_prefix, "_", month, year, "_Duration_", dur, ".svg")
    ggsave(file_name, p, width = 10, height = 6)
    
    plot_list[[paste0(month, "_", year, "_", dur, "h")]] <- p
  }
  
  return(plot_list)
}

# Panel plot for phylum --------------------------------------------------------
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

aug_label  <- wrap_elements(grid::textGrob("August 2023", gp = gpar(fontsize = 12)))
may_label  <- wrap_elements(grid::textGrob("May 2024", gp = gpar(fontsize = 12)))
june_label <- wrap_elements(grid::textGrob("June 2024", gp = gpar(fontsize = 12)))

top_row <- aug_label + may_label + june_label + plot_layout(ncol = 3)

plot_row1 <- plot_list$p_august_4h_2023 +
  plot_list$p_may_2h_2024 +
  plot_list$p_june_2h_2024 +
  plot_layout(ncol = 3)

plot_row2 <- plot_list$p_august_6h_2023 +
  plot_list$p_may_6h_2024 +
  plot_list$p_june_6h_2024 +
  plot_layout(ncol = 3)

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
    plot.margin = margin(t = 4, r = 6, b = 10, l = 10)
  )

print(final_plot)
grid.text("Collection Month and Time", x = 0.5, y = 0.02, gp = gpar(fontsize = 12))
grid.text("Relative Abundance of Phyla (>0.01%)", x = 0.01, y = 0.5, rot = 90, gp = gpar(fontsize = 12))

# Manually export as PDF 14 x 8 Graphs/phylum_panel

# Fungi genus panels -----------------------------------------------------------
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

top_f_row <- aug_label + may_label + june_label + plot_layout(ncol = 3)

plot_f_row1 <- fungi_plot_list$f_august_4h_2023 +
  fungi_plot_list$f_may_2h_2024 +
  fungi_plot_list$f_june_2h_2024 +
  plot_layout(ncol = 3)

plot_f_row2 <- fungi_plot_list$f_august_6h_2023 +
  fungi_plot_list$f_may_6h_2024 +
  fungi_plot_list$f_june_6h_2024 +
  plot_layout(ncol = 3)

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
grid.text("Collection Month and Time", x = 0.5, y = 0.02, gp = gpar(fontsize = 12))
grid.text("Relative Abundance of Fungi & Oomycetes (>0.02%)", x = 0.01, y = 0.5, rot = 90, gp = gpar(fontsize = 12))

# Manually export as PDF 14 x 8 Graphs/phylum_panel_2

# 24 hour plots only -----------------------------------------------------------
plot_single_month_24hr <- function(data, month, year,
                                   taxonomic_var = "phylum",
                                   fill_colours,
                                   fill_label = "Taxon",
                                   ylab_text = "Relative Abundance",
                                   save_prefix = "Taxon_24hr",
                                   save_path = "Graphs/") {
  
  data$Start_Time <- as.POSIXct(data$Start_Time, format = "%Y-%m-%d %H:%M", tz = "GMT")
  
  filtered_data <- data %>%
    filter(Year == year, Month == month, Duration_Hrs == 24, LPM == 200) %>%
    arrange(Start_Time, Repeat)
  
  filtered_data <- filtered_data %>%
    mutate(Sample_Label = paste0(format(Start_Time, "%b %d\n%H:%M"), " - 12:00\nRep ", Repeat))
  
  sample_order <- filtered_data %>%
    distinct(Sample, .keep_all = TRUE) %>%
    arrange(Start_Time, Repeat) %>%
    pull(Sample_Label)
  
  filtered_data$Sample_Label <- factor(filtered_data$Sample_Label, levels = sample_order)
  
  filtered_data$Taxon <- filtered_data[[taxonomic_var]]
  
  if (taxonomic_var == "phylum" && exists("ordered_phyla")) {
    filtered_data$Taxon <- factor(filtered_data$Taxon, levels = ordered_phyla)
  } else if (taxonomic_var == "genus" && exists("ordered_genera")) {
    filtered_data$Taxon <- factor(filtered_data$Taxon, levels = ordered_genera)
  }
  
  p <- ggplot(filtered_data, aes(x = Sample_Label, y = Abundance, fill = Taxon)) +
    geom_bar(stat = "identity") +
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
    ) +
    scale_fill_manual(values = fill_colours, drop = FALSE)   # <- keep shared legend
  
  file_name <- paste0(save_path, save_prefix, "_", month, "_", year, ".pdf")
  ggsave(file_name, p, width = 12, height = 4)
  
  return(p)
}

# 24 hour phylum ---------------------------------------------------------------
may_june_24hr <- samp_phylum %>%
  filter(Year == 2024,
         Month %in% c("May", "June"),
         Duration_Hrs == 24,
         LPM == 200)   # only the 24hr runs

present_phyla <- unique(may_june_24hr$phylum)


may_24hr_plot <- plot_single_month_24hr(
  data = samp_phylum,
  month = "May",
  year = 2024,
  taxonomic_var = "phylum",
  fill_colours = phy_colours_24,
  fill_label = "Phylum",
  ylab_text = "Relative Abundance of Phyla (>0.01%)",
  save_prefix = "Phylum_24hr"
)

june_24hr_plot <- plot_single_month_24hr(
  data = samp_phylum,
  month = "June",
  year = 2024,
  taxonomic_var = "phylum",
  fill_colours = phy_colours_24,
  fill_label = "Phylum",
  ylab_text = "Relative Abundance of Phyla (>0.01%)",
  save_prefix = "Phylum_24hr"
)

may_24hr_plot  <- may_24hr_plot  + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
june_24hr_plot <- june_24hr_plot + theme(axis.title.y = element_blank())

combined_24hr_plot <- may_24hr_plot / june_24hr_plot +
  plot_layout(ncol = 1, heights = c(1, 1), guides = "collect") &
  theme(legend.position = "right",
        plot.margin = margin(t = 4, r = 6, b = 10, l = 10))

combined_24hr_plot
grid.text("Relative Abundance of Phyla (>0.01%)", x = 0.01, y = 0.5, rot = 90, gp = gpar(fontsize = 12))

# Manually export 10 x 8 Graphs/phylum_24hr

# 24 hour fungi ---------------------------------------------------------------
may_june_24hr_fung <- fungi_data %>%
  filter(Year == 2024,
         Month %in% c("May", "June"),
         Duration_Hrs == 24,
         LPM == 200)   # only the 24hr runs

present_fung <- unique(may_june_24hr_fung$genus)

may_24hr_fung <- plot_single_month_24hr(
  data = fungi_data,
  month = "May",
  year = 2024,
  taxonomic_var = "genus",
  fill_colours = fung_colours_24,
  fill_label = "Genus",
  ylab_text = "Relative Abundance of Fungi (>0.02%)",
  save_prefix = "Fungi_Genus_24hr"
)

june_24hr_fung <- plot_single_month_24hr(
  data = fungi_data,
  month = "June",
  year = 2024,
  taxonomic_var = "genus",
  fill_colours = fung_colours_24,
  fill_label = "Genus",
  ylab_text = "Relative Abundance of Fungi (>0.02%)",
  save_prefix = "Fungi_Genus_24hr"
)

may_24hr_fung  <- may_24hr_fung  + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
june_24hr_fung <- june_24hr_fung + theme(axis.title.y = element_blank())

combined_24hr_fung <- may_24hr_fung / june_24hr_fung +
  plot_layout(ncol = 1, heights = c(1, 1), guides = "collect") &
  theme(legend.position = "right",
        plot.margin = margin(t = 4, r = 6, b = 10, l = 10))

combined_24hr_fung
grid.text("Relative Abundance of Fungi & Oomycete Genera (>0.02%)", x = 0.01, y = 0.5, rot = 90, gp = gpar(fontsize = 12))

# Manually export 10 x 8 Graphs/fungi_24hr

# Optional label helper --------------------------------------------------------
# geom_text(aes(
#   label = ifelse(is.na(Start_Time), "N/A", format(Start_Time, "%d")),
#   y = 1.05),
#   angle = 0,
#   size = 3,
#   color = "black")
