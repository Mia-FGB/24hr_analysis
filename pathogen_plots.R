# Script to create pathogen plots from MARTi read count count data

#Packages ------
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lubridate) # for date time
library(reshape2)  #to melt data
library(RColorBrewer)
library(patchwork)

# Set plotting theme
theme_set(theme_bw())

#Reading in & Processing data ----------

# Read counts ---
# Using summed data so I can work at Genus level 
data <-  read.delim("taxa_counts/all_24hr_data_summed.tsv", check.names = FALSE)

pathogens <- c("Puccinia", "Blumeria", "Fusarium", "Zymoseptoria", "Ustilago", "Magnaporthe",
                               "Claviceps", "Pyrenophora", "Parastagonospora", "Phaeosphaeria")
pathogen_data <- data %>% 
  filter(`NCBI Rank` == 'genus') %>% 
  filter(Name %in% pathogens)

# Melt data to long 
pathogen_data <- melt(pathogen_data, id.vars = c("Name", "NCBI ID", "NCBI Rank"), 
                      variable.name = "Sample", value.name = "Read_Count")

# Metadata ---
meta <- read.delim("metadata/All metadata - 24 hour collections (2023 & 2024).tsv")
# Correct type
meta$Duration_min <- as.integer(meta$Duration_min)
meta$Duration_Hrs <- as.integer(meta$Duration_Hrs)
#convert to date time
meta$Start_Time <- ymd_hm(meta$Start_Time)
meta$End_Time <- ymd_hm(meta$End_Time)
#Rename columns for merging
meta <- meta %>%
  rename(
    Sample = Sample_ID,
    Label = Name
  )

# Split the 'Label' text to drop the repeat
meta <- meta %>%
  mutate(
    Label = ifelse(
      Year == 2024,
      # If Year is 2024: strip trailing digits from Sample
      sub("_(\\d+)$", "", Sample),
      # If Year is not 2024: process only when Start and End times are not NA
      ifelse(
        !is.na(Start_Time) & !is.na(End_Time),
        paste0(
          format(Start_Time, "%B"), "_",
          format(Start_Time, "%H"), "_",
          format(End_Time, "%H")
        ),
        NA_character_  # If no start time then NA
      )
    )
  )

#  Read numbers ---
read_no <- read.delim("metadata/all_read_numbers.tsv")
read_no <- read_no  %>%
  rename(Sample = ID)

# Merging the data ------

# Merge on metadata & read numbers
pathogen_meta <- left_join(pathogen_data, meta, by = "Sample")
pathogen_meta <- left_join(pathogen_meta, read_no, by = "Sample")

# Calculate Hits_per_100k and other metrics 
pathogen_meta <- pathogen_meta %>%
  mutate(HPM = (Read_Count / ReadsPassBasecall) * 1000000,
         Hits_per_100k =  (Read_Count / ReadsPassBasecall) * 100000,
         Mid_Time = Start_Time + (Duration_Hrs/2) * 3600)

pathogen_meta <- pathogen_meta %>% 
  filter(!Duration_Hrs %in% c(0, 12, 24), !is.na(Duration_Hrs)) %>%  # don't want to plot the neg, 12 & 24 samples
  mutate(HPM = ifelse(HPM < 100, 0, HPM))  # Removing low read count values, but keeping as 0 for avg

# Plotting --------

fill_colours <- c(6 = "#496B61", 4 = "#BF99A3", 2 = "#F2522C")
line_colours <- c(6 = "#374A4A", 4 = "#80666D", 2 = "#A83A1E")


# Function to plot ----
# Define the function to generate the plot
generate_pathogen_plot <- function(month, year, name, data) {
  
  # Filter data
  filtered_data <- data %>%
    filter(Month == month,
           Year == year,
           Name == name)
  
    # Summary stats
    summary_data <- filtered_data %>%
      group_by(Label, Mid_Time, Duration_Hrs) %>%
      summarise(
        Avg_Hits_per_100k = mean(Hits_per_100k),
        SD_Hits_per_100k = sd(Hits_per_100k),
        SE_Hits_per_100k = SD_Hits_per_100k / sqrt(n()),
        .groups = "drop"
      ) %>%
      arrange(Mid_Time) %>% 
      # When this has all the levels the legend will contain them all - needed for the panels 
      mutate(Duration_Factor = factor(Duration_Hrs, levels = c(2, 4, 6)))
    
    min_time <- min(filtered_data$Start_Time, na.rm = TRUE)
    max_time <- max(filtered_data$End_Time, na.rm = TRUE)
    
    # plot
    pathogen_plot <- ggplot(summary_data, aes(
      x = Mid_Time, y = Avg_Hits_per_100k,
      width = Duration_Hrs * 3600,
      fill = Duration_Factor,
      colour = Duration_Factor)) +
      geom_bar(stat = "identity", position = "identity", alpha = 0.3) +
      geom_errorbar(aes(
        ymin = Avg_Hits_per_100k - SE_Hits_per_100k,
        ymax = Avg_Hits_per_100k + SE_Hits_per_100k),
        width = 0.4, color = "black") +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_datetime(
        date_labels = "%H:%M\n%d %b",
        date_breaks = "4 hours",
        limits = c(min_time, max_time)) +
      labs(
        title = paste(month, year),
        x = "Date and Time",
        y = "Average Hits per 100k",
        fill = "Sample Duration (Hrs)") +
      theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.8),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
      guides(colour = "none")
    
    # drop is false means all the plots will have all 3 colours in the legend
    pathogen_plot <- pathogen_plot +       
      scale_fill_manual(values = fill_colours, drop = FALSE) +
      scale_color_manual(values = line_colours, drop = FALSE)
  
  # Ensure directory exists
  out_dir <- paste0("Graphs/pathogen_graphs/", name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save the plot
  filename <- paste0(out_dir, "/", month, "_", year, "_both.svg")
  ggsave(filename, plot = pathogen_plot, device = "svg", bg = "transparent", width = 10, height = 6)
  
  return(pathogen_plot)
}


# Example function usage
# generate_pathogen_plot("June", 2024,"Zymoseptoria",pathogen_meta)

# To save the single plots ----
months_years <- list(c("June", 2024), c("May", 2024), c("August", 2023))

# Loop through each month/year combination and each pathogen name
for (month_year in months_years) {
  month <- month_year[[1]]
  year <- month_year[[2]]
  
  # Loop through each pathogen name
  for (name in pathogens) {
    # Call the function with the current month, year, and pathogen name
    generate_pathogen_plot(month = month, 
                           year = year, 
                           name = name, 
                           data = pathogen_meta)
  }
}

# Generate panelled plots for each pathogen 
dir.create("Graphs/pathogen_panels", recursive = TRUE, showWarnings = FALSE)

panel_order <- list(
  list("August", 2023),
  list("May", 2024),
  list("June", 2024)
)

for (name in pathogens) {
  plots <- list()
  
  for (entry in panel_order) {
    month <- entry[[1]]
    year <- entry[[2]]
    
    p <- generate_pathogen_plot(month, year, name, pathogen_meta)
    plots[[paste(month, year)]] <- p
  }
  
  combined_plot <- (plots[[1]] | plots[[2]] | plots[[3]]) +
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom")
  
  ggsave(filename = paste0("Graphs/pathogen_panels/", name, "_panel.svg"),
         plot = combined_plot, width = 18, height = 6, bg = "transparent")
}

# Old plot aesthetics code ----------

# To try and have the labels match the other graphs - Since there are two time points being plotted it doesn't work
# scale_x_datetime(
#   breaks = summary_data$Mid_Time,
#   labels = summary_data$Label,
#   limits = c(min_time, max_time)
# ) +