# Script to create radial plots from MARTi read count count data

#Packages ------
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lubridate) # for date time
library(reshape2)  #to melt data
library(RColorBrewer)

# Set plotting theme
theme_set(theme_bw())

#Reading in & Procesisng data ----------

# Read counts ---
# Using summed data so I can work at Genus level 
data <-  read.delim("taxa_counts/all_24hr_data_summed.tsv", check.names = FALSE)

pathogens <- c("Puccinia", "Blumeria", "Fusarium", "Zymoseptoria", "Ustilago", "Magnaporthe",
                               "Claviceps", "Pyrenophora", "Parastagonospora", "Phaeosphaeria")
pathogen_data <- data %>% 
  filter(`NCBI Rank` == 'genus') %>% 
  filter(Name %in% pathogens)

# Metadata ---
meta <- read.delim("metadata/All metadata - 24 hour collections (2023 & 2024).tsv")
# Correct type
meta$Duration_min <- as.integer(meta$Duration_min)
meta$Duration_Hrs <- as.integer(meta$Duration_Hrs)
#convert to date time
meta$Start_Time <- ymd_hm(meta$Start_Time)
meta$End_Time <- ymd_hm(meta$End_Time)
#Rename column
meta <- meta %>%
  rename(Sample = Sample_ID)

#  Read numbers ---
read_no <- read.delim("metadata/all_read_numbers.tsv")
read_no <- read_no  %>%
  rename(Sample = ID)

# Melting & merging the data ------

# Melt data to long 
pathogen_data <- melt(pathogen_data, id.vars = c("Name", "NCBI ID", "NCBI Rank"), 
     variable.name = "Sample", value.name = "Read_Count")

# Merge on metadata & read numbers
pathogen_meta <- left_join(pathogen_data, meta, by = "Sample")
pathogen_meta <- left_join(pathogen_meta, read_no, by = "Sample")

# Calculate HPM and other metrics 
pathogen_meta <- pathogen_meta %>%
  mutate(HPM = (Read_Count / ReadsPassBasecall) * 1000000,
         Mid_Time = Start_Time + (Duration_Hrs/2) * 3600)

# Create Label column so simultaneous samples can be grouped
# Only the 2024 samples were collected at the same time, 2023 should stay the same with repeat as 1
pathogen_meta <- pathogen_meta %>%
  mutate(
    Label = if_else(Year == 2024, 
                    str_extract(Sample, ".*(?=_[0-9]+$)"), 
                    Sample),
    Repeat = if_else(Year == 2024, 
                     str_extract(Sample, "(?<=_)[0-9]+$"), 
                     "1")
  )

pathogen_meta <- pathogen_meta %>% 
  filter(!Duration_Hrs %in% c(0, 12, 24), !is.na(Duration_Hrs)) %>%  # don't want to plot the neg, 12 & 24 samples
  mutate(HPM = ifelse(HPM < 100, 0, HPM))  # Removing low read count values, but keeping as 0 for avg

# Plotting --------

# Function to plot ----
# Define the function to generate the plot
generate_pathogen_plot <- function(month, year, name, data) {
  
  # Filter the data based on the input month, year, and name
  filtered_data <- data %>%
    filter(Month == month,
           Year == year,
           Name == name)
  
  # Check if there is any data with HPM > 0
  if (sum(filtered_data$HPM > 0, na.rm = TRUE) == 0) {
    message("No data with HPM > 0 for ", name, " in ", month, " ", year, ". Exiting function.")
    return()  # Exit the function if there is no HPM > 0 data
  }
  
  # Compute mean HPM for each Sample_Time
  summary_data <- filtered_data %>%
    group_by(Label, Mid_Time, Duration_Hrs) %>%
    summarise(
      Avg_HPM = mean(HPM),
      SD_HPM = sd(HPM),  # Standard deviation
      SE_HPM = SD_HPM / sqrt(n()),  # Standard error
      .groups = "drop"
    ) %>%
    arrange(Mid_Time)  # Sort the data by Mid_Time
  
  # Get min and max times to set limits for the filtered data
  min_time <- min(filtered_data$Start_Time, na.rm = TRUE)
  max_time <- max(filtered_data$End_Time, na.rm = TRUE)
  
  # Create the plot
  pathogen_plot <- ggplot(summary_data, aes(x = Mid_Time, y = Avg_HPM,
                                            width = Duration_Hrs * 3600,
                                            fill = factor(Duration_Hrs),
                                            colour = factor(Duration_Hrs))) +
    geom_bar(stat = "identity", position = "identity", alpha = 0.3) +
    scale_fill_manual(values = c("#496B61", "#F2522C")) +
    scale_color_manual(values = c('#374A4A', '#A83A1E')) +
    geom_errorbar(aes(ymin = Avg_HPM - SE_HPM, ymax = Avg_HPM + SE_HPM),
                  width = 0.4, color = "black") +
    scale_x_datetime(date_labels = "%H:%M\n%d %b",
                     date_breaks = "2 hours",
                     limits = c(min_time, max_time)) +  # Set limits based on the data
    labs(
      title = paste(name, month, year),
      x = "Start Time",
      y = "Average HPM",
      fill = "Sample Duration (Hrs)"
    ) +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black", size = 0.8),  # Make axis lines thicker and black
      panel.border = element_blank()  # Remove the panel border
    ) +
    scale_y_continuous(expand = c(0, 0)) +  # Set the y-axis to start at 0
    guides(colour = "none")
  
  # Save the plot
  filename <- paste0("Graphs/pathogen_graphs/", name, "/", month, "_", year, "_both.svg")
  ggsave(filename, plot = pathogen_plot, device = "svg", bg = "transparent", width = 10, height = 6,)
}

# Example function usage
# generate_pathogen_plot("June", 2024,"Zymoseptoria",pathogen_meta)

# Loop through the pathogens and submit functions
# List of months and years for these experiments
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

