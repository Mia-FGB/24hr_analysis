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

# Merge on metadata 
pathogen_meta <- left_join(pathogen_data, meta, by = "Sample")
pathogen_meta <- left_join(pathogen_meta, read_no, by = "Sample")

# Calculate HPM and other metrics 
pathogen_meta <- pathogen_meta %>%
  mutate(HPM = (Read_Count / ReadsPassBasecall) * 100000,
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

# Compute the midpoint of the time range for each sample
# pathogen_meta <- pathogen_meta %>%
#   mutate(
#     Start_Hour = as.numeric(format(Start_Time, "%H")) + as.numeric(format(Start_Time, "%M")) / 60,
#     End_Hour = as.numeric(format(End_Time, "%H")) + as.numeric(format(End_Time, "%M")) / 60,
#     Mid_Hour = (Start_Hour + End_Hour) / 2  # Midpoint between Start and End Hour
#   )

# Plotting --------
# Test plot with filtered data
filtered_data <- pathogen_meta %>%
  filter(Month == 'June',
         Year == 2024,
         Name == 'Zymoseptoria',
         Duration_Hrs <= 6 )

# Compute mean HPM for each Sample_Time
summary_data <- filtered_data %>%
  group_by(Label, Mid_Time, Duration_Hrs) %>%
  summarise(
    Avg_HPM = mean(HPM),
    SD_HPM = sd(HPM),   # Standard deviation
    SE_HPM = SD_HPM / sqrt(n()),  # Standard error
    .groups = "drop"
  )

# Sort the data by Start_Time
summary_data <- summary_data %>% arrange(Mid_Time)

# Get min and max Start_Time to set limits
min_start_time <- min(pathogen_meta$Start_Time)
max_start_time <- max(pathogen_meta$End_Time)

pathogen_plot <- ggplot(summary_data, aes(x = Mid_Time, y = Avg_HPM, width = Duration_Hrs * 3600, fill = factor(Duration_Hrs))) +
  geom_bar(stat = "identity", position = "identity", colour = 'black', alpha = 0.7) +
  geom_errorbar(aes(ymin = Avg_HPM - SE_HPM, ymax = Avg_HPM + SE_HPM), width = 0.2, color = "black") +
  scale_x_datetime(
    date_labels = "%d-%b\n%H:%M", 
    date_breaks = "2 hours", 
    limits = c(min_start_time, max_start_time)  # Set limits based on the data
  ) + 
  scale_fill_brewer(palette = "Set3", name = "Collection Duration") + 
  labs(
    title = "Zymoseptoria June",
    x = "Start Time",
    y = "Average HPM")
  # ) +
  # coord_polar(start = 0) #To make it radial

ggsave("/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/24_hour/Graphs/pathogen_graphs/Zymoseptoria/2_and_6_plot.png",
       pathogen_plot)
