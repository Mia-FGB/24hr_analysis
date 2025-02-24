#Script to plot bar charts showing number of species & proportion that are unique to that sample collection

#Packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(scales)
library(lubridate)

#-------------------------------------------------------------------------

# 2024 --------------------------------------------------------------------

#-------------------------------------------------------------------------

#Read in the data ----

#All metadata
meta <- read.csv("metadata/All metadata - 24 hour collections (2023 & 2024).tsv", sep ='\t')

meta$Duration_min <- as.integer(meta$Duration_min)
meta$Duration_Hrs <- as.integer(meta$Duration_Hrs)

meta <- meta %>%
  rename(Sample = Sample_ID)

# Split the 'Sample' column into two new columns 'Month_Time' and 'Repeat', while keeping the original 'Sample' column
meta <- meta %>%
  separate(Sample, into = c("Month_Time", "Repeat"), sep = "(?<=\\d)_(?=\\d+$)", extra = "drop", remove = FALSE) 

#convert to date time
meta$Start_Time <- ymd_hm(meta$Start_Time)
meta$End_Time <- ymd_hm(meta$End_Time)

meta <- meta %>% arrange(Start_Time)  

#Load in summarised data from MARTi 
# Generated the summed only read tsv from MARTi output using Scripts/split_marti_taxa.py 
marti_sum <- read.csv("taxa_counts/2024_data/11_feb_summed.tsv", sep ='\t')

marti_sum <- marti_sum %>% 
  rename(Rank = NCBI.Rank)


#Function to extract Unique & Shared species grouping repeats together ---------
process_and_plot_species <- function(year, month, duration) {
  # Filter data for specific year, month, duration, and species rank
  filtered_data <- marti_meta %>%
    filter(Duration_Hrs == duration, Month == month, Year == year, Rank == "species")
  
  # Convert to presence/absence table
  spe_presence_data <- filtered_data %>%
    mutate(`read count` = ifelse(`read count` > 0, 1, 0)) %>%
    select(Month_Time, Name, `read count`) %>%
    pivot_wider(names_from = Month_Time, values_from = `read count`, values_fn = list(`read count` = max))
  
  # Identify unique species (rows summing to 1)
  unique_spe <- spe_presence_data$Name[rowSums(spe_presence_data[, -1]) == 1]
  just_unique_spe <- spe_presence_data[spe_presence_data$Name %in% unique_spe, ]
  
  # Convert back to long format
  unique_sample_name <- melt(just_unique_spe, id.vars = "Name") %>%
    filter(value == 1)
  
  # Count unique species per sample
  num_uniq_spe <- unique_sample_name %>%
    group_by(variable) %>%
    summarise(count = n())
  
  # Ensure all samples are represented
  full_variables <- unique(filtered_data$Month_Time)
  missing_samples <- data.frame(variable = full_variables, count = 0)
  num_uniq_spe_comp <- merge(missing_samples, num_uniq_spe, by = "variable", all.x = TRUE) %>%
    mutate(count.y = ifelse(is.na(count.y), 0, count.y)) %>%
    select(variable, count.y) %>%
    rename(Sample = variable, Unique_All_Count = count.y)
  
  # Total species per sample
  numeric_col_presence_data <- spe_presence_data[, 2:ncol(spe_presence_data)]
  total_species_counts <- data.frame(Sample = colnames(numeric_col_presence_data),
                                     Total_Species_Count = colSums(numeric_col_presence_data))
  
  # Merge with metadata
  taxa_counts <- inner_join(num_uniq_spe_comp, total_species_counts, by = "Sample") %>%
    rename(Month_Time = Sample) %>%
    left_join(meta, by = "Month_Time") %>%
    filter(Repeat != 2) %>%
    arrange(Start_Time)
  
  # Create the plot
  uniq_species_plot <- ggplot(taxa_counts, aes(x = Month_Time)) +
    geom_bar(aes(y = Total_Species_Count, fill = "Total"), stat = "identity", width = 0.82) +
    geom_bar(aes(y = Unique_All_Count, fill = "Unique"), stat = "identity", width = 0.82) +
    labs(title = paste("Total and Unique Species per Sample -", month, duration, "hrs (", year, ")"),
         y = "Number of Species", x = "Collection Time (Hrs)", fill = "Species Count") +
    scale_fill_manual(values = c("Total" = "#F28241", "Unique" = "#19261C"),
                      labels = c("Total", "Unique")) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.line = element_line(color = "black"),
          strip.text = element_text(size = 12), 
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Save the plot
  filename <- paste0("Graphs/uniq_sp_", month, "_", duration, "h_", year, ".svg")
  ggsave(filename, plot = uniq_species_plot, device = "svg", bg = "transparent", width = 6, height = 4,)
  
  return(uniq_species_plot)
}

#Can run the function like this -
# process_and_plot_species(2024, "June", 6)

years <- c(2023, 2024)  # Define the years of interest

for (year in years) {
  # Ensure months exist for this year
  available_months <- unique(marti_meta$Month[marti_meta$Year == year])
  
  for (month in available_months) {
    # Ensure durations exist for this year & month, and drop NAs
    available_durations <- unique(marti_meta$Duration_Hrs[marti_meta$Year == year & marti_meta$Month == month])
    available_durations <- na.omit(available_durations)  # Remove NAs
    
    # Exclude 24h samples
    available_durations <- available_durations[available_durations != 24]
    
    for (duration in available_durations) {
      print(paste("Processing:", year, month, duration, "hrs"))
      process_and_plot_species(year, month, duration)  # Call your function
    }
  }
}

#This function is not working for 2023 data currently, but I think this is due to the df 
#which should be fixed when it is all exported from MARTi as one

#Ungrouped samples -------------
#Good for 24 hour samples - but also if I want to plot the others not grouped with repeats
process_and_plot_ungrouped_samples <- function(year, month, duration) {
  # Filter data for specific year, month, duration, and species rank
  filtered_data <- marti_meta %>%
    filter(Duration_Hrs == duration, Month == month, Year == year, Rank == "species")
  
  # Convert to presence/absence table
  spe_presence_data <- filtered_data %>%
    mutate(`read count` = ifelse(`read count` > 0, 1, 0)) %>%
    select(Sample, Name, `read count`) %>%
    pivot_wider(names_from = Sample, values_from = `read count`, values_fn = list(`read count` = max))
  
  # Identify unique species (rows summing to 1)
  unique_spe <- spe_presence_data$Name[rowSums(spe_presence_data[, -1]) == 1]
  just_unique_spe <- spe_presence_data[spe_presence_data$Name %in% unique_spe, ]
  
  # Convert back to long format
  unique_sample_name <- melt(just_unique_spe, id.vars = "Name") %>%
    filter(value == 1)
  
  # Count unique species per sample
  num_uniq_spe <- unique_sample_name %>%
    group_by(variable) %>%
    summarise(count = n())
  
  # Ensure all samples are represented
  full_variables <- unique(filtered_data$Sample)
  missing_samples <- data.frame(variable = full_variables, count = 0)
  num_uniq_spe_comp <- merge(missing_samples, num_uniq_spe, by = "variable", all.x = TRUE) %>%
    mutate(count.y = ifelse(is.na(count.y), 0, count.y)) %>%
    select(variable, count.y) %>%
    rename(Sample = variable, Unique_All_Count = count.y)
  
  # Total species per sample
  numeric_col_presence_data <- spe_presence_data[, 2:ncol(spe_presence_data)]
  total_species_counts <- data.frame(Sample = colnames(numeric_col_presence_data),
                                     Total_Species_Count = colSums(numeric_col_presence_data))
  
  # Merge with metadata
  taxa_counts <- inner_join(num_uniq_spe_comp, total_species_counts, by = "Sample") %>%
    left_join(meta, by = "Sample") %>%
    mutate(LPM = as.factor(LPM)) %>%  # Ensure LPM is treated as a factor
    mutate(LPM = factor(LPM, levels = c("50", "100", "200")))  # Order LPM levels explicitly
  
  # Order data by Start_Time 
  taxa_counts <- taxa_counts %>% arrange(Start_Time)
  # Ensure Samples are uniquely ordered by Start_Time
  taxa_counts$Sample <- factor(taxa_counts$Sample, levels = unique(taxa_counts$Sample))
  
  # Create the plot
  uniq_species_plot <- ggplot(taxa_counts, aes(x = Sample)) +
    geom_bar(aes(y = Total_Species_Count, fill = "Total"), stat = "identity", width = 0.82) +
    geom_bar(aes(y = Unique_All_Count, fill = "Unique"), stat = "identity", width = 0.82) +
    
    labs(title = paste("Total and Unique Species per Sample -", month, duration, "hrs (", year, ")"),
         y = "Number of Species", x = "Collection Time (Hrs)", fill = "Species Count") +
   
     scale_fill_manual(values = c("Total" = "#F28241", "Unique" = "#19261C"),
                      labels = c("Total", "Unique")) +
    #Facet by flow rate
    facet_wrap(~LPM, scales = "free_x") +  # Make x-axis free for each facet
    #Adding date labels
    geom_text(aes(
      label = ifelse(is.na(Start_Time), "N/A", format(Start_Time, "%d")),
      y = Total_Species_Count + 1),  # Place labels slightly above the bars
      angle = 0,  # Keep horizontal
      size = 4,
      color = "black") +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.line = element_line(color = "black"),
          strip.text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Save the plot
  filename <- paste0("Graphs/uniq_sp_", month, "_", duration, "h_", year, ".svg")
  ggsave(filename, plot = uniq_species_plot, device = "svg", bg = "transparent", width = 6, height = 4,)
  
  return(uniq_species_plot)
}

process_and_plot_ungrouped_samples(2024, "June", 24)
process_and_plot_ungrouped_samples(2024, "May", 24)

# 2023 data ---------------------------------------------------------------

#Older code from 2023, New code can do this quickre and simpler now

#Marti-------------
#sample meta data - needs more info
meta <- read.csv("metadata/2023_metadata/24hr_Cub_meta.csv")
#Load in summarised data from MARTi - this way can collapse to different taxonomic levels without losing reads
marti_sum <- read.csv("taxa_counts/2023_data/marti_summed_read_count_24hrCub.csv")

colnames(marti_sum)[colnames(marti_sum) == "NCBI.Rank"] <- "Rank"

#For this analysis remove the 12 & 24 hour data but retain the Negative
marti_sum <- marti_sum[, -c(4:6)]
marti_sum <-  marti_sum[, -ncol(marti_sum)]

#Create subsets of data
marti_species <- marti_sum %>% 
  filter(Rank == "species") %>% #filter to retain species only
  subset(select = -NCBI.ID)

marti_4h <- marti_sum %>% 
  select(1:3, starts_with("X4h")) %>% 
  filter(Rank == "species") %>% #filter to retain species only
  subset(select = -NCBI.ID)

marti_6h <- marti_sum %>% 
  select(1:3, starts_with("X6h")) %>% 
  filter(Rank == "species") %>% #filter to retain species only
  subset(select = -NCBI.ID)

#Presence / Absence dataset ----------------------------------------------------
# Convert non-zero counts to 1 for presence, for each taxa level

spe_presence_data <- marti_species 
spe_presence_data[, -c(1,2)] <-
  ifelse(marti_species[,-c(1,2)] > 0, 1, 0) #returns 1 if count > 0, else 0; ignores first 2 cols

spe_presence_data_4h <- marti_4h 
spe_presence_data_4h[, -c(1,2)] <-
  ifelse(marti_4h[,-c(1,2)] > 0, 1, 0) 

spe_presence_data_6h <- marti_6h 
spe_presence_data_6h[, -c(1,2)] <-
  ifelse(marti_6h[,-c(1,2)] > 0, 1, 0)

#Unique taxa count dataframe ----------
unique_spe <- spe_presence_data$Taxon[rowSums(spe_presence_data[, -c(1,2)]) == 1] #Sum row 1 = unique taxa.

just_unique_spe <- #subset df to only contain uniq species
  spe_presence_data[spe_presence_data$Taxon %in% unique_spe, ]

unique_sample_name <- 
  melt(just_unique_spe , id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

#Number of unique species per sample
num_uniq_spe <- unique_sample_name %>% 
  group_by(variable) %>% 
  summarise(count = n()) #no uniq species per sample

#This has removed a number of samples with no unique species, adding them back as follows
full_variables <- marti_sum %>% colnames() %>% .[-c(1:3)] #extracting the colnames from og dataset, ignoring first 3

# Create a data frame with full variables and counts initialized to 0
missing_samples <- data.frame(
  variable = full_variables,
  count = 0
)
# Merge the existing dataframe with the full dataframe
num_uniq_spe_comp <- merge(missing_samples, num_uniq_spe, by = "variable", all.x = TRUE)

# Replace NA counts with 0
num_uniq_spe_comp$count.y[is.na(num_uniq_spe_comp$count.y)] <- 0

num_uniq_spe_comp <- num_uniq_spe_comp %>%
  select(-count.x) %>%  # Remove the named column
  rename("Sample" = "variable",
         "Unique_All_Count" = "count.y") 

#4h---
unique_spe_4h <- spe_presence_data_4h$Taxon[rowSums(spe_presence_data_4h[, -c(1,2)]) == 1]

just_unique_spe_4h <- #subset df to only contain uniq species
  spe_presence_data_4h[spe_presence_data_4h$Taxon %in% unique_spe_4h, ]

unique_sample_name_4h <- 
  melt(just_unique_spe_4h , id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

#Number of unique species per sample
num_uniq_spe_4h <- unique_sample_name_4h %>% 
  group_by(variable) %>% 
  summarise(count = n()) #no uniq species per sample

#only missing the 12:00 - 16:00 sample
num_uniq_spe_4h <- 
  bind_rows(num_uniq_spe_4h, data.frame(variable = "X4h_12.00_16.00", count = 0)) %>% 
  rename("Sample" = "variable",
         "Unique_4h_Count" = "count")

#6h---
unique_spe_6h <- spe_presence_data_6h$Taxon[rowSums(spe_presence_data_6h[, -c(1,2)]) == 1]

just_unique_spe_6h <- #subset df to only contain uniq species
  spe_presence_data_6h[spe_presence_data_6h$Taxon %in% unique_spe_6h, ]

unique_sample_name_6h <- 
  melt(just_unique_spe_6h , id.vars = c("Taxon", "Rank")) %>% # Melt the data to long
  filter(value == 1) #Only keep present taxa

#Number of unique species per sample
num_uniq_spe_6h <- unique_sample_name_6h %>% 
  group_by(variable) %>% 
  summarise(count = n()) #no uniq species per sample
#None are missing

num_uniq_spe_6h <- num_uniq_spe_6h %>% 
  rename("Sample" = "variable",
         "Unique_6h_Count" = "count")

# Total Species in each Sample ---------------
numeric_col_presence_data <- spe_presence_data[, 3:ncol(spe_presence_data)]
total_species_counts <- data.frame(Sample = colnames(numeric_col_presence_data), 
                                   Total_Species_Count = colSums(numeric_col_presence_data))

#Merge with the dataframes
taxa_counts <- inner_join(num_uniq_spe_comp, total_species_counts, by = "Sample")

#Want to merge withe 4h/6h separately calculated unique to their time point
taxa_counts <- left_join(taxa_counts, num_uniq_spe_4h, by = "Sample")
taxa_counts <- left_join(taxa_counts, num_uniq_spe_6h, by = "Sample")

# Sum the 'Unique_4h_Count' and 'Unique_6h_Count' columns
taxa_counts$combined_unique_count <- rowSums(taxa_counts[, c("Unique_4h_Count", "Unique_6h_Count")], na.rm = TRUE)

#Cleaning up the data -----
#Changing the string so I can merge the metadata
taxa_counts$Sample <- sub("^X(\\d+h_\\d+\\.\\d+_\\d+\\.\\d+)$", "\\1", taxa_counts$Sample)
taxa_counts$Sample <- gsub("\\.", ":", taxa_counts$Sample)

taxa_counts$time_range <- sub("^\\d+h_(\\d+:\\d+)_(\\d+:\\d+)$", "\\1 - \\2",taxa_counts$Sample)

#Merging
taxa_count_data <- merge(meta, taxa_counts, by = "Sample", all = FALSE)

#Calculating start & end times -----

# Extract start and end times from the 'time_range' column
start_end_times <- strsplit(taxa_count_data$time_range, "_", fixed = TRUE)
start_times <- sapply(start_end_times, function(x) x[1])

# Convert start and end times to POSIXct format
taxa_count_data$start_datetime <- as.POSIXct(paste(taxa_count_data$Date_collected, start_times),
                                             format = "%d/%m/%Y %H:%M", tz = "GMT")

# Calculate end datetime using 'time_hrs'
taxa_count_data$end_datetime <- taxa_count_data$start_datetime + 
  as.difftime(taxa_count_data$Time_hrs, units = "hours")

#Plot aesthetics 
time_colours <- c('#F28241','#19261C')

# Plot
uniq_species_plot <- ggplot(taxa_count_data, aes(x = reorder(time_range, start_datetime))) + 
  geom_bar(aes(y = Total_Species_Count, fill = "Total"), stat = "identity", width = 0.82) +
  geom_bar(aes(
    #y = combined_unique_count, #To use the unique for the 4hr /6hr
    y = Unique_All_Count, #To use the unique species across all samples
    fill = "Unique"), stat = "identity", width = 0.82) +
  labs(title = "Total and Unique Species per Sample",
       y = "Number of Species", x = "Collection Time (Hrs)",
       fill = "Species Count") +  # Set the fill legend title
  scale_fill_manual(values = c( "Total" = "#F28241", "Unique" = "#19261C"),  # Define colors for legend
                    labels = c("Total", "Unique" )) +  # Define labels for legend
  theme_minimal() +
  facet_wrap(~Time_hrs, scales = "free_x", ncol = 1,
             labeller = labeller(Time_hrs = c("4" = "4hr Samples", "6" = "6hr Samples"))) +
  theme(
    axis.title.x = element_blank(),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 12) )

# Save the plot as an .svg file with transparent background
ggsave("Graphs/2023_taxa_graphs/unique_across_all.svg", plot = uniq_species_plot, device = "svg", bg = "transparent")

DNA_yield <- ggplot(taxa_count_data, aes(x = reorder(time_range, start_datetime))) + 
  geom_bar(aes(y = Yield), stat = "identity", fill = "#19261C") +
  labs(title = "DNA Yield per Sample",
       y = "DNA (ng/ul)", x = "Collection Time (Hrs)") +
  theme_minimal() +
  facet_wrap(~Time_hrs, scales = "free_x", ncol = 1,
             labeller = labeller(Time_hrs = c("4" = "4hr Samples", "6" = "6hr Samples"))) +
  theme(
    axis.title.x = element_blank(),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 12) )

# Save the plot as an .svg file with transparent background
ggsave("Graphs/2023_taxa_graphs/DNA_yield.svg", plot = DNA_yield, device = "svg", bg = "transparent")







