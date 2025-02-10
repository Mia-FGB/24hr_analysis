# 24-Hour Air Sample Analysis

This repository contains code for analyzing air samples collected over a 24-hour period in 2023 and 2024. The analysis involves processing, visualization, and statistical exploration of microbial taxa found in these air samples.

## Contents

### Notebooks and Scripts
- **2023 Analysis:**
  - `12_24_hr_2023_analysis.ipynb` – Jupyter notebook for analyzing 12-hour and 24-hour air samples from 2023.
  - `marti_4_6_hr_2023_analysis.ipynb` – Jupyter notebook for analyzing 4-hour and 6-hour samples from 2023 data.
  - `Phyloseq_24hrCub_2023.R` – R script for analyzing 2023 air samples.
  
- **2024 Analysis:**
  - `24hr_analysis.ipynb` – Jupyter notebook for analyzing the 2024 air samples.
  - `24hr_2024_stacked_bar.R` – R script to generate stacked bar charts of microbial taxa.

- **Supporting Scripts:**
  - `get_lineage_from_marti.sh` – Shell script to extract lineage data from MARTi input files.
  - `identity_from_json_working.ipynb` – Jupyter notebook for extracting identity information from MARTi JSON files and analyzing species-specific changes.
  - `split_marti_taxa.py` – Python script to process MARTi taxa count output by splitting it into assigned and summed data frames.
  - `Total_Uniq_Species.R` – R script to calculate and visualize the proportion of shared and unique species between groups using stacked bar charts.


