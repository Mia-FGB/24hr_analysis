#!/usr/bin/env python3

# A script that takes a MARTi tsv file as an argument 
# Then filters it to select columns which are assigned / summed 
# outputs new tsv files with the filtered columns

# Usage:
# In 24hr_env conda environment
# python filter_marti_taxa.py <input_file> <output_prefix>

# Example usage:
# python Scripts/filter_marti_taxa.py 
# taxa_counts/2024_data/marti_assignments_noML_lca_0.1_all_levels_2025-FEB-10_10-36-18.tsv 
# metadata/2024_metadata/10_feb

import argparse 

def filter_marti_taxa(marti_file):
    with open(marti_file, "r") as f:
        lines = f.readlines()
        header = lines[0]
        header = header.strip().split("\t")
        
        #Pick out the different columns
        first_3_cols = header[:3]
        summed_cols = [col for col in header if "Summed" in col]
        identity_cols = [col for col in header if "identity" in col]
        assigned_cols = [col for col in header if col not in first_3_cols and col not in summed_cols and col not in identity_cols]
       

        #Create new data frames
        summed_df = [first_3_cols + summed_cols]
        assigned_df = [first_3_cols + assigned_cols]

        # Add the data rows
        for line in lines[1:]:
            values = line.strip().split("\t")
            summed_values = [values[header.index(col)] for col in summed_df[0]]
            assigned_values = [values[header.index(col)] for col in assigned_df[0]]
            summed_df.append(summed_values)
            assigned_df.append(assigned_values)

        # Modify the headers for summed_df and assigned_df
        summed_df[0] = first_3_cols + [col.split(" ")[0] for col in summed_cols]
        assigned_df[0] = first_3_cols + [col.split(" ")[0] for col in assigned_cols]

        return summed_df, assigned_df

def write_output(data, output_file):
    with open(output_file, "w") as f:
        for row in data:
            f.write("\t".join(row) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split MARTi taxa count file into separate assigned & summed reads.")
    parser.add_argument("marti_file", help="MARTi taxa count tsv file")
    parser.add_argument("output_prefix", help="Prefix to the output files")
    args = parser.parse_args()

    # Build the output file paths using the provided prefix
    summed_output_file = f"{args.output_prefix}_summed.tsv"
    assigned_output_file = f"{args.output_prefix}_assigned.tsv"

    summed_df, assigned_df = filter_marti_taxa(args.marti_file)

    write_output(summed_df, summed_output_file)
    write_output(assigned_df, assigned_output_file)