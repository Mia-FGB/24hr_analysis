#!/bin/bash
#https://bioinf.shenwei.me/taxonkit/usage/
# Run in taxonkit_environment
# Run from  taxa_counts/2024_data

#MARTi file name
sample=marti_assignments_noML_lca_0.1_all_levels_2025-FEB-10_10-36-18

# Extract just the 2nd col which is the NCBI IDs
cut -f 2 ${sample}.tsv > ${sample}_taxaID.txt

#Get lineages
taxonkit lineage ${sample}_taxaID.txt > ${sample}_taxaID_lineage.txt

#Fill in blanks
taxonkit reformat ${sample}_taxaID_lineage.txt -r Unassigned | cut -f 1,3 > ${sample}_taxaID_lineage_clean.txt

#Change delimiter
sed 's/;/\t/g' ${sample}_taxaID_lineage_clean.txt | awk -F'\t' 'BEGIN {OFS=","} { print $1, $2, $3, $4, $5, $6, $7, $8 }' > ${sample}_taxaID_lineage_sep.csv

#Add headers
echo "taxid,kingdom,phylum,class,order,family,genus,species" > header.txt
cat header.txt ${sample}_taxaID_lineage_sep.csv > ${sample}_taxaID_lineage.csv

#Remove all the intermediate files generated
rm ${sample}_taxaID.txt ${sample}_taxaID_lineage.txt ${sample}_taxaID_lineage_clean.txt ${sample}_taxaID_lineage_sep.csv header.txt
