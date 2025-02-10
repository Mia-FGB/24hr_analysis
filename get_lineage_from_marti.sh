#!/bin/bash
#https://bioinf.shenwei.me/taxonkit/usage/
 #I may have a better script than this which I can use now - will need to look

#Name of sample, will change this with each analysis
sample=species_4h_1_2_5

#Prep the taxaIDs_counts file
cut -f 2 -d $',' ${sample}.csv > ${sample}_taxaID.txt

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


