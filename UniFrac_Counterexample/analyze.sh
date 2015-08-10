#!/bin/bash -l
# Make sure you launch macqiime before invoking this script*
# Dependencies:
#	1) MacQIIME 1.8.0-20140103
#	2) R 3.2.1	w/ packages and their respective dependencies: 
#		i) GUniFrac
#		ii) ALDEx2
#		iii) vegan
#		iv) ape

# Download the OTU table, Metadata table and Tree from these links respectively:
#	http://downloads.hmpdacc.org/data/HMQCP/otu_table_psn_v35.txt.gz
#	http://downloads.hmpdacc.org/data/HMQCP/v35_map_uniquebyPSN.txt.bz2
#	http://downloads.hmpdacc.org/data/HMQCP/rep_set_v35.fna.gz
# Put the 3 files in a folder called data/
# Analysis and scripts should be in the parent directory of data/ 
#	i.e
#	FooBar/	(scripts go here)
#	FooBar/data/ (data goes here)
#	
# Simply invoking ./analyze.sh will run the files specified
gunzip rep_set_v35.fna.gz	# Unzips to rep_set_v35.fna 

# Format tree file for use in R to select files
sed -e 's/>//g' data/rep_set_v35.fna > data/rep_set_v35_modified.fna
tr ' ' '\t' < data/rep_set_v35_modified.fna > data/rep_set_v35_modified2.fna
sed 'N;s/\n/+/g' data/rep_set_v35_modified2.fna > data/rep_set_v35_modified.fna
rm data/rep_set_v35_modified2.fna
# data/rep_set_v35_modified.fna is the file that may be modified


# Subsample from the original set of files
# DEP represents the minimum sample depth for rarefaction purposes
# Subset the number of features to use in generating a tree
subsite="Stool"
subsamples=60
DEP=$(R --slave --no-save --args $subsite $subsamples < random_samples_hmp.R)
echo $DEP


# Convert back to the original format for conversion into tree format
sed -e 's/O/>O/g' data/rep_set_v35_minimal_final_100.fna > data/rep_set_v35_modified2.fna
tr '+' '\n' < data/rep_set_v35_modified2.fna > data/rep_set_v35_modified.fna
tr '\t' ' ' < data/rep_set_v35_modified.fna > data/rep_set_v35_modified2.fna
mv data/rep_set_v35_modified2.fna data/rep_set_v35_modified.fna


# Call muscle to create the tree
mkdir muscle
muscle -in data/rep_set_v35_modified.fna -out muscle/all_seed_OTUs_bad.mfa
awk '/^>/{gsub(/lcl\|[0-9]+\|num\|[0-9]+\|OTU\|/, "")}; 1' muscle/all_seed_OTUs_bad.mfa > muscle/all_seed_OTUs.mfa
rm muscle/all_seed_OTUs_bad.mfa
make_phylogeny.py -i muscle/all_seed_OTUs.mfa -o muscle/fasttree_all_seed_OTUs.tre

# Call generate CLR PCoA's to populate ait_all_i folders with PC1 values
R CMD BATCH generate_aitchison_data.R


# meta_sample_subsamples.txt needs #SampleID prepended to it
# reads_sample_subsamples.txt needs #OTU_ID prepended to it
dir="data/"
red="reads"
met="meta"
und="_"
ext=".txt"

reads=$(cat $dir$red$und$subsite$und$subsamples$ext)
echo -en "#SampleID\t$reads" > $dir$red$und$subsite$und$subsamples$ext

meta=$(cat $dir$met$und$subsite$und$subsamples$ext)
echo -en "#OTU_ID\t$meta" > $dir$met$und$subsite$und$subsamples$ext



# single_rarefaction.py
# Generate the biom file from the set of data 
ot="otu_table_"
org="original"
bm=".biom"

biom convert -i $dir$red$und$subsite$und$subsamples$ext -o $dir$ot$org$bm --table-type="otu table"


# Create separate rarefactions for each set in the data
# Rarefy the data 10 times 
# Generate the tables
rar="rarefied"

COUNTER=0
while [ $COUNTER -lt 11 ]; do
        echo $COUNTER
        let COUNTER=COUNTER+1
        single_rarefaction.py -i $dir$ot$org$bm -o $dir$ot$rar$COUNTER$bm -d $DEP
done


bdiv="uni_all"
brc="brc_all"
# Beta/core analysis
COUNTER=1
while [ $COUNTER -lt 11 ]; do
	echo $COUNTER	# Display progress
	beta_diversity_through_plots.py -i $dir$ot$rar$COUNTER$bm -m $dir$met$und$subsite$und$subsamples$ext -o $dir$bdiv$und$COUNTER -t muscle/fasttree_all_seed_OTUs.tre -f
	core_diversity_analyses.py -i $dir$ot$rar$COUNTER$bm -m $dir$met$und$subsite$und$subsamples$ext -o $dir$brc$und$COUNTER --nonphylogenetic_diversity -e $DEP
    let COUNTER=COUNTER+1
done


# Call respective diversity analysis on the set of data
# Use aitchison to generate the data from aitchison 

# Plot the data w/ max vs med plot and save the files in a folder called outputs
# Clean the directories and all the separate files 

R --slave --no-save --args $DEP < max_vs_med_plot.R





