################################################################################
#
#	random_samples_hmp.R
#	Useful for pulling random samples of a site of interest from the HMP dataset
# 	Version 1.0
#	
#	Author: Jia Rong Wu
#	jwu424@gmail.com
#
#	random_samples_hmp.R is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	any later version.
#	
#	random_samples_hmp.R is distributed in the hope that it will be useful,
# 	but WITHOUT ANY WARRANTY; without even the implied warranty of
# 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	USAGE: 
#	DEP=$(R --slave --no-save '--args subsite subsamples' < random_samples_hmp.R)
#		subsite is the site you are extracting from
#		subsamples is the number of samples you're drawing (must be even)
#
#	Example:
#	DEP=$(R --slave --no-save '--args Tongue_dorsum 60' < random_samples_hmp.R)
################################################################################

library(GUniFrac)
args=(commandArgs(TRUE))

subsite <- as.character(args[1])
subsamples <- as.numeric(args[2])

# Read in files from HMP dataset
id <- read.table("data/v35_map_uniquebyPSN.txt", header=TRUE, sep="\t", row.names=1, comment.char="")
otu <- t( read.table("data/otu_table_psn_v35.txt", header=T, comment.char="", skip=1, sep="\t", row.names=1, check.names=F) )
rownames(otu) <- sub("^X", "", rownames(otu))

# Remove the taxonomic information
otu  <- as.data.frame(t(otu))
ind <- which(colnames(otu)=="Consensus Lineage")
otu[,ind] <- NULL
otu <- as.data.frame(t(otu))



# Extract reads from the main table of the subsite
reads.rn <- rownames(id)[which(id$HMPbodysubsite==subsite)]
reads <- otu[reads.rn,]
meta <- data.frame(id[reads.rn,],check.names=F)

# Subsample the reads randomly to the number provided
reads.sub <- reads[sample(1:nrow(reads), subsamples, replace=FALSE),]
reads.sub <- t(apply(reads.sub, 1, function(x){as.numeric(x)}))
meta.sub <- meta[rownames(reads.sub),]

colnames(reads.sub) <- colnames(reads)


# Set up false-positive vector 
loc <- c(rep("0",subsamples/2), rep("1",subsamples/2))
meta.sub$Location <- loc

# Strip the OTUs/features that are less than 100 across samples
reads.sub <- as.data.frame(t(reads.sub))
reads.sub.final <- reads.sub[which(apply(reads.sub, 1, sum) > 99),]
reads.sub.final <- as.data.frame(t(reads.sub.final))

# Extract the minimum count PER SAMPLE for rarefaction purposes
mn <- rowSums(reads.sub.final)
cat(min(mn))	# Echo the result of minimum sample depth for rarefaction 

dat <- "data/"
und <- "_"
ext <- ".txt"
met <- "meta"
red <- "reads"

readNames <- paste(dat,red,und,subsite,und,subsamples,ext, sep="")
metaNames <- paste(dat,met,und,subsite,und,subsamples,ext,sep="")

# Write original files
write.table( t(reads.sub.final), readNames, sep="\t", quote=F)
write.table(meta.sub, metaNames, sep="\t", quote=F)

# Implement a shell script to format the rep-set tree in order to generate for qiime
minimal <- as.data.frame(read.table("data/rep_set_v35_modified.fna", sep="\t", row.names=1, header=F))
minimal.sub <- as.data.frame(minimal[colnames(reads.sub.final),])
rownames(minimal.sub) <- colnames(reads.sub.final)
write.table(minimal.sub, "data/rep_set_v35_minimal_final_100.fna", sep="\t", quote=F, col.names=F)


# Create directories for the output of files
at <- "ait_all_"
for (i in 1:11)
{
	dirName <- paste(dat,at,i,sep="")
	dir.create(dirName)	
}

