################################################################################
#
#	generate_aitchison_data.R
#	Ancillary function to max_vs_med_plot.R 
#	Generates the Aitchison CLR-transformed data and the respective dirichlet 
#	instances for a set of OTU's
# 	Version 1.0
#	
#	Author: Jia Rong Wu
#	jwu424@gmail.com
#
#	generate_aitchison_data.R is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	any later version.
#	
#	generate_aitchison_data.R is distributed in the hope that it will be useful,
# 	but WITHOUT ANY WARRANTY; without even the implied warranty of
# 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
################################################################################

library(ALDEx2)
library(GUniFrac)
args=(commandArgs(TRUE))

subsite <- as.character(args[1])
subsamples <- as.numeric(args[2])

filename <- paste("data/reads_",subsite,"_",subsamples,".txt",sep="")

otus <- t(read.table(filename, sep="\t", header=T, row.names=1))
rownames(otus) <- gsub("^X","", rownames(otus))	# Clean up formatting
otus <- as.data.frame(otus)
clr.otus <- aldex.clr(otus, mc.samples=128)		# Generate CLR Data

# Set up the folder names
folder <- "data/ait_all_"
name <- "/atc_pc"
ext <- ".txt"

for (i in 1:11)
{
	filename <- paste(folder,i,name,ext,sep="")	
	clr.otu.inst <- sapply(getMonteCarloInstances(clr.otus), function(y){y[,i]})	
	ait.dist <- dist(clr.otu.inst)
	ait.matrix <- as.matrix(ait.dist)
	ait.p <- pcoa(ait.matrix)
	ait.pcoa <- ait.p[4]
	ait.pcoa <- as.data.frame(ait.pcoa)
	a.1 <- ait.pcoa$vectors.Axis.1
	
	write.table(a.1, file=filename, sep="\t", quote=F)
}
