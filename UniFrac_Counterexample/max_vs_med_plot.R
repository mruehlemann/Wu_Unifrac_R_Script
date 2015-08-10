################################################################################
#
#	max_vs_med_plot.R
#	Useful for visualizing differences between different analysis
# 	Version 1.0
#	
#	Author: Jia Rong Wu
#	jwu424@gmail.com
#
#	max_vs_med_plot.R is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	any later version.
#	
#	max_vs_med_plot.R is distributed in the hope that it will be useful,
# 	but WITHOUT ANY WARRANTY; without even the implied warranty of
# 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
################################################################################

library(ALDEx2)
library(vegan)
library(ape)
args=(commandArgs(TRUE))
raref <- as.numeric(args[1])
subsite <- as.character(args[2])

component <- 1			# This is the component being examined (i.e. PC1)
minRange <- 0			# Lower limit of values across each rarefaction
maxRange <- 0			# Upper limit of value  across each rarefaction
xMax <- 0				# Maximum x-axis of the plot
yMax <- 0				# Maximum y-axis of the plot

numFiles <- 10			# Edit this number based on the number of rarefactions
numSamples <- 60		# Edit this number based on number of samples you have
numAnalysis <- 4		# Specifies the number of DIFFERENT analysis
folders <- c("data/uni_all_", "data/uni_all_", "data/brc_all_", "data/ait_all_")
matricies <- c("/unweighted_unifrac_pc", "/weighted_unifrac_pc", paste("/bdiv_even",raref,"/bray_curtis_pc",sep=""), "/atc_pc")
# Create list of values in order to hold the plot information
listPlot <- as.list (c(1:numAnalysis))
for (i in 1:numAnalysis)
{
	listPlot[[i]] <- data.frame(matrix("", ncol=2, nrow=numFiles), stringsAsFactors=FALSE)
}

################################################################################
#####								Analysis								####
################################################################################

for (j in 1: numAnalysis)
{
	# Fetch the analysis data being plotted 
	folder <- folders[j]
	matrix <- matricies[j]
	extension <- ".txt"
	
	# Iterate through set of files to determine the range of all the PC1's
	for (i in 1:numFiles)
	{
		# Dynamically read in set of files
		fileName <- paste(folder,i,matrix,extension, sep="")	
		analysis <- read.table(fileName, sep="\t", header=T, row.names=1)
		
		# Determine the extremes of the set of components (PC1)
		# Using the extremes prevents any values from being > 1
		tempMin <- min(analysis[1:numSamples,component])
		minRange <- if(tempMin < minRange) tempMin else minRange
		tempMax <- max(analysis[1:numSamples,component])
		maxRange <- if(tempMax > maxRange) tempMax else maxRange
	}	
	
	weighted.range <- abs(minRange - maxRange)
	
	# Fetch the reference matrix to plot against for each analysis
	# By default, the first rarefaction/dirichlet instance  is the default reference
	fileName = paste(folder,1,matrix,extension,sep="")
	ref.plot <- read.table(fileName, sep="\t", header=T, row.names=1)
	ref.plot <- ref.plot[1:numSamples,component]
	
	# Fetch every other value to plot
	for (i in 2:numFiles)
	{
		fileName = paste(folder,i,matrix,extension,sep="")
		new.plot <- read.table(fileName, sep="\t", header=T, row.names=1)
		new.plot <- new.plot[1:numSamples,component]
		
		new.plot <- as.vector(procrustes(ref.plot, new.plot)$Yrot)
		med.plot <- ( abs(ref.plot - new.plot)/weighted.range  ) 
		
		listPlot[[j]][i,1] <- median(med.plot)
		listPlot[[j]][i,2] <- max(med.plot)
		
		xMax <- if(xMax < median(med.plot)) median(med.plot)+.05 else xMax
		yMax <- if(yMax < max(med.plot)) max(med.plot)+0.05 else yMax
	}
}

################################################################################
#####								Plotting								####
################################################################################

tla <- "Max Relative Deviation of Rarefactions vs Median Deviation"
xla <- "Median Relative Deviation (Vres)"
yla <- "Max Relative Deviation (Vres)"

# Increment here for insertion of values
shapes <- c(20,20,20,20)
colours <- c(rgb(1,0,0,0.5), rgb(0,1,0,0.5), rgb(0,0,1,0.5), rgb(0,0,0,0.5) )

#dir.create("Figures")

pdf(paste("Figures/MRD_vs_MD_",subsite,".pdf",sep=""))

plot(0, xlim=c(0,xMax), ylim=c(0,yMax), xlab=xla, ylab=yla, main=tla)
legend(x=xMax-0.03,y=0.1, legend=c("un_uniFrac", "w_uniFrac", "bray_curtis", "aitchison"), pch=shapes,
col=colours)

# Plot the points
for (i in 1:numAnalysis)
{
	points(listPlot[[i]][,1], listPlot[[i]][,2], pch=shapes[i], col=colours[i], cex=2)
}
dev.off()
