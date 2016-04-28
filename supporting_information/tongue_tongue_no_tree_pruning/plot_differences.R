#!/usr/bin/Rscript

## FOR GG
# Plots a set of 30 tongue samples vs 30 tongue samples
# Should be no differences, but there's a difference in the sets of 
library(vegan)
source("../../UniFrac.r")

# read data
original.tongue.data <- read.table("../../data/tongue_dorsum/tongue_vs_tongue_30_forR.txt",sep="\t",check.names=FALSE,quote="",comment.char="", header=TRUE,row.names=1)
# tongue.tree <- read.tree("../data/tongue_dorsum/tongue_vs_tongue.tre")
tongue.tree <- read.tree("../../tongue_sample_migration/fasttree_all_seed_OTUs.tre")

# remove all OTUs with less than 100 counts across all samples
tongue.otu.sum <- apply(original.tongue.data,1,sum)
original.tongue.data <- original.tongue.data[which(tongue.otu.sum >= 100),]
tongue.otu.sum <- tongue.otu.sum[which(tongue.otu.sum>= 100)]

# make sure tree tip names match OTU names by taking out single quotes
tongue.tree$tip.label <- gsub("'","",tongue.tree$tip.label)

# remove OTUs that are not in the tree
original.tongue.data <- original.tongue.data[which(rownames(original.tongue.data) %in% tongue.tree$tip.label),]
original.tongue.data <- t(original.tongue.data)

# remove extra taxa from tree
absent <- tongue.tree$tip.label[!(tongue.tree$tip.label %in% colnames(original.tongue.data))]
if (length(absent) != 0) {
		tongue.tree <- drop.tip(tongue.tree, absent)
}

# root tree (rooted tree is required)
tongue.tree <- midpoint(tongue.tree)

# tongue and cheek data have more read counts per sample, so we're rarefying to the lowest number of per sample counts in tongue only data
d.tongue.data <- read.table("../../data/tongue_dorsum/A_otu_table_tab_r659_original.txt",sep="\t",quote="",header=TRUE,row.names=1,check.names=FALSE)
e.tongue.data <- read.table("../../data/tongue_dorsum/B_otu_table_tab_r659_original.txt",sep="\t",quote="",header=TRUE,row.names=1,check.names=FALSE)

d.tongue.data <- t(d.tongue.data)
e.tongue.data <- t(e.tongue.data)

# remove OTUs with all zero counts after rarefaction from tree
d.tongue.otu.sum <- apply(d.tongue.data,2,sum)
d.tongue.data <- d.tongue.data[,which(d.tongue.otu.sum > 0)]
d.tongue.tree <- tongue.tree
absent <- d.tongue.tree$tip.label[!(d.tongue.tree$tip.label %in% colnames(d.tongue.data))]
if (length(absent) != 0) {
		d.tongue.tree <- drop.tip(d.tongue.tree, absent)
}

e.tongue.otu.sum <- apply(e.tongue.data,2,sum)
e.tongue.data <- e.tongue.data[,which(e.tongue.otu.sum > 0)]
e.tongue.tree <- tongue.tree
absent <- e.tongue.tree$tip.label[!(e.tongue.tree$tip.label %in% colnames(e.tongue.data))]
if (length(absent) != 0) {
		e.tongue.tree <- drop.tip(e.tongue.tree, absent)
}

d.tongue.unifrac <- getDistanceMatrix(d.tongue.data,d.tongue.tree,method="unweighted",verbose=TRUE,pruneTree=FALSE)
e.tongue.unifrac <- getDistanceMatrix(e.tongue.data,e.tongue.tree,method="unweighted",verbose=TRUE,pruneTree=FALSE)

d.tongue <- pcoa(d.tongue.unifrac)
e.tongue <- pcoa(e.tongue.unifrac)

#function to get variance explained for the PCOA component labels
getVarExplained <- function(vector) {
	rawVarEx <- apply(vector,2,function(x) sd(x)*sd(x))
	totalVarExplained <- sum(rawVarEx)
	varEx <- rawVarEx/totalVarExplained
	return(varEx)
}

plotMigration <- function(d,e) {
	d.varEx <- getVarExplained(d$vectors)
	# convert to percentage
	d.varEx <- d.varEx * 100

	# Setup axis labelling conventions
	per <- "%"
	x.pc1.explained <- d.varEx[1]
	x.pc1.explained <- round(x.pc1.explained, digits=1)
	xlabel <- "Principal Component1 Eigenvalues: "
	xla <- paste(xlabel, x.pc1.explained,per, sep="", collapse=NULL)

	x.pc2.explained <- d.varEx[2]
	x.pc2.explained <- round(x.pc2.explained, digits=1)
	ylabel <- "Principal Component2 Eigenvalues: "
	yla <- paste(ylabel, x.pc2.explained,per, sep="", collapse=NULL)

	#perform procrustes fit
	fit <- procrustes(d$vectors, e$vectors)

	rtl <- NULL
	ltr <- NULL

	component <- 1

	for (i in 1:60)	# Iterate through all 60 samples
	{
		if ( (d$vectors[i,component] > 0) & (fit$Yrot[i,component] < 0))
		{	
			print("Right Cluster to Left Cluster Movement: ")	# Moved from right to left
			print(rownames(d$vectors)[i])
			rtl <- c(rtl, rownames(d$vectors)[i])	# Get list of Right to Left Movement
		}	
	}
	for (i in 1:60)
	{
		if ( (d$vectors[i,component] < 0) & (fit$Yrot[i,component] > 0))
		{
			print("Left Cluster to Right Cluster Movement: ")	# Moved from left to right
			print(rownames(d$vectors)[i])
			ltr <- c(ltr, rownames(d$vectors)[i])	# Get list of Left to Right movement
		}
	}

	# First/second refer to the columns the data is being read from 
	first <- 1
	second <- 2

	plot(d$vectors[1:60,first], d$vectors[1:60,second],main="Sample Migration between Rarefactions",xlab=xla,ylab=yla, pch=19, col=rgb(0,0,0,0.4))	# Plot the 1st rarefation

	shapes <- c(19,15,19,15,19)
	colours <- c(rgb(1,0,0,0.4), rgb(1,0,0,0.4), rgb(0,0,1,0.4), rgb(0,0,1,0.4), rgb(0,0,0,0.4) )

	#legend(x=-0.25,y=0.31,title="Sample Movement", legend = c("Left->Right","Origin","Right->Left","Origin", "No Change"), pch=shapes, col=colours)

	# Points that have moved from right to left are red
	# Points that have moved from left to right are blue

	# Points represent a plot of Rarefaction 1, but are coloured based on sample movement of Rare 2
	points(d$vectors[rtl,first][fit$Yrot[rtl,first] < 0], d$vectors[rtl,second][fit$Yrot[rtl,first] < 0], pch=19, col=rgb(1,0,0,0.4))	# Red
	points(d$vectors[ltr,first][fit$Yrot[ltr,first] > 0], d$vectors[ltr,second][fit$Yrot[ltr,first] > 0], pch=19, col=rgb(0,0,1,0.4))	# Blue

	# The squares indicate WHERE they have moved from 
	# Can comment this out to remove square drawing
	points(fit$Yrot[rtl,first], fit$Yrot[rtl,second], pch = 15, col=rgb(1,0,0,0.4))
	points(fit$Yrot[ltr,first], fit$Yrot[ltr,second], pch = 15, col=rgb(0,0,1,0.4))

	# Draw line segments between the movements
	arrows(fit$Yrot[rtl,first], fit$Yrot[rtl, second], d$vectors[rtl, first], d$vectors[rtl, second], col = rgb(1,0,0,0.5), length=0.1)

	arrows(fit$Yrot[ltr,first], fit$Yrot[ltr, second], d$vectors[ltr, first], d$vectors[ltr, second], col = rgb(0,0,1,0.5), length=0.1)
}



pdf("UniFrac_tvst_movement.pdf")	# Comment out if not plotting

plotMigration(d.tongue,e.tongue)

dev.off()







