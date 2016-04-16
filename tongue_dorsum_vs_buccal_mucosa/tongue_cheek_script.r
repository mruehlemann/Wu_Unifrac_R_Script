#!/usr/bin/env Rscript 
options(error=recover)

#this script prints out PDF pcoa plots and distance matrices, given an OTU table, phylogenetic tree, and metadata
# run: nohup Rscript ./tongue_cheek_script.r > tongue_cheek_script_nohup.out 2>&1&

source("../UniFrac.r")
library(ape)
library(phangorn)
library(vegan)
library(stringr)

otu.tab <- read.table("../data/tongue_dorsum_vs_buccal_mucosa/hmp_tongue_cheek_data.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

#remove taxonomy column to make otu count matrix numeric
taxonomy <- rownames(otu.tab)
#make row names samples, and col names OTUs
otu.tab <- t(as.matrix(otu.tab))

# filter out taxa that are less than 1% abundant in all samples
sample.sum <- apply(otu.tab,1,sum)
one.percent <- sample.sum*0.01
otus.greater <- apply(otu.tab,2,function(x) length(which(x > one.percent)))
otu.tab <- otu.tab[,which(otus.greater > 0)]

tree <- read.tree("../data/tongue_dorsum_vs_buccal_mucosa/hmp_tongue_cheek_subtree.tre")

# remove extra taxa from tree
absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
if (length(absent) != 0) {
		tree <- drop.tip(tree, absent)
}

# root tree (rooted tree is required)
tree <- midpoint(tree)

#sort taxa from most to least abundant
taxaOrder <- rev(order(apply(otu.tab,2,sum)))
taxonomy <- taxonomy[taxaOrder]
otu.tab <- otu.tab[,taxaOrder]

# read metadata
MyMetaOrdered <- rownames(otu.tab)
MyMetaOrdered <- gsub("_.*$","",MyMetaOrdered)

#rarefy data for unweighted unifrac
otu.tab.rarefy <- rrarefy(otu.tab, min(apply(otu.tab,1,sum)))

# do Unweighted UniFrac separately with rarefied data
unweighted <- getDistanceMatrix(otu.tab.rarefy,tree,method="unweighted",verbose=TRUE)
write.table(unweighted,file="output/unweighted_distance_matrix.txt",sep="\t",quote=FALSE)

#calculate distance matrix
all_distance_matrices <- getDistanceMatrix(otu.tab,tree,method="all",verbose=TRUE)

#output distance matrices
weighted <- all_distance_matrices[["weighted"]]
write.table(weighted,file="output/weighted_distance_matrix.txt",sep="\t",quote=FALSE)
information <- all_distance_matrices[["information"]]
write.table(information,file="output/information_distance_matrix.txt",sep="\t",quote=FALSE)
ratio_no_log <- all_distance_matrices[["ratio_no_log"]]
write.table(ratio_no_log,file="output/ratio_distance_matrix.txt",sep="\t",quote=FALSE)
# create bray curtis dist object using vegan, and turn into distance matrix
braycurtis.vegdist <- vegdist(otu.tab,method="bray")
braycurtis <- matrix(nrow=nrow(otu.tab),ncol=nrow(otu.tab))
braycurtis[lower.tri(braycurtis)] <- braycurtis.vegdist
diag(braycurtis) <- 0
braycurtis.vegdist <- vegdist(otu.tab,method="bray",upper=TRUE)
braycurtis[upper.tri(braycurtis)] <- braycurtis.vegdist
write.table(braycurtis,file="output/bray_curtis_distance_matrix.txt",sep="\t",quote=FALSE)

# unweighted <- read.table("tongue_cheek_output/unweighted_distance_matrix.txt", sep = "\t", header = TRUE, row.names=1, quote = "")
# weighted <- read.table("tongue_cheek_output/weighted_distance_matrix.txt", sep = "\t", header = TRUE, row.names=1, quote = "")
# information <- read.table("tongue_cheek_output/information_distance_matrix.txt", sep = "\t", header = TRUE, row.names=1, quote = "")
# ratio_no_log <- read.table("tongue_cheek_output/ratio_no_log_normalize_distance_matrix.txt", sep = "\t", header = TRUE, row.names=1, quote = "")

unweighted.pcoa <- pcoa(unweighted)
weighted.pcoa <- pcoa(weighted)
information.pcoa <- pcoa(information)
ratio_no_log.pcoa <- pcoa(ratio_no_log)
braycurtis.pcoa <- list()
braycurtis.pcoa$vectors <- cmdscale(braycurtis,k=nrow(otu.tab)-1,add=TRUE)$points
rownames(braycurtis.pcoa$vectors) <- rownames(otu.tab)
colnames(braycurtis.pcoa$vectors) <- paste("Axis",c(1:ncol(braycurtis.pcoa$vectors)),sep='.')


#function to get variance explained for the PCOA component labels
getVarExplained <- function(vector) {
	rawVarEx <- apply(vector,2,function(x) sd(x)*sd(x))
	totalVarExplained <- sum(rawVarEx)
	varEx <- rawVarEx/totalVarExplained
	return(varEx)
}


unweighted.varEx <- getVarExplained(unweighted.pcoa$vectors)
weighted.varEx <- getVarExplained(weighted.pcoa$vectors)
information.varEx <- getVarExplained(information.pcoa$vectors)
ratio_no_log.varEx <- getVarExplained(ratio_no_log.pcoa$vectors)
braycurtis.varEx <- getVarExplained(braycurtis.pcoa$vectors)

#get vector version of distance matrices for correlation plots below
unweighted.vector <- unlist(unweighted[lower.tri(unweighted,diag=TRUE)])
weighted.vector <- unlist(weighted[lower.tri(weighted,diag=TRUE)])
information.vector <- unlist(information[lower.tri(information,diag=TRUE)])
ratio_no_log.vector <- unlist(ratio_no_log[lower.tri(ratio_no_log,diag=TRUE)])
braycurtis.vector <- unlist(braycurtis[lower.tri(braycurtis,diag=TRUE)])

# replace abbreviations with full body site names (there aren't actually any dominant taxa in this data set)
taxonomyGroups <- as.factor(c("Buccal Mucosa", "Tongue Dorsum"))
groups <- str_extract(rownames(otu.tab),"^[a-z]*")
groups[which(groups=="bm")] <- "Buccal Mucosa"
groups[which(groups=="td")] <- "Tongue Dorsum"
groups <- as.factor(groups)

palette(c("black", "blue", "red"))
dev.off()

pdf("output/pcoa_plots.pdf")
par(oma=c(1,1,1,7))
#plot pcoa plots
plot(unweighted.pcoa$vectors[,1],unweighted.pcoa$vectors[,2], col=groups,main="Unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Coordinate", round(unweighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Coordinate", round(unweighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(taxonomyGroups), pch=rep(19,length(taxonomyGroups)), col=palette()[1:length(taxonomyGroups)], xpd=NA, inset=c(-0.4,0))
plot(weighted.pcoa$vectors[,1],weighted.pcoa$vectors[,2], col=groups,main="Weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Coordinate", round(weighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Coordinate", round(weighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(taxonomyGroups), pch=rep(19,length(taxonomyGroups)), col=palette()[1:length(taxonomyGroups)], xpd=NA, inset=c(-0.4,0))
plot(information.pcoa$vectors[,1],information.pcoa$vectors[,2], col=groups,main="Information UniFrac\nprincipal coordinates analysis",xlab=paste("First Coordinate", round(information.varEx[1],digits=3),"variance explained"),ylab=paste("Second Coordinate", round(information.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(taxonomyGroups), pch=rep(19,length(taxonomyGroups)), col=palette()[1:length(taxonomyGroups)], xpd=NA, inset=c(-0.4,0))
plot(ratio_no_log.pcoa$vectors[,1],ratio_no_log.pcoa$vectors[,2], col=groups,main="Ratio UniFrac\nprincipal coordinates analysis",xlab=paste("First Coordinate", round(ratio_no_log.varEx[1],digits=3),"variance explained"),ylab=paste("Second Coordinate", round(ratio_no_log.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(taxonomyGroups), pch=rep(19,length(taxonomyGroups)), col=palette()[1:length(taxonomyGroups)], xpd=NA, inset=c(-0.4,0))
plot(braycurtis.pcoa$vectors[,1],braycurtis.pcoa$vectors[,2], col=groups,main="Bray Curtis Dissimilarity\nprincipal coordinates analysis",xlab=paste("First Coordinate", round(braycurtis.varEx[1],digits=3),"variance explained"),ylab=paste("Second Coordinate", round(braycurtis.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(taxonomyGroups), pch=rep(19,length(taxonomyGroups)), col=palette()[1:length(taxonomyGroups)], xpd=NA, inset=c(-0.4,0))

#plot correlation between different UniFrac modes
plot(unweighted.vector,information.vector,main="Unweighted vs. information UniFrac")
plot(weighted.vector,information.vector,main="Weighted vs. information UniFrac")
plot(unweighted.vector,weighted.vector,main="Unweighted vs. weighted UniFrac")
plot(ratio_no_log.vector,weighted.vector,main="Ratio vs. weighted UniFrac")
plot(braycurtis.vector,weighted.vector,main="Bray Curtis dissimilarity vs. weighted UniFrac")

dev.off()
