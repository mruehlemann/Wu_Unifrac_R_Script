#!/usr/bin/Rscript
library(ape)
library(phangorn)
library(vegan)
library(car)
source("../UniFrac.r")

# read data
d <- read.table("../data/tongue_dorsum/tongue_vs_tongue_30_forR.txt",sep="\t",check.names=FALSE,quote="",comment.char="", header=TRUE,row.names=1)
tree <- read.tree("../data/tongue_dorsum/fasttree_all_seed_OTUs.tre")

# remove all OTUs with zero counts across all samples (this data set should have a min. of 100 counts per OTU, but this code is here just in case.)
otu.sum <- apply(d,1,sum)
d <- d[which(otu.sum > 0),]
otu.sum <- otu.sum[which(otu.sum>0)]

# make sure tree tip names match OTU names by taking out single quotes
tree$tip.label <- gsub("'","",tree$tip.label)
# remove OTUs that are not in the tree
d <- d[which(rownames(d) %in% tree$tip.label),]

# remove extra OTUs from tree
absent <- tree$tip.label[!(tree$tip.label %in% rownames(d))]
if (length(absent) != 0) {
		tree <- drop.tip(tree, absent)
}

# root tree (required for UniFrac)
tree <- midpoint(tree)

# generate 100 rarefactions
d <- t(d)
rarefactions <- list()
for (i in c(1:100)) {
  rarefactions[[i]] <- rrarefy(d, min(otu.sum))
}

# calculate unifrac distance matrices
unifrac <- list()
for (i in c(1:100)) {
  unifrac[[i]] <- list()
  rare.otu.sum <- apply(rarefactions[[i]],2,sum)
  rarefactions[[i]] <- rarefactions[[i]][,which(rare.otu.sum > 0)]
  empty.tips <- tree$tip.label[which(!(tree$tip.label %in% colnames(rarefactions[[i]])))]
  temp.tree <- tree
  temp.tree <- drop.tip(temp.tree, empty.tips)
  all_distance_matrices <- getDistanceMatrix(rarefactions[[i]],temp.tree,method="all",verbose=TRUE)
  unifrac[[i]][["unweighted"]] <- all_distance_matrices[["unweighted"]]
  unifrac[[i]][["weighted"]] <- all_distance_matrices[["weighted"]]
  unifrac[[i]][["information"]] <- all_distance_matrices[["information"]]
  unifrac[[i]][["ratio"]] <- all_distance_matrices[["ratio"]]
  # get bray curtis distance through vegan, convert into distance matrix format
  braycurtis.vegdist <- vegdist(rarefactions[[i]],method="bray")
  unifrac[[i]][["braycurtis"]] <- matrix(nrow=nrow(rarefactions[[i]]),ncol=nrow(rarefactions[[i]]))
  unifrac[[i]][["braycurtis"]][lower.tri(unifrac[[i]][["braycurtis"]])] <- braycurtis.vegdist
  diag(unifrac[[i]][["braycurtis"]]) <- 0
  braycurtis.vegdist <- vegdist(rarefactions[[i]],method="bray",upper=TRUE)
  unifrac[[i]][["braycurtis"]][upper.tri(unifrac[[i]][["braycurtis"]])] <- braycurtis.vegdist
}

# make pcoa
axis.1 <- list()
for (i in c(1:100)) {
  axis.1[[i]] <- list()
  axis.1[[i]][["unweighted"]] <- pcoa(unifrac[[i]][["unweighted"]])$vectors[,1]
  axis.1[[i]][["weighted"]] <- pcoa(unifrac[[i]][["weighted"]])$vectors[,1]
  axis.1[[i]][["information"]] <- pcoa(unifrac[[i]][["information"]])$vectors[,1]
  axis.1[[i]][["ratio"]] <- pcoa(unifrac[[i]][["ratio"]])$vectors[,1]
  # fake the pcoa for bray curtis using cmdscale (allows for non euclidean input)
  axis.1[[i]][["braycurtis"]] <- cmdscale(unifrac[[i]][["braycurtis"]],k=nrow(rarefactions[[i]])-1,add=TRUE)$points[,1]
}

# perform procrustes fit against first rarefaction
target <- 1
fit.plots <- list()
fit.plots[[target]] <- list()
fit.plots[[target]][["unweighted"]] <- axis.1[[target]][["unweighted"]]
fit.plots[[target]][["weighted"]] <- axis.1[[target]][["weighted"]]
fit.plots[[target]][["information"]] <- axis.1[[target]][["information"]]
fit.plots[[target]][["ratio"]] <- axis.1[[target]][["ratio"]]
fit.plots[[target]][["braycurtis"]] <- axis.1[[target]][["braycurtis"]]

for (i in c(2:100)) {
  fit.plots[[i]] <- list()
  fit.plots[[i]][["unweighted"]] <- as.vector(procrustes(fit.plots[[target]][["unweighted"]], axis.1[[i]][["unweighted"]]))
  fit.plots[[i]][["weighted"]] <- as.vector(procrustes(fit.plots[[target]][["weighted"]], axis.1[[i]][["weighted"]]))
  fit.plots[[i]][["information"]] <- as.vector(procrustes(fit.plots[[target]][["information"]], axis.1[[i]][["information"]]))
  fit.plots[[i]][["ratio"]] <- as.vector(procrustes(fit.plots[[target]][["ratio"]], axis.1[[i]][["ratio"]]))
  fit.plots[[i]][["braycurtis"]] <- as.vector(procrustes(fit.plots[[target]][["braycurtis"]], axis.1[[i]][["braycurtis"]]))
}

# get mean/max deviation
mean.dev <- list()
max.dev <- list()
minRange <- list()
maxRange <- list()
minRange[["unweighted"]] <- min(fit.plots[[target]][["unweighted"]])
minRange[["weighted"]] <- min(fit.plots[[target]][["weighted"]])
minRange[["information"]] <- min(fit.plots[[target]][["information"]])
minRange[["ratio"]] <- min(fit.plots[[target]][["ratio"]])
minRange[["braycurtis"]] <- min(fit.plots[[target]][["braycurtis"]])

maxRange[["unweighted"]] <- max(fit.plots[[target]][["unweighted"]])
maxRange[["weighted"]] <- max(fit.plots[[target]][["weighted"]])
maxRange[["information"]] <- max(fit.plots[[target]][["information"]])
maxRange[["ratio"]] <- max(fit.plots[[target]][["ratio"]])
maxRange[["braycurtis"]] <- max(fit.plots[[target]][["braycurtis"]])

for (i in c(2:100)) {
  mean.dev[[i]] <- list()
  mean.dev[[i]][["unweighted"]] <- mean(abs(fit.plots[[i]]$unweighted$Yrot - fit.plots[[target]][["unweighted"]]))
  mean.dev[[i]][["weighted"]] <- mean(abs(fit.plots[[i]]$weighted$Yrot - fit.plots[[target]][["weighted"]]))
  mean.dev[[i]][["information"]] <- mean(abs(fit.plots[[i]]$information$Yrot - fit.plots[[target]][["information"]]))
  mean.dev[[i]][["ratio"]] <- mean(abs(fit.plots[[i]]$ratio$Yrot - fit.plots[[target]][["ratio"]]))
  mean.dev[[i]][["braycurtis"]] <- mean(abs(fit.plots[[i]]$braycurtis$Yrot - fit.plots[[target]][["braycurtis"]]))
  
  max.dev[[i]] <- list()
  max.dev[[i]][["unweighted"]] <- max(abs(fit.plots[[i]]$unweighted$Yrot - fit.plots[[target]][["unweighted"]]))
  max.dev[[i]][["weighted"]] <- max(abs(fit.plots[[i]]$weighted$Yrot - fit.plots[[target]][["weighted"]]))
  max.dev[[i]][["information"]] <- max(abs(fit.plots[[i]]$information$Yrot - fit.plots[[target]][["information"]]))
  max.dev[[i]][["ratio"]] <- max(abs(fit.plots[[i]]$ratio$Yrot - fit.plots[[target]][["ratio"]]))
  max.dev[[i]][["braycurtis"]] <- max(abs(fit.plots[[i]]$braycurtis$Yrot - fit.plots[[target]][["braycurtis"]]))
	
	minRange[["unweighted"]] <- min(minRange[["unweighted"]], min(fit.plots[[i]]$unweighted$Yrot))
	minRange[["weighted"]] <- min(minRange[["weighted"]], min(fit.plots[[i]]$weighted$Yrot))
	minRange[["information"]] <- min(minRange[["information"]], min(fit.plots[[i]]$information$Yrot))
	minRange[["ratio"]] <- min(minRange[["ratio"]], min(fit.plots[[i]]$ratio$Yrot))
	minRange[["braycurtis"]] <- min(minRange[["braycurtis"]], min(fit.plots[[i]]$braycurtis$Yrot))

	maxRange[["unweighted"]] <- max(maxRange[["unweighted"]], max(fit.plots[[i]]$unweighted$Yrot))
	maxRange[["weighted"]] <- max(maxRange[["weighted"]], max(fit.plots[[i]]$weighted$Yrot))
	maxRange[["information"]] <- max(maxRange[["information"]], max(fit.plots[[i]]$information$Yrot))
	maxRange[["ratio"]] <- max(maxRange[["ratio"]], max(fit.plots[[i]]$ratio$Yrot))
	maxRange[["braycurtis"]] <- max(maxRange[["braycurtis"]], max(fit.plots[[i]]$braycurtis$Yrot))
}

# get mean and max into vectors for plotting
unweighted.mean <- unlist(lapply(mean.dev,function(x) { return(x[["unweighted"]]) }))
weighted.mean <- unlist(lapply(mean.dev,function(x) { return(x[["weighted"]]) }))
information.mean <- unlist(lapply(mean.dev,function(x) { return(x[["information"]]) }))
ratio.mean <- unlist(lapply(mean.dev,function(x) { return(x[["ratio"]]) }))
braycurtis.mean <- unlist(lapply(mean.dev,function(x) { return(x[["braycurtis"]]) }))

unweighted.max <- unlist(lapply(max.dev,function(x) { return(x[["unweighted"]]) }))
weighted.max <- unlist(lapply(max.dev,function(x) { return(x[["weighted"]]) }))
information.max <- unlist(lapply(max.dev,function(x) { return(x[["information"]]) }))
ratio.max <- unlist(lapply(max.dev,function(x) { return(x[["ratio"]]) }))
braycurtis.max <- unlist(lapply(max.dev,function(x) { return(x[["braycurtis"]]) }))

normalize.on.range <- function(myvalues,minRange,maxRange,name) {
	return( (myvalues)/abs(maxRange[[name]]-minRange[[name]]))
}

unweighted.mean <- normalize.on.range(unweighted.mean, minRange, maxRange, "unweighted")
weighted.mean <- normalize.on.range(weighted.mean, minRange, maxRange, "weighted")
information.mean <- normalize.on.range(information.mean, minRange, maxRange, "information")
ratio.mean <- normalize.on.range(ratio.mean, minRange, maxRange, "ratio")
braycurtis.mean <- normalize.on.range(braycurtis.mean, minRange, maxRange, "braycurtis")

unweighted.max <- normalize.on.range(unweighted.max, minRange, maxRange, "unweighted")
weighted.max <- normalize.on.range(weighted.max, minRange, maxRange, "weighted")
information.max <- normalize.on.range(information.max, minRange, maxRange, "information")
ratio.max <- normalize.on.range(ratio.max, minRange, maxRange, "ratio")
braycurtis.max <- normalize.on.range(braycurtis.max, minRange, maxRange, "braycurtis")

# plot ellipses
pdf("tongue_dorsum_ellipses_log_ratio.pdf")
xMax <- max(c(unweighted.mean, weighted.mean, information.mean, ratio.mean, braycurtis.mean))
xMax <- xMax * 1.2
yMax <- max(c(unweighted.max, weighted.max, information.max, ratio.max, braycurtis.max))
yMax <- yMax * 1.2
plot(c(0,xMax),c(0,yMax),col="white",main="Max vs. mean relative deviation of rarefactions",xlab="Maximum relative deviation per rarefaction",ylab="Mean relative deviation per rarefaction")
dataEllipse(unweighted.mean, unweighted.max,xlim=c(0,xMax), ylim=c(0,yMax), plot.points=FALSE, col = "red", levels=c(0.95), center.pch=0)
par(new=TRUE)
dataEllipse(weighted.mean, weighted.max,xlim=c(0,xMax), ylim=c(0,yMax), plot.points=FALSE, col = "green", levels=c(0.95), center.pch=0)
par(new=TRUE)
dataEllipse(information.mean, information.max,xlim=c(0,xMax), ylim=c(0,yMax), plot.points=FALSE, col = "blue", levels=c(0.95), center.pch=0)
par(new=TRUE)
dataEllipse(ratio.mean, ratio.max,xlim=c(0,xMax), ylim=c(0,yMax), plot.points=FALSE, col = "black", levels=c(0.95), center.pch=0)
par(new=TRUE)
dataEllipse(braycurtis.mean, braycurtis.max,xlim=c(0,xMax), ylim=c(0,yMax), plot.points=FALSE, col = "purple", levels=c(0.95), center.pch=0)
par(new=TRUE)

dev.off()
