library(ape)
library(phangorn)
library(vegan)
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
  unifrac[[i]][["ratio"]] <- all_distance_matrices[["ratio_no_log"]]
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
for (i in c(2:100)) {
  mean.dev[[i]] <- list()
  mean.dev[[i]][["unweighted"]] <- mean(abs(fit.plots[[i]][["unweighted"]] - fit.plots[[target]][["unweighted"]]))
  mean.dev[[i]][["weighted"]] <- mean(abs(fit.plots[[i]][["weighted"]] - fit.plots[[target]][["weighted"]]))
  mean.dev[[i]][["information"]] <- mean(abs(fit.plots[[i]][["information"]] - fit.plots[[target]][["information"]]))
  mean.dev[[i]][["ratio"]] <- mean(abs(fit.plots[[i]][["ratio"]] - fit.plots[[target]][["ratio"]]))
  mean.dev[[i]][["braycurtis"]] <- mean(abs(fit.plots[[i]][["braycurtis"]] - fit.plots[[target]][["braycurtis"]]))
  
  max.dev[[i]] <- list()
  max.dev[[i]][["unweighted"]] <- max(abs(fit.plots[[i]][["unweighted"]] - fit.plots[[target]][["unweighted"]]))
  max.dev[[i]][["weighted"]] <- max(abs(fit.plots[[i]][["weighted"]] - fit.plots[[target]][["weighted"]]))
  max.dev[[i]][["information"]] <- max(abs(fit.plots[[i]][["information"]] - fit.plots[[target]][["information"]]))
  max.dev[[i]][["ratio"]] <- max(abs(fit.plots[[i]][["ratio"]] - fit.plots[[target]][["ratio"]]))
  max.dev[[i]][["braycurtis"]] <- max(abs(fit.plots[[i]][["braycurtis"]] - fit.plots[[target]][["braycurtis"]]))
}

# get mean and max into vectors for plotting
unweighted.mean <- unlist(lapply(mean.dev),function(x) { return(x[["unweighted"]]) })
weighted.mean <- unlist(lapply(mean.dev),function(x) { return(x[["weighted"]]) })
information.mean <- unlist(lapply(mean.dev),function(x) { return(x[["information"]]) })
ratio.mean <- unlist(lapply(mean.dev),function(x) { return(x[["ratio"]]) })
braycurtis.mean <- unlist(lapply(mean.dev),function(x) { return(x[["braycurtis"]]) })

unweighted.max <- unlist(lapply(max.dev),function(x) { return(x[["unweighted"]]) })
weighted.max <- unlist(lapply(max.dev),function(x) { return(x[["weighted"]]) })
information.max <- unlist(lapply(max.dev),function(x) { return(x[["information"]]) })
ratio.max <- unlist(lapply(max.dev),function(x) { return(x[["ratio"]]) })
braycurtis.max <- unlist(lapply(max.dev),function(x) { return(x[["braycurtis"]]) })

# plot ellipses
pdf("tongue_dorsum_ellipses.pdf")
xMax <- max(c(unweighted.mean, weighted.mean, information.mean, ratio.mean, braycurtis.mean))
yMax <- max(c(unweighted.max, weighted.max, information.max, ratio.max, braycurtis.max))

dataEllipse(unweighted.mean, unweighted.max,xlim=c(0,xMax), ylim=c(0,yMax), plot.points=FALSE, col = palette()[i], levels=c(0.95), center.pch=0)
par(new=TRUE)
dataEllipse(weighted.mean, weighted.max,xlim=c(0,xMax), ylim=c(0,yMax), plot.points=FALSE, col = palette()[i], levels=c(0.95), center.pch=0)
par(new=TRUE)
dataEllipse(information.mean, information.max,xlim=c(0,xMax), ylim=c(0,yMax), plot.points=FALSE, col = palette()[i], levels=c(0.95), center.pch=0)
par(new=TRUE)
dataEllipse(ratio.mean, ratio.max,xlim=c(0,xMax), ylim=c(0,yMax), plot.points=FALSE, col = palette()[i], levels=c(0.95), center.pch=0)
par(new=TRUE)
dataEllipse(braycurtis.mean, braycurtis.max,xlim=c(0,xMax), ylim=c(0,yMax), plot.points=FALSE, col = palette()[i], levels=c(0.95), center.pch=0)
par(new=TRUE)

dev.off()
