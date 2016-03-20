#!/usr/bin/env Rscript 
options(error=recover)

library(ape)
library(phangorn)
library(vegan)

# commenting out incorrect clr dirichlet

# get default par
plotParameters <- par()

source("../UniFrac.r")

# read OTU table and format appropriately for input into UniFrac methods
breastmilk.otu.tab <- read.table("./data/td_OTU_tag_mapped_lineage.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

#remove taxonomy column to make otu count matrix numeric
taxonomy <- breastmilk.otu.tab$taxonomy
breastmilk.otu.tab <- breastmilk.otu.tab[-length(colnames(breastmilk.otu.tab))]
breastmilk.otu.tab <- t(as.matrix(breastmilk.otu.tab))

#sort taxa from most to least abundant
taxaOrder <- rev(order(apply(breastmilk.otu.tab,2,sum)))
taxonomy <- taxonomy[taxaOrder]
breastmilk.otu.tab <- breastmilk.otu.tab[,taxaOrder]
breastmilk.otu.tab.rarefy <- rrarefy(breastmilk.otu.tab, min(apply(breastmilk.otu.tab,1,sum)))

# read and root tree (rooted tree is required)
breastmilk.tree <- read.tree("./data/fasttree_all_seed_OTUs.tre")
breastmilk.tree <- midpoint(breastmilk.tree)

# read metadata
MyMeta<- read.table("./data/metadata.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

#remove infected sample S38I
#MyMeta <- MyMeta[(which(rownames(MyMeta)!="S38I")),]

# filter OTU table and metadata so that only samples which appear in both are retained
otu_indicies <- match(rownames(MyMeta),rownames(breastmilk.otu.tab))
otu_indicies <- otu_indicies[!is.na(otu_indicies)]
breastmilk.otu.tab <- breastmilk.otu.tab[otu_indicies,]
MyMetaOrdered <- MyMeta[match(rownames(breastmilk.otu.tab),rownames(MyMeta)),]

unweighted <- getDistanceMatrix(breastmilk.otu.tab.rarefy,breastmilk.tree,method="unweighted",verbose=TRUE)

all_distance_matrices <- getDistanceMatrix(breastmilk.otu.tab,breastmilk.tree,method="all",verbose=TRUE)

weighted <- all_distance_matrices[["weighted"]]
information <- all_distance_matrices[["information"]]
ratio_no_log <- all_distance_matrices[["ratio_no_log"]]

groups <- rep("Not Infected",length(MyMetaOrdered$Gestation))
groups[which(rownames(MyMetaOrdered)=="S38I")] <- "Infected"
groups <- as.factor(groups)

otuSum <- apply(breastmilk.otu.tab,1,sum)

# caculate pcoa vectors
unweighted.pcoa <- pcoa(unweighted)
weighted.pcoa <- pcoa(weighted)
information.pcoa <- pcoa(information)
ratio_no_log.pcoa <- pcoa(ratio_no_log)

# calculate total variance explained
unweighted.varExplained <- sum(apply(unweighted.pcoa$vector,2,function(x) sd(x)*sd(x)))
weighted.varExplained <- sum(apply(weighted.pcoa$vector,2,function(x) sd(x)*sd(x)))
information.varExplained <- sum(apply(information.pcoa$vector,2,function(x) sd(x)*sd(x)))
ratio_no_log.varExplained <- sum(apply(ratio_no_log.pcoa$vector,2,function(x) sd(x)*sd(x)))

# calculate proportion of variance explained by First Coordinate
unweighted.pc1.varEx <- sd(unweighted.pcoa$vector[,1])*sd(unweighted.pcoa$vector[,1])/unweighted.varExplained
#calculate proportion of variance explained by Second Coordinate
unweighted.pc2.varEx <- sd(unweighted.pcoa$vector[,2])*sd(unweighted.pcoa$vector[,2])/unweighted.varExplained

weighted.pc1.varEx <- sd(weighted.pcoa$vector[,1])*sd(weighted.pcoa$vector[,1])/weighted.varExplained
weighted.pc2.varEx <- sd(weighted.pcoa$vector[,2])*sd(weighted.pcoa$vector[,2])/weighted.varExplained

information.pc1.varEx <- sd(information.pcoa$vector[,1])*sd(information.pcoa$vector[,1])/information.varExplained
information.pc2.varEx <- sd(information.pcoa$vector[,2])*sd(information.pcoa$vector[,2])/information.varExplained

ratio_no_log.pc1.varEx <- sd(ratio_no_log.pcoa$vector[,1])*sd(ratio_no_log.pcoa$vector[,1])/ratio_no_log.varExplained
ratio_no_log.pc2.varEx <- sd(ratio_no_log.pcoa$vector[,2])*sd(ratio_no_log.pcoa$vector[,2])/ratio_no_log.varExplained

#save plots as PDF
pdf("output/breastmilk_pcoa_plots.pdf")

# MAKE BAR PLOTS

#convert to dist structure
unweighted.dist <- as.dist(unweighted)
weighted.dist <- as.dist(weighted)
information.dist <- as.dist(information)

#"average" is most similar to UPGMA, apparently
unweighted.dendo <- hclust(unweighted.dist, method="average")
weighted.dendo <- hclust(weighted.dist, method="average")
information.dendo <- hclust(information.dist, method="average")

#get otu proportions for barplot
prop <- t(apply(breastmilk.otu.tab,1,function(x) x/sum(x)))

# plot dendogram with bar plots

#GG legacy code. Fix size, margins, position
par(mfrow=c(2,1), mar=c(1, 3, 2, 1) + 0.1,cex=0.3)

plot(unweighted.dendo, axes=F, ylab=NULL, ann=F, hang=-1)

#order the barplot 
colors <- c("steelblue3","skyblue1", "indianred1", "mediumpurple1", "olivedrab3", "pink", "#FFED6F", "mediumorchid3", "ivory2", "tan1", "aquamarine3", "#C0C0C0", "royalblue4", "mediumvioletred", "#999933", "#666699", "#CC9933", "#006666", "#3399FF", "#993300", "#CCCC99", "#666666", "#FFCC66", "#6699CC", "#663366", "#9999CC", "#CCCCCC", "#669999", "#CCCC66", "#CC6600", "bisque", "#9999FF", "#0066CC", "#99CCCC", "#999999", "#FFCC00", "#009999", "#FF9900", "#999966", "#66CCCC", "#339966", "#CCCC33", "#EDEDED")
barplot(t(prop[unweighted.dendo$order,]), space=0,col=colors, las=2)

plot(weighted.dendo, axes=F, ylab=NULL, ann=F, hang=-1)
#order the barplot 
barplot(t(prop[weighted.dendo$order,]), space=0,col=colors, las=2)

plot(information.dendo, axes=F, ylab=NULL, ann=F, hang=-1)
#order the barplot 
barplot(t(prop[information.dendo$order,]), space=0,col=colors, las=2)


par(plotParameters)



#choose colors for each condition
palette(c("red","black","cyan","dodgerblue","blue","orange"))

#plot pcoa plots with legend
plot(unweighted.pcoa$vectors[,1],unweighted.pcoa$vectors[,2], type="p",col=groups,main="Unweighted UniFrac\nprincipal coordinate analysis",xlab=paste("First Coordinate", round(unweighted.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Coordinate", round(unweighted.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
#placement with S38I included
#legend(-0.1,-0.055,levels(groups),col=palette(),pch=19)
# #placement with S38I excluded
# legend(0.055,0.15,levels(groups),col=palette(),pch=19)

plot(weighted.pcoa$vectors[,1],weighted.pcoa$vectors[,2], col=groups,main="Weighted UniFrac\nprincipal coordinate analysis",xlab=paste("First Coordinate", round(weighted.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Coordinate", round(weighted.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
#legend(-0.1,-0.12,levels(groups),col=palette(),pch=19)

plot(information.pcoa$vectors[,1],information.pcoa$vectors[,2], col=groups,main="Information UniFrac\nprincipal coordinate analysis",xlab=paste("First Coordinate", round(information.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Coordinate", round(information.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
#placement with S38I included
#legend(-0.15,-0.4,levels(groups),col=palette(),pch=19)
# #placement with S38I excluded
# legend(0.4,-0.15,levels(groups),col=palette(),pch=19)

plot(ratio_no_log.pcoa$vectors[,1],ratio_no_log.pcoa$vectors[,2], col=groups,main="Centered Ratio UniFrac\nprincipal coordinate analysis",xlab=paste("First Coordinate", round(ratio_no_log.pc1.varEx,digits=3),"variance explained"),ylab=paste("Second Coordinate", round(ratio_no_log.pc2.varEx,digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)

dev.off()
