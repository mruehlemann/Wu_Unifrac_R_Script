#!/usr/bin/env Rscript 
options(error=recover)

library(ape)
library(phangorn)
library(vegan)
library(randomcoloR)

# commenting out incorrect clr dirichlet

# get default par
plotParameters <- par()

# read OTU table and format appropriately for input into UniFrac methods
breastmilk.otu.tab <- read.table("../../data/breastmilk/td_OTU_tag_mapped_lineage.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

#remove taxonomy column to make otu count matrix numeric
taxonomy <- breastmilk.otu.tab$taxonomy
breastmilk.otu.tab <- breastmilk.otu.tab[-length(colnames(breastmilk.otu.tab))]
breastmilk.otu.tab <- t(as.matrix(breastmilk.otu.tab))

#sort taxa from most to least abundant
taxaOrder <- rev(order(apply(breastmilk.otu.tab,2,sum)))
taxonomy <- taxonomy[taxaOrder]
breastmilk.otu.tab <- breastmilk.otu.tab[,taxaOrder]

# extract genus + species name from taxonomy
taxonomy <- as.character(taxonomy)
for (i in c(1:length(taxonomy))) {
  taxonomy[i] <- paste(strsplit(taxonomy[i],c(";"))[[1]][6],strsplit(taxonomy[i],c(";"))[[1]][7])
}

# rarefy before performing unweighted UniFrac
breastmilk.otu.tab.rarefy <- rrarefy(breastmilk.otu.tab, min(apply(breastmilk.otu.tab,1,sum)))

# read and root tree (rooted tree is required)
breastmilk.tree <- read.tree("../../data/breastmilk/fasttree_all_seed_OTUs.tre")
breastmilk.tree <- midpoint(breastmilk.tree)

# read metadata
MyMeta<- read.table("../../data/breastmilk/metadata.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

#remove infected sample S38I
#MyMeta <- MyMeta[(which(rownames(MyMeta)!="S38I")),]

# filter OTU table and metadata so that only samples which appear in both are retained
otu_indicies <- match(rownames(MyMeta),rownames(breastmilk.otu.tab))
otu_indicies <- otu_indicies[!is.na(otu_indicies)]
breastmilk.otu.tab <- breastmilk.otu.tab[otu_indicies,]
MyMetaOrdered <- MyMeta[match(rownames(breastmilk.otu.tab),rownames(MyMeta)),]

infected.sample <- breastmilk.otu.tab[c("S38I"),]
names(infected.sample) <- taxonomy

# permute taxa from the infected sample to generate random monocultures
monocultures <- data.frame(matrix(nrow=60,ncol=length(infected.sample)))
rownames(monocultures) <- c(paste("Pasteurella",c(1:20),sep="_"),paste("Staphylococcus",c(1:20),sep="_"),paste("Pseudomonas",c(1:20),sep="_"))
colnames(monocultures) <- taxonomy
for (i in c(1:60)) {
  permutation <- sample(c(1:length(infected.sample)),length(infected.sample))
  monocultures[i,] <- infected.sample[permutation]
}

# make 20 of the monocultures Pasteurella dominated, 20 Staphylococcus dominated, and 20 Pseudomonas dominated (these are the taxa with the top abundances in this data set)
max.otu <- max(monocultures[1,])
for (i in c(1:20)) {
  monocultures[i,which(monocultures[i,]==max.otu)] <- monocultures[i,"Pasteurella |93"]
  monocultures[i,"Pasteurella |93"] <- max.otu
}
for (i in c(21:40)) {
  monocultures[i,which(monocultures[i,]==max.otu)] <- monocultures[i,"Staphylococcus |77" ]
  monocultures[i,"Staphylococcus |77" ] <- max.otu
}
for (i in c(41:60)) {
  monocultures[i,which(monocultures[i,]==max.otu)] <- monocultures[i,"Pseudomonas |92"]
  monocultures[i,"Pseudomonas |92"] <- max.otu
}

#flip table so that sample names are in column headers, to be consistent with other data sets
monocultures <- t(monocultures)
write.table(monocultures,file="monocultures.txt",sep="\t",quote=FALSE,row.names=TRUE)