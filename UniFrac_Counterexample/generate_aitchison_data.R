otus <- t(read.table("RANDOM_otus_tongue_vs_tongue_60.txt", sep="\t", header=T, row.names=1))
rownames(otus) <- gsub("^X","", rownames(otus))	# Clean up formatting
otus <- as.data.frame(otus)
clr.otus <- aldex.clr(otus, mc.samples=128)

clr.otu.inst <- sapply(getMonteCarloInstances(clr.otus), function(y){y[,1]})
ait.dist <- dist(clr.otu.inst)
ait.matrix <- as.matrix(ait.dist)
ait.p <- pcoa(ait.matrix)
ait.pcoa <- ait.p[4]
ait.pcoa <- as.data.frame(ait.pcoa)
a.1 <- ait.pcoa$vectors.Axis.1

folder <- "ac_all_"

name <- "/atc_pc"
ext <- ".txt"

for (i in 1:10)
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