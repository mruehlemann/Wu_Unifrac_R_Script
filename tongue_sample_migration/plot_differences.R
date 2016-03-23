## FOR GG
# Plots a set of 30 tongue samples vs 30 tongue samples
# Should be no differences, but there's a difference in the sets of 
library(vegan)
source("../UniFrac.r")

# read data
original.data <- read.table("../data/tongue_dorsum/tongue_vs_tongue_30_forR.txt",sep="\t",check.names=FALSE,quote="",comment.char="", header=TRUE,row.names=1)
tree <- read.tree("../data/tongue_dorsum_vs_buccal_mucosa/hmp_tongue_cheek_subtree.tre")

# remove all OTUs with zero counts across all samples (this data set should have a min. of 100 counts per OTU, but this code is here just in case.)
otu.sum <- apply(original.data,1,sum)
original.data <- original.data[which(otu.sum > 0),]
otu.sum <- otu.sum[which(otu.sum>0)]

# make sure tree tip names match OTU names by taking out single quotes
tree$tip.label <- gsub("'","",tree$tip.label)
# remove OTUs that are not in the tree
original.data <- original.data[which(rownames(original.data) %in% tree$tip.label),]

original.data <- t(original.data)

# remove extra taxa from tree
absent <- tree$tip.label[!(tree$tip.label %in% colnames(original.data))]
if (length(absent) != 0) {
		tree <- drop.tip(tree, absent)
}

# root tree (rooted tree is required)
tree <- midpoint(tree)

d.data <- rrarefy(original.data, min(otu.sum))
e.data <- rrarefy(original.data, min(otu.sum))

d.otu.sum <- apply(d.data,2,sum)
d.data <- d.data[,which(d.otu.sum > 0)]
d.tree <- tree
absent <- d.tree$tip.label[!(d.tree$tip.label %in% colnames(d.data))]
if (length(absent) != 0) {
		d.tree <- drop.tip(d.tree, absent)
}

e.otu.sum <- apply(e.data,2,sum)
e.data <- e.data[,which(e.otu.sum > 0)]
e.tree <- tree
absent <- e.tree$tip.label[!(e.tree$tip.label %in% colnames(e.data))]
if (length(absent) != 0) {
		e.tree <- drop.tip(e.tree, absent)
}

d.unifrac <- getDistanceMatrix(d.data,d.tree,method="unweighted",verbose=TRUE)
e.unifrac <- getDistanceMatrix(e.data,e.tree,method="unweighted",verbose=TRUE)

d <- pcoa(d.unifrac)
e <- pcoa(e.unifrac)

#function to get variance explained for the PCOA component labels
getVarExplained <- function(vector) {
	rawVarEx <- apply(vector,2,function(x) sd(x)*sd(x))
	totalVarExplained <- sum(rawVarEx)
	varEx <- rawVarEx/totalVarExplained
	return(varEx)
}


d.varEx <- getVarExplained(d$vectors)

# Setup axis labelling conventions
per <- "% Explained"
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

pdf("UniFrac_tvst_movement.pdf")	# Comment out if not plotting
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



dev.off()







