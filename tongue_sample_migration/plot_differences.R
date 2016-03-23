## FOR GG
# Plots a set of 30 tongue samples vs 30 tongue samples
# Should be no differences, but there's a difference in the sets of 

d <- read.table("Rare1_659/bdiv_all_r659_original/unweighted_unifrac_pc.txt", header=T, row.names=1, sep="\t")
e <- read.table("Rare2_659/bdiv_all_r659_original/unweighted_unifrac_pc.txt", header=T, row.names=1, sep="\t")

# Rudimentary procrustes analysis 
# Correct for the shape (eigenvalues) (ONLY FOR DATASETS THAT ARE FLIPPED ABOUT X AND Y AXIS
# i.e. for the 100 dataset leave lines 9 and 10 as is, for 200 or 3, comment them out
e[,1] <- (e[,1]*-1)
e[,2] <- (e[,2]*-1)




# Setup axis labelling conventions
per <- "% Explained"
x.pc1.explained <- d[nrow(d),1]
x.pc1.explained <- round(x.pc1.explained, digits=1)
xlabel <- "Principal Component1 Eigenvalues: "
xla <- paste(xlabel, x.pc1.explained,per, sep="", collapse=NULL)

x.pc2.explained <- d[nrow(d),2]
x.pc2.explained <- round(x.pc2.explained, digits=1)
ylabel <- "Principal Component2 Eigenvalues: "
yla <- paste(ylabel, x.pc2.explained,per, sep="", collapse=NULL)

rtl <- NULL
ltr <- NULL

component <- 1

for (i in 1:60)	# Iterate through all 60 samples
{
	if ( (d[i,component] > 0) & (e[i,component] < 0))
	{	
		print("Right Cluster to Left Cluster Movement: ")	# Moved from right to left
		print(rownames(d[i,]))
		rtl <- c(rtl, rownames(d[i,]))	# Get list of Right to Left Movement
	}	
}
for (i in 1:60)
{
	if ( (d[i,component] < 0) & (e[i,component] > 0))
	{
		print("Left Cluster to Right Cluster Movement: ")	# Moved from left to right
		print(rownames(d[i,]))
		ltr <- c(ltr, rownames(d[i,]))	# Get list of Left to Right movement
	}
}

# First/second refer to the columns the data is being read from 
first <- 1
second <- 2

pdf("UniFrac_tvst_movement.pdf")	# Comment out if not plotting
plot(d[1:60,first], d[1:60,second],main="Sample Migration between Rarefactions",xlab=xla,ylab=yla, pch=19, col=rgb(0,0,0,0.4))	# Plot the 1st rarefation

shapes <- c(19,15,19,15,19)
colours <- c(rgb(1,0,0,0.4), rgb(1,0,0,0.4), rgb(0,0,1,0.4), rgb(0,0,1,0.4), rgb(0,0,0,0.4) )

#legend(x=-0.25,y=0.31,title="Sample Movement", legend = c("Left->Right","Origin","Right->Left","Origin", "No Change"), pch=shapes, col=colours)

# Points that have moved from right to left are red
# Points that have moved from left to right are blue

# Points represent a plot of Rarefaction 1, but are coloured based on sample movement of Rare 2
points(d[rtl,first][e[rtl,first] < 0], d[rtl,second][e[rtl,first] < 0], pch=19, col=rgb(1,0,0,0.4))	# Red
points(d[ltr,first][e[ltr,first] > 0], d[ltr,second][e[ltr,first] > 0], pch=19, col=rgb(0,0,1,0.4))	# Blue

# The squares indicate WHERE they have moved from 
# Can comment this out to remove square drawing
points(e[rtl,first], e[rtl,second], pch = 15, col=rgb(1,0,0,0.4))
points(e[ltr,first], e[ltr,second], pch = 15, col=rgb(0,0,1,0.4))

# Draw line segments between the movements
arrows(e[rtl,first], e[rtl, second], d[rtl, first], d[rtl, second], col = rgb(1,0,0,0.5), length=0.1)

arrows(e[ltr,first], e[ltr, second], d[ltr, first], d[ltr, second], col = rgb(0,0,1,0.5), length=0.1)



dev.off()







