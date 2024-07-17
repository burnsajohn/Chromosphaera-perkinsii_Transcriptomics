###run correlation analyses on "basals" as defined however you would like. In this case, "basals" is C. perkinsii, NV, AQ, ML. Other options are correlations between all animals, or coreelations between just the "basal" animals without C.perkinsii.

#libraries needed:
library(topGO)
library(plyr)
library(dendextend)
library(pheatmap)

###First you must run: 
#1. Cperk_TimeCourse_DEGanalyses.R
#2. Figure_2D.R

basals<-OG_matList[[1]][c("Cperk","NV","AQ","ML")]
#basals<-OG_matList[[1]][1:11]
#basals<-OG_matList[[1]][c("NV","AQ","ML")]

#basals<-OG_matList[[1]]
# Assuming basals is your list of matrices
n <- length(basals)
if (n < 2) {
  stop("Need at least two matrices for correlation")
}

# Initialize a vector to store the average correlations
average_correlations <- numeric(nrow(basals[[1]]))

####compute pairwise correlations across every row for all pairs of animals in "basals":
nrow(basals[[1]])
# Loop over each row
for (i in 1:nrow(basals[[1]])) {
  print(i)
  correlations <- c()

  # Compute pairwise correlations
  for (j in 1:(n-1)) {
    for (k in (j+1):n) {
      row1 <- basals[[j]][i, ]
      row2 <- basals[[k]][i, ]
	  row1[is.na(row1)] <- 0
	  row2[is.na(row2)] <- 0
	  #print(row1)
	  #print(row2)
	  #print(diss.CDM(row1,row2))
      correlations <- c(correlations, cor(row1, row2,method="pearson"))
    }
  }

  # Compute the average correlation for the row
  average_correlations[i] <- mean(correlations, na.rm = TRUE)
}

# Set the names of the average_correlations vector to the row names
names(average_correlations) <- rownames(basals[[1]])

# The result is in average_correlations
#average_correlations

#Figure S2G
hist(na.omit(average_correlations),breaks=100,main="Early diverging animals plus C. perkinsii average genewise correlations",xlab="Correlation coefficient")
abline(v=0.5,col="blue",lty=2,lwd=2)
abline(v=-0.5,col="red",lty=2,lwd=2)
###

###collect highly correlated genes:
highCor<-as.character(na.omit(names(average_correlations[(average_correlations)>0.5])))
highCorBaseID <- unique(gsub("\\..*", "", highCor))
Cperk_highCorGenes<-Cperk_OGs2[match(highCorBaseID,Cperk_OGs2[,1]),]
Cperk_highCorGenes<-Cperk_highCorGenes[Cperk_highCorGenes[,2] %in% timeDEs,]

