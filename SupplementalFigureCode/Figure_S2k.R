#Figure S2k pairwise correlation plots for each cluster across animals.

#Library needed:
library(corrplot)

###First you must run: 
#source("Cperk_TimeCourse_DEGanalyses.R")
#source("Figure_2D.R")

#look at clusters vs. nottimeDEs
Cperk_timeOGs<-unique(na.omit(Cperk_OGs2[match(timeDEs,Cperk_OGs2[,2]),1]))
Cperk_notTimeOGs<-setdiff(Cperk_OGs2[,2],Cperk_timeOGs)

# Define your desired cluster order
desired_cluster_order <- c(4,3,1,2,5)

GOogs.nt<-unique(Cperk_notTimeOGs) ###not DE expressed in development
Cperk.k1<-names(kClusters[kClusters==4]) 
GOogs.k1<-unique(na.omit(Cperk_OGs2[match(Cperk.k1,Cperk_OGs2[,2]),1])) ##Cluster 1
Cperk.k2<-names(kClusters[kClusters==3])
GOogs.k2<-unique(na.omit(Cperk_OGs2[match(Cperk.k2,Cperk_OGs2[,2]),1])) ##Cluster 2
Cperk.k3<-names(kClusters[kClusters==1])
GOogs.k3<-unique(na.omit(Cperk_OGs2[match(Cperk.k3,Cperk_OGs2[,2]),1])) ##Cluster 3
Cperk.k4<-names(kClusters[kClusters==2])
GOogs.k4<-unique(na.omit(Cperk_OGs2[match(Cperk.k4,Cperk_OGs2[,2]),1])) ##Cluster 4
Cperk.k5<-names(kClusters[kClusters==5])
GOogs.k5<-unique(na.omit(Cperk_OGs2[match(Cperk.k5,Cperk_OGs2[,2]),1])) ##Cluster 5

###define which cluster you would like to plot? nt, k1, k2, k3, k4, k5 
GOogs<-GOogs.nt

# Assuming OG_matList is a list of lists of matrices, and we're focusing on the first list of matrices
matrices <- OG_matList[[1]]  # Extract the first list of matrices

###lets make sure to look at all patterns per OG!
base_ids <- ifelse(grepl("\\.", GOogs), sub("\\..*", "", GOogs), GOogs)
patterns <- paste0("^", base_ids, "(\\.\\d+)?$")
OGpats <- unlist(lapply(patterns, function(pat) grep(pat, rownames(matrices[[1]]), value = TRUE)))

# Subset all matrices by the selected rownames
matrices2 <- lapply(matrices, function(x) x[OGpats, ])

# Function to flatten each matrix into a single column
flatten_matrix <- function(mat) {
  return(as.vector(t(mat)))  # Transpose and convert to vector to stack rows
}

# Apply the function to each matrix and combine them into a single matrix
flattened_matrices <- do.call(cbind, lapply(matrices2, flatten_matrix))
#flattened_matrices<-na.omit(flattened_matrices)
#flattened_matrices[is.na(flattened_matrices)]<-0#rnorm(length(which(is.na(flattened_matrices))),0,1)


myorder<-c("PT","FW","CE","HD","DM","SU","ZF","NV","AQ","ML","Cperk")

flattened_matrices<-flattened_matrices[,myorder]

myCorDat<- Hmisc::rcorr(flattened_matrices, type="pearson")

#myCorDat<-rcorr.test(na.omit(flattened_matrices), plot = TRUE, table = TRUE, var.names = NULL,scale.font = 1)

# Create a copy of the correlation matrix and set the diagonal to NA
modified_corr_matrix <- myCorDat$r
diag(modified_corr_matrix) <- NA

#Figure S2k--one plot. Change the parameter in line 27 to plot other cluster correlations: GOogs<-GOogs.nt
# Plot the correlation matrix using the calculated range without diagonal values
corrplot(corr=modified_corr_matrix, p.mat = myCorDat$P, is.corr = FALSE, col.lim=c(-0.5,0.5), method = 'color', col = COL2('PiYG'), cl.pos = 'b', diag=F, type = 'lower',tl.srt = 45,order="original",sig.level = 0.05, addgrid.col = 'white')

