#Figure S2j empirical statistics calculation and plotting.

###First you must run: 
#source("Cperk_TimeCourse_DEGanalyses.R")
#source("Figure_2D.R")

################################################################################################
####empirical stats#############################################################################
################################################################################################
#10,000 random samples
cormeds<-vector(length=10000)

for(stat in 1:10000){
print(stat)
GOogs<-sample(allOGs,length(devGOogs))
matrices <- OG_matList[[1]]  # Extract the first list of matrices
###lets make sure to look at all patterns per OG!
base_ids <- ifelse(grepl("\\.", GOogs), sub("\\..*", "", GOogs), GOogs)
patterns <- paste0("^", base_ids, "(\\.\\d+)?$")
OGpats <- unlist(lapply(patterns, function(pat) grep(pat, rownames(matrices[[1]]), value = TRUE)))

if(length(OGpats)<2){next}

# Subset all matrices by the selected rownames
basals <- lapply(matrices, function(x) x[OGpats, ])


# Initialize a vector to store the average correlations
correlations <- c()

####compute pairwise correlations across every row for Cperk vs. each animal in "basals":
nrow(basals[[1]])
# Loop over each row
for (i in 1:nrow(basals[[1]])) {
  #print(i)
  # Compute pairwise correlations
  for (j in 2:11) {
      row1 <- basals[[1]][i, ]
      row2 <- basals[[j]][i, ]
	  if(all(is.na(row1)) || all(is.na(row2))){
	    #print(row1)
		#print(row2)
		correlations <- c(correlations, NA)
	  }else{
	  #print(row1)
	  #print(row2)
	  #print(cor(row1, row2,method="pearson"))
	  #print(diss.CDM(row1,row2))
	  if(cor(row1, row2,method="pearson")==0){
	  print(row1)
	  print(row2)
	  }
      correlations <- c(correlations, cor(row1, row2,method="pearson"))
	  }
  }
}

cormeds[stat]<-median(correlations,na.rm=T)
}

#boxplot(correlations)


###single_sample

Cperk_develGenes<-c("Nk52_evm15s1360","Nk52_evm78s1737","Nk52_evm1s1535","Nk52_evm57s1073","Nk52_evm44s1360","Nk52_evm60s152","Nk52_evm4s2657","Nk52_evm64s2039","Nk52_evm98s2192","Nk52_evm16s267","Nk52_evm7s490","Nk52_evm24s1967","Nk52_evm90s1737","Nk52_evm65s270","Nk52_evm36s621","Nk52_evm4s293","Nk52_evm20s2474","Nk52_evm62s1020","Nk52_evm18s1967","Nk52_evm48s153","Nk52_evm11s256","Nk52_evm45s1569","Nk52_evm29s2367","Nk52_evm50s352","Nk52_evm4s217","Nk52_evm19s307","Nk52_evm18s164","Nk52_evm7s2340","Nk52_evm14s2596","Nk52_evm88s270","Nk52_evm40s292","Nk52_evm1s2634","Nk52_evm50s236","Nk52_evm8s2462","Nk52_evm55s226","Nk52_evm6s217")
Cperk_develOGs<-na.omit(Cperk_OGs2[match(Cperk_develGenes,Cperk_OGs2[,2]),1])
GOogs<-Cperk_develOGs

###GO:0006096 -- glycolysis
GOgenes<-ls(attributes(attributes(attributes(GOdataBP)$graph)$nodeData)$data$`GO:0006096`$genes)
GOogs<-sample(GOgenes,34)


matrices <- OG_matList[[1]]  # Extract the first list of matrices
###lets make sure to look at all patterns per OG!
base_ids <- ifelse(grepl("\\.", GOogs), sub("\\..*", "", GOogs), GOogs)
patterns <- paste0("^", base_ids, "(\\.\\d+)?$")
OGpats <- unlist(lapply(patterns, function(pat) grep(pat, rownames(matrices[[1]]), value = TRUE)))

if(length(OGpats)<2){next}

# Subset all matrices by the selected rownames
basals <- lapply(matrices, function(x) x[OGpats, ])


# Initialize a vector to store the average correlations
correlations <- c()

####compute pairwise correlations across every row for Cperk vs. each animal in "basals":
nrow(basals[[1]])
# Loop over each row
for (i in 1:nrow(basals[[1]])) {
  #print(i)
  # Compute pairwise correlations
  for (j in 2:11) {
      row1 <- basals[[1]][i, ]
      row2 <- basals[[j]][i, ]
	  if(all(is.na(row1)) || all(is.na(row2))){
	    #print(row1)
		#print(row2)
		correlations <- c(correlations, NA)
	  }else{
	  #print(row1)
	  #print(row2)
	  #print(cor(row1, row2,method="pearson"))
	  #print(diss.CDM(row1,row2))
	  if(cor(row1, row2,method="pearson")==0){
	  print(row1)
	  print(row2)
	  }
      correlations <- c(correlations, cor(row1, row2,method="pearson"))
	  }
  }
}




#for dev genes:
devcor<-median(correlations,na.rm=T)
#for glycolysis genes:
glycor<-median(correlations,na.rm=T)

##Figure S2j
hist(cormeds,breaks=100,main="",xlab="Correlation between C. perkinsii and Animal development")
abline(v=glycor,col="blue",lty=2,lwd=2)
abline(v=devcor,col="red",lty=2,lwd=3)


