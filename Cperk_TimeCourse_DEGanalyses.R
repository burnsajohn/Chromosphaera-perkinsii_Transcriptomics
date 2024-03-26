###libraries needed:
library(tximport)
library(edgeR)
library(gplots)
library(RColorBrewer)
library("e1071")
library(ggplot2)
library(reshape)
library(dplyr)


###files needed#########################
#sample control file: "Cperk_samples.txt"
#quant directories from Salmon output
#species specific gene ontology annotation file: "Nk52_Cperk_UID.goa"
########################################

#set some plot color parameters
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[sampleTable[,1]]

###sample table for analyses
sampletable <- read.table("Cperk_samples.txt", header = TRUE)
sampletable

###read in all of the quantification files
files <- file.path("quants", sampletable$quant, "quant.sf")
filenames<-matrix(unlist(strsplit(files,"/")),ncol=24)[2,]

namearr<-vector(length=length(filenames))
for(i in 1:length(filenames)){
namearr[i]<-sampletable[match(filenames[i],sampletable[,3]),2]
}

names(files) <- namearr
txi.salmon <- tximport(files, type = "salmon", txOut=T)
head(txi.salmon$counts)

cts <- txi.salmon$counts
normMat <- txi.salmon$length

###collect only A and C replicates due to evidence of infection in B reps.
cts<-cts[,c(1,3,4,6,7,9,10,12,13,15,16,18,19,21,22,24)]
normMat<-normMat[,c(1,3,4,6,7,9,10,12,13,15,16,18,19,21,22,24)]

cts<-cts[,3:16]
normMat<-normMat[,3:16]
##########################################################

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
y <- DGEList(cts)
y <- scaleOffset(y, normMat)

# filtering using the design information
sampleTable <- data.frame(condition = factor(rep(c("T54", "T72", "T84", "T90", "T96", "T104", "T120"), each = 2)))

design <- model.matrix(~0+condition, data = sampleTable)
keep <- filterByExpr(y, design)
y <- y[keep, ]

# y is now ready for estimate dispersion functions see edgeR User's Guide
y <- estimateDisp(y,design)

###multiple dimensional scaling plot to visualize replicates. 
plotMDS(y)

#Fit a quasi-likelihood negative binomial generalized log-linear model to the data
fit <- glmQLFit(y,design)

#[1] "conditionF"    "conditionT104" "conditionT120" "conditionT54" 
#[5] "conditionT72"  "conditionT84"  "conditionT90"  "conditionT96" 

###run statistical tests: All vs. T52 using contrasts:
#########################################
#parameters in case you want to set FDR or log fold change thresholds. For the manuscript we had no log fold change threshold and used an FDR of 0.05
thresh<-0.05
lfc <- -10000000000000000
########################################

qlf.54vs72 <- glmQLFTest(fit,contrast=c(0,0,-1,1,0,0,0))
summary(decideTests(qlf.54vs72))
plotMD(qlf.54vs72)
abline(h=c(-1, 1), col="blue")
top72 <- topTags(qlf.54vs72, n=Inf)
top72ids <- rownames(top72$table[top72$table$FDR<thresh & top72$table$logFC>=lfc,])

qlf.54vs84 <- glmQLFTest(fit,contrast=c(0,0,-1,0,1,0,0))
summary(decideTests(qlf.54vs84))
plotMD(qlf.54vs84)
abline(h=c(-1, 1), col="blue")
top84 <- topTags(qlf.54vs84, n=Inf)
top84ids <- rownames(top84$table[top84$table$FDR<thresh & top84$table$logFC>=lfc,])

qlf.54vs90 <- glmQLFTest(fit,contrast=c(0,0,-1,0,0,1,0))
summary(decideTests(qlf.54vs90))
plotMD(qlf.54vs90)
abline(h=c(-1, 1), col="blue")
top90 <- topTags(qlf.54vs90, n=Inf)
top90ids <- rownames(top90$table[top90$table$FDR<thresh & top90$table$logFC>=lfc,])

qlf.54vs96 <- glmQLFTest(fit,contrast=c(0,0,-1,0,0,0,1))
summary(decideTests(qlf.54vs96))
plotMD(qlf.54vs96)
abline(h=c(-1, 1), col="blue")
top96 <- topTags(qlf.54vs96, n=Inf)
top96ids <- rownames(top96$table[top96$table$FDR<thresh & top96$table$logFC>=lfc,])

qlf.54vs104 <- glmQLFTest(fit,contrast=c(1,0,-1,0,0,0,0))
summary(decideTests(qlf.54vs104))
plotMD(qlf.54vs104)
abline(h=c(-1, 1), col="blue")
top104 <- topTags(qlf.54vs104, n=Inf)
top104ids <- rownames(top104$table[top104$table$FDR<thresh & top104$table$logFC>=lfc,])

qlf.54vs120 <- glmQLFTest(fit,contrast=c(0,1,-1,0,0,0,0))
summary(decideTests(qlf.54vs120))
plotMD(qlf.54vs120)
abline(h=c(-1, 1), col="blue")
top120 <- topTags(qlf.54vs120, n=Inf)
top120ids <- rownames(top120$table[top120$table$FDR<thresh & top120$table$logFC>=lfc,])

###compile DEGs from all timepoints into a single list:
timeDEs<-unique(c(top72ids,top84ids,top90ids,top96ids,top104ids,top120ids))

###convert data to cpms and logcounts:
cpms <- edgeR::cpm(y, offset = y$offset, log = FALSE)
logcounts <- edgeR::cpm(y, offset = y$offset, log = TRUE)

###scale the data using differentially expressed genes only:
scaledata <- t(scale(t(logcounts[timeDEs,])))
#########################################################################################
####cluster data#########################################################################
###calculate clustering stats to determine cluster number: WSS=The sum distance within the centroids
wss <- (nrow(scaledata)-1)*sum(apply(scaledata,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(scaledata,
                                     centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

####Use k=5 based on wss plot:
set.seed(20)
kClust <- cmeans(as.matrix(scaledata), centers = 5, iter.max = 1000, m=2)
kClusters <- kClust$cluster

# function to find centroid in cluster i
clust.centroid = function(i, dat, clusters) {
  myname<- names(clusters[clusters == i])
  colMeans(dat[myname,])
}

kClustcentroids <- sapply(levels(factor(kClusters)), clust.centroid, scaledata, kClusters)

#get in long form for plotting
Kmolten <- melt(kClustcentroids)
colnames(Kmolten) <- c('sample','cluster','value')

#plot centroids:
p1 <- ggplot(Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() + 
  geom_line() +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster Expression by Time",color = "Cluster")
p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

####DEG and cluster data is used for further analyses and visualization of the data.
