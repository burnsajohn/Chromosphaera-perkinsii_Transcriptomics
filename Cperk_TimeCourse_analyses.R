###libraries needed:

library(tximport)
library(edgeR)
library(gplots)
library(RColorBrewer)
library("e1071")
library(ggplot2)
library(reshape)
library(dplyr)
library(topGO)


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

###set working directory to location of quant directories from Salmon
setwd("quants")
list.files()

###sample table for analyses
sampletable <- read.table("Cperk_samples.txt", header = TRUE)
sampletable

###read in all of the quantification files
files <- file.path(".", sampletable$quant, "quant.sf")
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

plotMDS(y)

fit <- glmQLFit(y,design)

#[1] "conditionF"    "conditionT104" "conditionT120" "conditionT54" 
#[5] "conditionT72"  "conditionT84"  "conditionT90"  "conditionT96" 

###run statistical tests: All vs. T52 using contrasts:
qlf.54vs72 <- glmQLFTest(fit,contrast=c(0,0,-1,1,0,0,0))
summary(decideTests(qlf.54vs72))
plotMD(qlf.54vs72)
abline(h=c(-1, 1), col="blue")

thresh<-0.05
lfc <- -10000000000000000
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

timeDEs<-unique(c(top72ids,top84ids,top90ids,top96ids,top104ids,top120ids))

cpms <- edgeR::cpm(y, offset = y$offset, log = FALSE)
logcounts <- edgeR::cpm(y, offset = y$offset, log = TRUE)

###scale the data using differentially expressed genes only:
scaledata <- t(scale(t(logcounts[timeDEs,])))

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


#######facet plot all clusters#########################
###plot equivalent to Figure S2B in manuscript#########
#######################################################

i=1
#Subset the cores molten dataframe so we can plot the core
core2 <- Kmolten[Kmolten$cluster==i,]
core2$CLUSTER<-i

#get cluster 1
K2 <- (scaledata[names(kClusters[kClusters==i]),])
#calculate the correlation with the core
corscore <- function(x){cor(x,core2$value)}
score <- apply(K2, 1, corscore)
#get the data frame into long format for plotting
K2molten <- melt(K2)
colnames(K2molten) <- c('gene','sample','value')
#add the score
K2molten <- merge(K2molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K2molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K2molten$order_factor <- 1:length(K2molten$gene)
#order the dataframe by score
K2molten <- K2molten[order(K2molten$score),]
#set the order by setting the factors
K2molten$order_factor <- factor(K2molten$order_factor , levels = K2molten$order_factor)
K2molten$CLUSTER<-i
plotclusters<-K2molten
coreclusters<-core2

#get the remaining clusters
for(i in 2:length(unique(Kmolten$cluster))){
#Subset the cores molten dataframe so we can plot the core
core2 <- Kmolten[Kmolten$cluster==i,]
core2$CLUSTER<-i
#get cluster 2
K2 <- (scaledata[names(kClusters[kClusters==i]),])
#calculate the correlation with the core
corscore <- function(x){cor(x,core2$value)}
score <- apply(K2, 1, corscore)
#get the data frame into long format for plotting
K2molten <- melt(K2)
colnames(K2molten) <- c('gene','sample','value')
#add the score
K2molten <- merge(K2molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K2molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K2molten$order_factor <- 1:length(K2molten$gene)
#order the dataframe by score
K2molten <- K2molten[order(K2molten$score),]
#set the order by setting the factors
K2molten$order_factor <- factor(K2molten$order_factor , levels = K2molten$order_factor)
K2molten$CLUSTER<-i

plotclusters<-rbind(plotclusters,K2molten)
coreclusters<-rbind(coreclusters,core2)

}

###plot clusters as barplots over time:
plotclusters2<-plotclusters
plotclusters2$time <- as.numeric(sapply(strsplit(as.character(plotclusters2$sample), "_"), function(x) as.numeric(x[2])))
# Calculate medians
medians <- plotclusters2 %>% dplyr::group_by(CLUSTER, time) %>% dplyr::summarise(median_value = median(value))
# Merge medians back with original data for plotting
plot_data <- merge(plotclusters2, medians, by = c("CLUSTER", "time"))

# Create the plot
p2 <- ggplot(plot_data, aes(x = as.factor(time), y = value)) +
  geom_boxplot() +
  geom_line(data = medians, aes(x = as.factor(time), y = median_value, group = CLUSTER), color = "blue") +
  facet_grid(cols = vars(CLUSTER)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Print the plot Figure S2B
print(p2)
####################################################################

##############################################################################################################
###GO analyses of each cluster################################################################################
##############################################################################################################

###UID mapping based GO annotations:
geneID2GO <- readMappings(file = "Nk52_Cperk_UID.goa", IDsep=";|,")


#order clusters based on expression profiles:

#Old Cluster 5 becomes Cluster 1
#Old Cluster 2 remains Cluster 2
#Old Cluster 1 becomes Cluster 3
#Old Cluster 3 becomes Cluster 4
#Old Cluster 4 becomes Cluster 5

clust1<-names(kClusters[kClusters==5])
clust2<-names(kClusters[kClusters==2])
clust3<-names(kClusters[kClusters==1])
clust4<-names(kClusters[kClusters==3])
clust5<-names(kClusters[kClusters==4])

myInterestingGenes <-list(clust1, clust2, clust3, clust4, clust5)

###Complete GO enrichment analyses for each cluster:
for(i in 1:length(myInterestingGenes)){
pat<-i

geneNames<-names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes[[pat]]))
names(geneList) <- geneNames
sum(as.numeric(as.character(geneList)))
length(myInterestingGenes[[pat]])

GOdataBP <- tryCatch(new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO), error=function(e) "error")
resultFisherBP.weight <- runTest(GOdataBP, algorithm = "weight", statistic = "fisher") 
allResBP<- GenTable(GOdataBP, weight = resultFisherBP.weight, orderBy = "weight", topNodes = 200, numChar=1000)
allResWeighted<- allResBP[as.numeric(gsub("<","",allResBP[,"weight"]))<0.05,]
#head(allResBP,10)
allResWeighted


gotblnm<-paste(pat,"_GOPBtbls.txt",sep="")

write.table(allResWeighted,gotblnm,quote=F,sep="\t")

library(plyr)
###function to pull genes associated with each GOterm.
getGOgenes<-function(allResTable,updown,GOdata){
	GOtermGenes<-list()
	for (i in 1:length(allResTable[,1])){
		cmd<-paste("ls(attributes(attributes(attributes(GOdata)$graph)$nodeData)$data$`",allResTable[i,1],"`$genes)",sep="")
		GOgenes<-eval(parse(text=cmd))
		updowngenes<-intersect(GOgenes,updown)
		GOtermGenes[[i]]<-updowngenes
	}
	names(GOtermGenes)<-allResTable[,1]
	gogenesDF<-ldply(GOtermGenes,cbind)
	colnames(gogenesDF)<-c("GO.ID","UNIgene")
	return(gogenesDF)
}

allResupdown<-allResBP
#filter GOterm by pvalue --column6 is weighted fisher test
updownres<-allResupdown[as.numeric(gsub("<","",allResupdown[,6]))<0.05,]
###get genes associated with sig GO terms

goupdown<-getGOgenes(updownres,myInterestingGenes[[pat]],GOdataBP)

gognnm<-goupdown
for(i in 1:length(gognnm[,1])){
gognnm[i,3]<-updownres[grep(gognnm[i,1],updownres[,1]),2]
}

###save GO enrichment tables
gognfilenm<-paste(pat,"_GOdf.txt",sep="")
write.table(gognnm,gognfilenm,quote=F,sep="\t")
}

