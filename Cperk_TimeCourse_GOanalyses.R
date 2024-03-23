##first have to run Cperk_TimeCourse_DEGanalyses.R to get objects needed for GO analyses

#library needed:
library(topGO)


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

