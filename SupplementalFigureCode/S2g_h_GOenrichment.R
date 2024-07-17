###This script will run GO enrichment analysis on the output of the correlation analysis of figure S2G or S2H. 

#libraries needed:
library(topGO)
library(plyr)

#first source the data:

#for S2G--positively correlated genes in early branching animals
source("Figure_S2G.R")

#for S2H--positively correlated genes across all animals
#source("Figure_S2H.R")

####################GOenrichment#####################################

###UID mapping based GO annotations:
geneID2GO <- readMappings(file = "Annotations/OGcat.UID.goa", IDsep=";|,")
geneID2GO<-geneID2GO[grep("^OG",names(geneID2GO))]

myInterestingGenes <-Cperk_highCorGenes[,1]

geneNames<-names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
sum(as.numeric(as.character(geneList)))
length(myInterestingGenes)

GOdataBP <- tryCatch(new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO), error=function(e) "error")
resultFisherBP.weight <- runTest(GOdataBP, algorithm = "weight", statistic = "fisher") 
allResBP<- GenTable(GOdataBP, weight = resultFisherBP.weight, orderBy = "weight", topNodes = 200, numChar=1000)
allResWeighted<- allResBP[as.numeric(gsub("<","",allResBP[,"weight"]))<0.05,]
#head(allResBP,10)
head(allResWeighted)

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

goupdown<-getGOgenes(updownres,myInterestingGenes,GOdataBP)

gognnm<-goupdown
for(i in 1:length(gognnm[,1])){
gognnm[i,3]<-updownres[grep(gognnm[i,1],updownres[,1]),2]
}

##If the analysis was done on S2G--positively correlated genes in early branching animals, the output will be relevant to:
##Extended Data Table 2--Tab: Fig S2G
#write.table(allResWeighted,"GOenrichment_HighCor_basalAnimals.txt",sep="\t",quote=F)


##If the analysis was done on S2H--positively correlated genes across all animals, the output will be relevant to:
##Extended Data Table 2--Tab: Fig S2H
#write.table(allResWeighted,"GOenrichment_HighCor_allAnimals.txt",sep="\t",quote=F)



