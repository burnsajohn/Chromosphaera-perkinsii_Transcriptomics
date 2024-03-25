###This code will the plots for Figure 2C from the manuscript. 
###It uses count data for C. perkinsii generated for this project and incoporates count data generated for Sphaeroforma arctica previously: https://elifesciences.org/articles/49801
###Figure 3—source data 2: https://cdn.elifesciences.org/articles/49801/elife-49801-fig3-data2-v4.txt


#libraries needed:
library(dplyr)
library(tibble)
library(zoo)
library(matrixStats)
library(scales)
library(tximport)
library(edgeR)
library(matrixStats)
library(igraph)
library(topGO)

###required code collections:
source("Scripts/run_oneOG_at_atime.v2.r")
source("Scripts/animal_analysis2.rscript")

####read in the C. perkinsii data:
sampletable <- read.table("Cperk_samples.txt", header = TRUE)
sampletable

files <- file.path("quants", sampletable$quant, "quant.sf")
filenames<-matrix(unlist(strsplit(files,"/")),ncol=21)[2,]

namearr<-vector(length=length(filenames))
for(i in 1:length(filenames)){
namearr[i]<-sampletable[match(filenames[i],sampletable[,3]),2]
}

names(files) <- namearr
txi.salmon <- tximport(files, type = "salmon", txOut=T)
head(txi.salmon$abundance)

###use TPMs to be compatible with 2012 paper.
Cperk.tpm<-txi.salmon$abundance

###remove "B" replicates
Cperk.tpm<-Cperk.tpm[,c(1,3,4,6,7,9,10,12,13,15,16,18,19,21)]

minsum<-min(Cperk.tpm[Cperk.tpm>0])/2
Cperk_logcounts<-log(Cperk.tpm+minsum, base=2)

###open map of C.perkinsii gene names to orthogroups defined across animals:
Cperk_OGs<-read.table("Orthogroups/Nk52Cperk_TRX_OGdf.txt",header=T)
Cperk_OGs2<-Cperk_OGs

#convert logcount data to data frame and split out A and C replicates
Cperk_OGplot<-as.data.frame(Cperk_logcounts)
CperkA_OGplot<-Cperk_OGplot[,c(1,3,5,7,9,11,13)]
CperkC_OGplot<-Cperk_OGplot[,c(2,4,6,8,10,12,14)]

###format the dataframes for mapping OGs to count data
CperkA_OGplot <- CperkA_OGplot %>% rownames_to_column("Nk52Cperk_TRX")
CperkA_OGplot <- CperkA_OGplot %>% left_join(Cperk_OGs2, by = "Nk52Cperk_TRX")
CperkC_OGplot <- CperkC_OGplot %>% rownames_to_column("Nk52Cperk_TRX")
CperkC_OGplot <- CperkC_OGplot %>% left_join(Cperk_OGs2, by = "Nk52Cperk_TRX")

rownames(CperkA_OGplot)<-CperkA_OGplot$Nk52Cperk_TRX
CperkA_OGplot_OG<-CperkA_OGplot[, !(colnames(CperkA_OGplot) %in% c("Nk52Cperk_TRX", "vars"))]
rownames(CperkC_OGplot)<-CperkC_OGplot$Nk52Cperk_TRX
CperkC_OGplot_OG<-CperkC_OGplot[, !(colnames(CperkC_OGplot) %in% c("Nk52Cperk_TRX", "vars"))]
   
CperkA_OGplot<-t(scale(t(CperkA_OGplot_OG[, !(colnames(CperkA_OGplot_OG) %in% c("Orthogroup"))])))
CperkC_OGplot<-t(scale(t(CperkC_OGplot_OG[, !(colnames(CperkC_OGplot_OG) %in% c("Orthogroup"))])))

CperkA_OGplot<-as.data.frame(CperkA_OGplot)
CperkC_OGplot<-as.data.frame(CperkC_OGplot)

colnames(CperkA_OGplot)<-c("T54","T72","T84","T90","T96","T104","T120")
colnames(CperkC_OGplot)<-c("T54","T72","T84","T90","T96","T104","T120")

CperkA_OGplot$Orthogroup<-CperkA_OGplot_OG$Orthogroup
CperkC_OGplot$Orthogroup<-CperkC_OGplot_OG$Orthogroup

###########################################################################################
#####READ in the S. arctica data##########################################################
############################################################################################

Sarc.tpm<-read.table("elife-49801-fig3-data2-v4.txt",header=T)

row_values <- Sarc.tpm[,"gene"]
# Remove the row from the data frame
Sarc.tpm <-  Sarc.tpm[, -which(names(Sarc.tpm) == "gene")]
# Set the extracted values as the new rowname
rownames(Sarc.tpm) <- row_values

minsum<-min(Sarc.tpm[Sarc.tpm>0])/2
Sarc_logcounts<-log(Sarc.tpm+minsum, base=2)

#> colnames(Sarc_logcounts)
# [1] "BT12" "BT18" "BT24" "BT30" "BT36" "BT42" "BT48" "BT54" "BT60" "BT66"
#[11] "CT12" "CT18" "CT24" "CT30" "CT36" "CT42" "CT54" "CT60" "CT66" "CT72"

#For comparison of the two developmental time courses, select replicates as follows, based on how data was used in Dudin et. al 2019:
#Rep 1 : 18, 30, 42, 48, 54, 60, 66
#Rep2 : 18, 30, 42, 54, 60 ,66 ,72
#With these timepoints you keep all the key timepoints covering the whole growth type etc.
SarcB<-Sarc_logcounts[,c(2,4,6,7,8,9,10)]
SarcC<-Sarc_logcounts[,c(12,14,16,17,18,19,20)]
colnames(SarcB)<-c("T18","T30","T42","T48","T54","T60","T66")
colnames(SarcC)<-c("T18","T30","T42","T54","T60","T66","T72")

####read in the S. arctica gene to OG map:
Sarc_OGs<-read.table("Orthogroups/Sarc4_TRX_OGdf.txt",header=T)

###Here we are collecting common orthogroups between S. arctical and C. perkinsii:
Sarc_OGs_timeDEs <- Sarc_OGs[Sarc_OGs[,1] %in% Cperk_OGs2[,1], ]
SarcCons <- unique(Sarc_OGs_timeDEs[,2])
  
#get expressed genes (some might not be even if they have an ortholog!)
SarcConsPres <- rownames(Sarc_logcounts)[rownames(Sarc_logcounts) %in% SarcCons]
SarcBConsPres <- rownames(SarcB)[rownames(SarcB) %in% SarcCons]
SarcCConsPres <- rownames(SarcC)[rownames(SarcC) %in% SarcCons]

mygenes <- intersect(SarcBConsPres,SarcCConsPres)

Sarc_OGs2 <- Sarc_OGs[Sarc_OGs[,2] %in% mygenes, ]
#because had to remove info from gene, and some genes had multiple proteins, end up with duplicate rows per gene. Remove here without consequence for expression analyses
Sarc_OGs2 <- Sarc_OGs2[!duplicated(Sarc_OGs2$Sarc4_TRX), ]
colnames(Sarc_OGs2)[2]<-"TRX_name"
  
SarcB_OGplot <- as.data.frame(SarcB[mygenes,])
SarcC_OGplot <- as.data.frame(SarcC[mygenes,])

ntimes="all"
if(ntimes=="all"){ntimes=length(mygenes)}

Cperk_numpatterns<-length(Cperk_OGplot[1,])-1
 
plotMDS(SarcB_OGplot)
plotMDS(SarcC_OGplot)

SarcB_OGplot_scale <- t(scale(t(SarcB_OGplot)))
SarcC_OGplot_scale <- t(scale(t(SarcC_OGplot)))

###add column of OG names to data frame
SarcB_OGplot_scale <- as.data.frame(SarcB_OGplot_scale) %>% rownames_to_column("TRX_name")
SarcC_OGplot_scale <- as.data.frame(SarcC_OGplot_scale) %>% rownames_to_column("TRX_name")

SarcB_OGplot_scale <- SarcB_OGplot_scale %>% left_join(Sarc_OGs2, by = "TRX_name")
SarcC_OGplot_scale <- SarcC_OGplot_scale %>% left_join(Sarc_OGs2, by = "TRX_name")
  
###rename rows to gene name, remove unneeded columns:
rownames(SarcB_OGplot_scale)<-SarcB_OGplot_scale$TRX_name
genecol<-"Sarc4_TRX"
SarcB_OGplot_scale<-SarcB_OGplot_scale[, !(colnames(SarcB_OGplot_scale) %in% c("TRX_name", genecol, "vars"))]

rownames(SarcC_OGplot_scale)<-SarcC_OGplot_scale$TRX_name
genecol<-"Sarc4_TRX"
SarcC_OGplot_scale<-SarcC_OGplot_scale[, !(colnames(SarcC_OGplot_scale) %in% c("TRX_name", genecol, "vars"))]

SarcB_OGplot<-SarcB_OGplot_scale
SarcC_OGplot<-SarcC_OGplot_scale

######Now we have all data loaded and orthogroups assigned for both species####################
###need common colnames for downstream steps
colnames(CperkA_OGplot)<-c("T1","T2","T3","T4","T5","T6","T7","Orthogroup")
matrices_list<-list(CperkA_OGplot,CperkC_OGplot,SarcB_OGplot,SarcC_OGplot)
matrices_list<-rename_columns(matrices_list)

###get all orthogroups:
allogs<-vector()
for(i in 1:length(matrices_list)){
	allogs<-c(allogs,matrices_list[[i]]$Orthogroup)
}
allogs<-unique(allogs)
allogs<-allogs[order(allogs)]

# Count the occurrences of each Orthogroup across all matrices
all_orthogroups <- unlist(lapply(matrices_list, function(dt) dt$Orthogroup))
orthogroups_table <- table(all_orthogroups)

# Set the tolerance - the number of matrices that can be missing the orthogroup
tolerance <- 4 #for just these two species, looking at individual replicates a tolerance of 4 will retain all orthogroups for downstream analyses. This code is not strictly necessary here, but is used in other analyses to set a threshold when additional species are being compared.

# If you have N matrices/lists, you want strings present in at least (N - tolerance)
required_presence <- length(matrices_list) - tolerance

# Filter orthogroups based on the threshold
common_orthogroups <- names(orthogroups_table[orthogroups_table >= required_presence])
length(common_orthogroups)

matrices_common<-list()
for(i in 1:length(matrices_list)){
matrices_common[[i]]<-matrices_list[[i]][matrices_list[[i]]$Orthogroup %in% common_orthogroups, ]
}

names(matrices_common)<-c("CperkA","CperkC","SarcB","SarcC")
names(matrices_list)<-c("CperkA","CperkC","SarcB","SarcC")

# Function to append list element name to row names
append_list_name_to_rownames <- function(matrix, list_element_name) {
  rownames(matrix) <- paste(list_element_name, rownames(matrix), sep = "_")
  return(matrix)
}

# Apply the function to each matrix in the list
matrices_common2 <- lapply(names(matrices_common), function(name) {
  append_list_name_to_rownames(matrices_common[[name]], name)
})


# Apply the function to each matrix in the list
matrices_list2 <- lapply(names(matrices_list), function(name) {
  append_list_name_to_rownames(matrices_list[[name]], name)
})

mylines<-vector()
for(j in 1:length(matrices_common2)){
	print(j)
	mylines<-rbind(mylines,matrices_common2[[j]])#[grep(allogs[i],matrices_common[[j]]$Orthogroup),])
	
}

OGgene<-cbind(rownames(mylines),mylines$Orthogroup)
allOGs<-unique(OGgene[,2])

allOGs<-allOGs[order(allOGs)]
length(allOGs)


###use a network structure to split orthogroups by expression pattern: Each orthogroup might contain several distinct expression patterns across developmental time the following code breaks those patterns up into distinct orthogroup patterns. So: OG1 can become OG1.1, OG1.2, OG1.3 if it has enough genes that have distinct expression patterns.

#initialize OG_matList 
Cperk_numpatterns<-7
thresh=0.2
bccsthresh=2
degreethresh=2
pattern_mat<-matrix(nrow=0,ncol=(Cperk_numpatterns+3))
colnames(pattern_mat)<-c(colnames(mylines)[1:Cperk_numpatterns],"numSpec","CperkPA","SarcPA")
OG_matList<-getOGpat(allOGs[2],thresh,degreethresh,bccsthresh)

#length(allOGs)
# Loop through OGs
for (i in 3:length(allOGs)) {
  print(i)
  thisOGmatList<-getOGpat(allOGs[i],thresh,degreethresh,bccsthresh)
  # Loop through each matrix in the list
  if(length(rownames(thisOGmatList[[1]][[1]]))>0){
   for (j in 1:length(OG_matList[[1]])) {
     # Add new rows to the existing matrix using rbind
     OG_matList[[1]][[j]] <- rbind(OG_matList[[1]][[j]], thisOGmatList[[1]][[j]])
   }
   OG_matList[[2]]<-rbind(OG_matList[[2]], thisOGmatList[[2]])
   OG_matList[[3]]<-rbind(OG_matList[[3]], thisOGmatList[[3]])
  }   
}

for(i in 1:length(OG_matList[[1]])){
class(OG_matList[[1]][[i]])<-"numeric"
}

names(OG_matList[[1]])<-c("CperkA","CperkC","SarcB","SarcC")

###################Initialize a GO structure for the orthogroups#####################################
###UID mapping based GO annotations:
geneID2GO <- readMappings(file = "Annotations/OGcat.UID.goa", IDsep=";|,")
geneID2GO<-geneID2GO[grep("^OG",names(geneID2GO))]
myInterestingGenes <- as.character(na.omit(sample(matrices_list[[1]]$Orthogroup,200)))
geneNames<-names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
sum(as.numeric(as.character(geneList)))
length(myInterestingGenes)
GOdataBP <- tryCatch(new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO), error=function(e) "error")


########################GOs for comparison########################################################################
###Cell Cycle--GO:0007049
GOgenes1<-ls(attributes(attributes(attributes(GOdataBP)$graph)$nodeData)$data$`GO:0007049`$genes)
###Cytokenisis--GO:0000910
GOgenes2<-ls(attributes(attributes(attributes(GOdataBP)$graph)$nodeData)$data$`GO:0000910`$genes)
###DNA Replication--GO:0006260
GOgenes3<-ls(attributes(attributes(attributes(GOdataBP)$graph)$nodeData)$data$`GO:0006260`$genes)
###Cell Polarity--GO:0007163
GOgenes4<-ls(attributes(attributes(attributes(GOdataBP)$graph)$nodeData)$data$`GO:0007163`$genes)
###Vacuolar transport--GO:0007034
GOgenes5<-ls(attributes(attributes(attributes(GOdataBP)$graph)$nodeData)$data$`GO:0007034`$genes)
###Cell adhesion--GO:0098609
GOgenes6<-ls(attributes(attributes(attributes(GOdataBP)$graph)$nodeData)$data$`GO:0098609`$genes)

GOgenes<-list(GOgenes1,GOgenes2,GOgenes3,GOgenes4,GOgenes5,GOgenes6)
###################################################################################################################

###Map OGs between species for each GO gene and collect plotting order based on clustering in C. perkinsii:

		     

OGorder_cperk<-list()
OGorder_animal<-list()
for(i in 1:length(GOgenes)){
GOogs<-GOgenes[[i]]
animal<-3
cperkplotA<-matrices_list[[1]][matrices_list[[1]]$Orthogroup %in% GOogs,]
cperkplotB<-matrices_list[[2]][matrices_list[[2]]$Orthogroup %in% GOogs,]

#ANIMAL genes
animalplotB<-matrices_list[[animal]][matrices_list[[animal]]$Orthogroup %in% GOogs,]
animalplotC<-matrices_list[[animal+1]][matrices_list[[animal+1]]$Orthogroup %in% GOogs,]

matches <- unique(unlist(lapply(GOogs, function(pattern) grep(pattern, rownames(OG_matList[[1]][[1]]), value = TRUE, fixed = FALSE))))
cperkmatches<-rownames(na.omit(as.matrix(OG_matList[[1]][[1]][matches,1:Cperk_numpatterns])))


##which ones are not present in the animal OGs?
notAnimalB<-setdiff(cperkmatches,rownames(na.omit(as.matrix(OG_matList[[1]][[animal]][cperkmatches,1:Cperk_numpatterns]))))
notAnimalmatB<-matrices_list[[animal]][matrices_list[[animal]]$Orthogroup %in% notAnimalB,]

notAnimalC<-setdiff(cperkmatches,rownames(na.omit(as.matrix(OG_matList[[1]][[animal+1]][cperkmatches,1:Cperk_numpatterns]))))
notAnimalmatC<-matrices_list[[animal+1]][matrices_list[[animal+1]]$Orthogroup %in% notAnimalC,]

library(dplyr)
# Calculate the mean for each Orthogroup
getAnimalOGsB <- notAnimalmatB %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::summarize(across(starts_with("T"), ~mean(.x, na.rm = TRUE)), .groups = 'drop')

getAnimalOGsC <- notAnimalmatC %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::summarize(across(starts_with("T"), ~mean(.x, na.rm = TRUE)), .groups = 'drop')

getAnimalOGsB<-as.data.frame(getAnimalOGsB)
getAnimalOGsC<-as.data.frame(getAnimalOGsC)

# Set the first column as row names
rownames(getAnimalOGsB) <- getAnimalOGsB[,1]
rownames(getAnimalOGsC) <- getAnimalOGsC[,1]

# Remove the first column from the data frame
getAnimalOGsB <- getAnimalOGsB[,-1]
getAnimalOGsC <- getAnimalOGsC[,-1]

AnimalmatB<-as.data.frame(na.omit(as.matrix(OG_matList[[1]][[animal]][cperkmatches,1:Cperk_numpatterns])))
AnimalmatC<-as.data.frame(na.omit(as.matrix(OG_matList[[1]][[animal+1]][cperkmatches,1:Cperk_numpatterns])))

plotcperkA<-na.omit(as.matrix(OG_matList[[1]][[1]][cperkmatches,1:Cperk_numpatterns]))
plotcperkC<-na.omit(as.matrix(OG_matList[[1]][[2]][cperkmatches,1:Cperk_numpatterns]))

if(exists("getAnimalOGsB")){
plotAnimalB<-as.matrix(rbind(AnimalmatB,getAnimalOGsB))
plotAnimalC<-as.matrix(rbind(AnimalmatC,getAnimalOGsC))
}else{
plotAnimalB<-AnimalmatB
plotAnimalC<-AnimalmatC
}

plottable<-intersect(rownames(plotAnimalC),intersect(rownames(plotcperkC),intersect(rownames(plotcperkA),rownames(plotAnimalB))))


plotcperkA<-plotcperkA[plottable,]
plotcperkA<-rescale(plotcperkA, to = c(-2, 2))   

plotcperkC<-plotcperkC[plottable,]
plotcperkC<-rescale(plotcperkC, to = c(-2, 2))   

plotAnimalB<-plotAnimalB[plottable,]
plotAnimalB<-rescale(plotAnimalB, to = c(-2, 2))   

plotAnimalC<-plotAnimalC[plottable,]
plotAnimalC<-rescale(plotAnimalC, to = c(-2, 2))   


# Perform clustering on the first data frame
hc <- hclust(dist(plotcperkA))  # Using transposed data for row clustering
# Get the order of the rows from the clustering
row_order <- hc$order
OGorder_cperk[[i]]<-rownames(plotcperkA)[row_order]

hc <- hclust(dist(plotAnimalB))  # Using transposed data for row clustering
# Get the order of the rows from the clustering
row_order <- hc$order
OGorder_animal[[i]]<-rownames(plotAnimalB)[row_order]
}

		     
