###This code will produce the raw plots used in Figure 2D from the manuscript. 

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
library(dendextend)
library(pheatmap)

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

#convert logcount data to data frame
Cperk_OGplot<-as.data.frame(Cperk_logcounts)

###average biological replicates to get a single expression value per timepoint
Cperk_OGplot_noreps<-average_columns_pairs(Cperk_OGplot)

###instead of replacing rownames, let's add a column that gives the OG name.
#Cperk_OGplot<-replace_rownames(Cperk_OGplot_noreps, Cperk_OGs3, "Nk52Cperk_TRX", "Orthogroup")

Cperk_OGplot_noreps <- Cperk_OGplot_noreps %>% rownames_to_column("Nk52Cperk_TRX")
Cperk_OGplot_noreps <- Cperk_OGplot_noreps %>% left_join(Cperk_OGs2, by = "Nk52Cperk_TRX")

rownames(Cperk_OGplot_noreps)<-Cperk_OGplot_noreps$Nk52Cperk_TRX
Cperk_OGplot_noreps<-Cperk_OGplot_noreps[, !(colnames(Cperk_OGplot_noreps) %in% c("Nk52Cperk_TRX", "vars"))]
  
Cperk_OGplot<-t(scale(t(Cperk_OGplot_noreps[, !(colnames(Cperk_OGplot_noreps) %in% c("Orthogroup"))])))
###rescale between -2,2
Cperk_OGplot <-rescale(Cperk_OGplot, to = c(-2, 2))

Cperk_OGplot<-as.data.frame(Cperk_OGplot)

###which columns/timepoints are we using here???
###leaving off the last timepoint because it corresponds to the multicellular stage breaking up which diverges from animal/multicellular development programs.
Cperk_OGplot<-Cperk_OGplot[,c(1:6)]


###rename columns to match across animals
colnames(Cperk_OGplot)<-paste0("Avg_", seq_len(length(Cperk_OGplot[1,])))
Cperk_OGplot$Orthogroup<-Cperk_OGplot_noreps$Orthogroup

##########################################################################
###read in data for animal developmental time courses from: https://www.nature.com/articles/nature16994
####Data is averaged across time using a sliding window to convert gene expression from real time to developmental time--mapping gene changes across developmental events rather than time points. Time averaging is implemented to match the number of timpoints collected for C. perkinsii using sliding window rules similar to those of Levin et. al 2016. Expression files of animal development formatted for these analyses are provided in directory "Animal_Data". They were downloaded from GEO database GSE70185. If you reuse, please cite the original Levin et al. 2016 paper. (https://doi.org/10.1038/nature16994).

###Blastula numbers
Cperk_numpatterns<-length(Cperk_OGplot[1,])-1
PT_OGplot <- animal_analysis("PT",14,Cperk_numpatterns)
FW_OGplot <- animal_analysis("FW",9,Cperk_numpatterns)
CE_OGplot <- animal_analysis("CE",20,Cperk_numpatterns)
HD_OGplot <- animal_analysis("HD",15,Cperk_numpatterns)
DM_OGplot <- animal_analysis("DM",9,Cperk_numpatterns)
SU_OGplot <- animal_analysis("SU",22,Cperk_numpatterns)
ZF_OGplot <- animal_analysis("ZF",6,Cperk_numpatterns)
NV_OGplot <- animal_analysis("NV",10,Cperk_numpatterns)
AQ_OGplot <- animal_analysis("AQ",11,Cperk_numpatterns)
ML_OGplot <- animal_analysis("ML",11,Cperk_numpatterns)

matrices_list <- list(Cperk_OGplot,PT_OGplot,FW_OGplot,CE_OGplot,HD_OGplot,DM_OGplot,SU_OGplot,ZF_OGplot,NV_OGplot,AQ_OGplot,ML_OGplot)
names(matrices_list)<-c("Cperk", "PT","FW","CE","HD","DM","SU","ZF","NV","AQ","ML")

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
tolerance <- 4 # All but one, adjust this value accordingly

# If you have N matrices/lists, you want strings present in at least (N - tolerance)
required_presence <- length(matrices_list) - tolerance

# Filter orthogroups based on the threshold
common_orthogroups <- names(orthogroups_table[orthogroups_table >= required_presence])
length(common_orthogroups)

matrices_common<-list()
for(i in 1:length(matrices_list)){
matrices_common[[i]]<-matrices_list[[i]][matrices_list[[i]]$Orthogroup %in% common_orthogroups, ]
}

names(matrices_common)<-c("Cperk", "PT","FW","CE","HD","DM","SU","ZF","NV","AQ","ML")
names(matrices_list)<-c("Cperk", "PT","FW","CE","HD","DM","SU","ZF","NV","AQ","ML")

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
thresh=0.2
bccsthresh=3
degreethresh=3
pattern_mat<-matrix(nrow=0,ncol=(Cperk_numpatterns+2))
colnames(pattern_mat)<-c(colnames(mylines)[1:Cperk_numpatterns],"numSpec","CperkPA")
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

names(OG_matList[[1]])<-c("Cperk", "PT","FW","CE","HD","DM","SU","ZF","NV","AQ","ML")

###############################################################
######Plot Developmental Gene Sets#############################
###############################################################

dgenes<-read.table("Annotations/DevelGenes_MAP.txt",header=T,sep="\t")

Gene.x.OG.patternMAP<-as.data.frame(na.omit(OG_matList$OGpatternMAP))
Gene.x.OG.patternMAP[,1] <- sub("^[^_]*_", "", Gene.x.OG.patternMAP[,1])

###Here can select which functional class of gene to plot: TF, Adhesion, or Signaling. They are read from file DevelGenes_MAP.txt. Other OGs can be selected from annotation files for comparison. A subset is shown in Figure 2D in the manuscript.
GOgenes<-dgenes[dgenes[,3]=="TF",]

GOogs<-Gene.x.OG.patternMAP[Gene.x.OG.patternMAP[,1] %in% GOgenes[,2],2]
GOogs<-unique(GOogs)

merged_df <- merge(GOgenes, Gene.x.OG.patternMAP, by.x = "Cperk", by.y = "clust3")  
# Extract the base identifier
merged_df$BaseID <- gsub("\\..*", "", merged_df$V2)
# Aggregate based on the BaseID
combined_df <- aggregate(. ~ BaseID, data = merged_df, FUN = function(x) paste(unique(x), collapse = ","))
merged_df<-combined_df

###lets make sure to look at all patterns per OG!
base_ids <- ifelse(grepl("\\.", GOogs), sub("\\..*", "", GOogs), GOogs)
patterns <- paste0("^", base_ids, "(\\.\\d+)?$")

# Find all matching patterns in Gene.x.OG.patternMAP to get all expression patterns for a given OG/function
matched_rows <- lapply(patterns, function(p) Gene.x.OG.patternMAP[grep(p, Gene.x.OG.patternMAP[,2]), , drop = FALSE]) %>% bind_rows()

# Initialize an empty dataframe
expanded_df <- data.frame(V2 = character(), Cperk = character(), Animals = character(), Category = character(), stringsAsFactors = FALSE)

# Loop through each row of merged_df
for (i in 1:nrow(merged_df)) {
    current_row <- merged_df[i, ]
    base_id <- current_row$BaseID
    
    # Find matching rows in matched_rows that start with "NK52_"
    matches <- matched_rows[grepl(base_id, matched_rows$V2) & grepl("^Nk52_", matched_rows$clust3), ]
    
 if (nrow(matches) > 0) {
        unique_v2s <- unique(matches$V2)
        # Iterate through each unique V2 string
        for (v2 in unique_v2s) {
            new_row <- current_row
            new_row$V2 <- v2
            new_row$Cperk <- matches$clust3[matches$V2 == v2][1]  # Take the first NK52_ corresponding to the V2
            suffix <- ifelse(grepl("\\.", v2), gsub(".*(\\..*)", "\\1", v2), "")
            new_row$Animals <- paste0(current_row$Animals, suffix)
            expanded_df <- rbind(expanded_df, new_row)
        }
    }
}

# Print the expanded dataframe
print(expanded_df)

#GOgenes %in% Gene.x.OG.patternMAP[,1] 
#GOgenes %in% rownames(Cperk_OGplot[,1])

####manual dendrogram--This dendrogram follows the tree in Leven et al. 2016.
# Manually create an hclust object
my_hclust <- list(
  merge = matrix(c(-1, -2, -3, -4, -5, -6, 2, -7, 4, 1, 5, 3, 6, -8, 7, -9, 8, -10, 9, -11), ncol = 2, byrow = TRUE),
  height = 1:10,
  order = 1:11,
  labels = c("PT","FW","DM","HD","ZF","SU","CE","NV","AQ","ML","Cperk"),
  method = "complete",
  call = match.call()
)
class(my_hclust) <- "hclust"

# Convert to dendrogram
my_dendrogram <- as.dendrogram(my_hclust)
my_dendrogram_rotated <- rotate(my_hclust,c("PT","FW","CE","HD","DM","SU","ZF","NV","AQ","ML","Cperk"))  
#convert back to hclust 
my_dendrogram_rotated<-as.hclust(my_dendrogram_rotated)

GOogs<-expanded_df$V2
###heatmaps
geneheatmap.l<-list()
for(i in 1:length(GOogs)){
target_rowname <- GOogs[i] # Replace with your actual row name
# Function to extract a row and rename it
extract_and_rename_row <- function(matrix, rowname, new_name) {
  extracted_row <- matrix[rowname, , drop = FALSE]
  rownames(extracted_row) <- new_name
  return(extracted_row)
}
# Apply the function to each matrix in the list
extracted_rows <- lapply(names(OG_matList[[1]]), function(name) {
  extract_and_rename_row(OG_matList[[1]][[name]], target_rowname, name)
})
gene_plot<-do.call(rbind,extracted_rows)

gene_plot<-rescale(gene_plot, to = c(-2, 2))

gene_plot[is.na(gene_plot)] <- -2

#orderRows<-c("PT","FW","CE","HD","DM","SU","ZF","NV","AQ","ML","Cperk","Sarc")
#orderRows<-c("PT","FW","CE","HD","DM","SU","ZF","NV","AQ","ML","Cperk")
gene_plot<-gene_plot[my_dendrogram_rotated$labels,]


color_ramp <- viridis::cividis(20)
mymin <- -2
mymax <- 2
mybreaks <- seq(mymin,mymax,by=(abs((mymin-mymax)/20))) #c(-2.5, -2, -1.5,-1,-0.5,0,0.5,1,1.5,2,2.5)


mtitle<-expanded_df[match(GOogs[i],expanded_df$V2),"Animals"]
geneheatmap.l[[i]]<-pheatmap(as.matrix(gene_plot),cluster_cols = FALSE,cluster_rows=my_dendrogram_rotated,scale="none",silent=TRUE,legend=F,main=mtitle,fontsize=5,treeheight_row = 5,breaks=mybreaks,color = color_ramp)$gtable
}

###plots OG across animal development as heatmap ordered by phylogeny--Fig 2D. 
n <- length(geneheatmap.l)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(geneheatmap.l, ncol=nCol))

###make one plot to print scale
#pheatmap(as.matrix(gene_plot),cluster_cols = FALSE,cluster_rows=my_dendrogram_rotated,scale="none",main=mtitle,fontsize=5,treeheight_row = 5,breaks=mybreaks,color = color_ramp,cellwidth=40, cellheight=20)


