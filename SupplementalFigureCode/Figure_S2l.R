#####plot highly correlated genes###########################################
###this will make a plot of individual genes within an enriched GO category, as defined by the set of highly correlated genes. Makes figures relevant to Figure S2L.

#libraries needed:
library(dendextend)
library(pheatmap)

#First source the data and GO analyses needed:
source("SupplementalFigureCode/Figure_S2g.R") #alternatively: source("Figure_S2H.R")
source("SupplementalFigureCode/S2g_h_GOenrichment.R")

whichGO<-1 ###plot expression value heatmaps for genes associated with this GO term (ranked from 1--lowest enrichment p-value to n--highest enrichment p-value).
sigGOs<-goupdown[grep(unique(goupdown[,1])[whichGO],goupdown[,1]),2]
Cperk_highCorGenes<-Cperk_OGs2[match(sigGOs,Cperk_OGs2[,1]),]

GOgenes<-Cperk_highCorGenes

names(GOgenes)<-c("Animals","Cperk")
GOgenes$Category<-rep("highCor",length(GOgenes[,1]))

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

####manual dendrogram
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
my_dendrogram_rotated <- rotate(my_dendrogram,c("PT","FW","CE","HD","DM","SU","ZF","NV","AQ","ML","Cperk"))  
#convert back to hclust 
my_dendrogram_rotated<-as.hclust(my_dendrogram_rotated)

GOogs<-expanded_df$V2

OGnameMap<-read.table("Annotations/OG_NAME_MAP.txt", sep="\t",header=T,quote="")

color_ramp <- viridis::cividis(20)
mymin <- -2
mymax <- 2
mybreaks <- seq(mymin,mymax,by=(abs((mymin-mymax)/20)))

###heatmaps--make a list of heatmap structures, then plot later using grids.
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
gene_plot <- rescale(gene_plot, to = c(-2, 2))   

gene_plot[is.na(gene_plot)] <- -2

#gene_plot<- t(apply(gene_plot,1,function(x){rescale(x,to=c(-2,2))}))

#orderRows<-c("PT","FW","CE","HD","DM","SU","ZF","NV","AQ","ML","Cperk","Sarc")
#orderRows<-c("PT","FW","CE","HD","DM","SU","ZF","NV","AQ","ML","Cperk")
gene_plot<-gene_plot[my_dendrogram_rotated$labels,]

GObase<-ifelse(grepl("\\.", GOogs[i]), sub("\\..*", "", GOogs[i]), GOogs[i])
#mtitle<-expanded_df[match(GOogs[i],expanded_df$V2),"Animals"]
if(!is.na(OGnameMap[match(GObase,OGnameMap$OG),"NAME"])){
mtitle<-OGnameMap[match(GObase,OGnameMap$OG),"NAME"]
}else{mtitle<-expanded_df[match(GOogs[i],expanded_df$V2),"Animals"]}

geneheatmap.l[[i]]<-pheatmap(as.matrix(gene_plot),cluster_cols = FALSE,cluster_rows=my_dendrogram_rotated,scale="none",silent=TRUE,legend=F,main=mtitle,fontsize=5,treeheight_row = 5,breaks=mybreaks,color = color_ramp)$gtable
}

n <- length(geneheatmap.l)
n

###Figure S2L
###if n is large (too many genes to plot), can subset geneheatmap.l and plot sets of genes (9 is a good number):
geneheatmap.l.plot<-geneheatmap.l[1:9]
n <- length(geneheatmap.l.plot)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(geneheatmap.l.plot, ncol=nCol))

