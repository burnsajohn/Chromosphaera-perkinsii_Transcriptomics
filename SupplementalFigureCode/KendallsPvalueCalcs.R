#Extended Data Table 3: Tab: Kendall's p-value Dev genes

###First you must run: 
#source("Cperk_TimeCourse_DEGanalyses.R")
#source("Figure_2D.R")

#Libraries Needed
library(irr)

Cperk_develGenes<-c("Nk52_evm15s1360","Nk52_evm78s1737","Nk52_evm1s1535","Nk52_evm57s1073","Nk52_evm44s1360","Nk52_evm60s152","Nk52_evm4s2657","Nk52_evm64s2039","Nk52_evm98s2192","Nk52_evm16s267","Nk52_evm7s490","Nk52_evm24s1967","Nk52_evm90s1737","Nk52_evm65s270","Nk52_evm36s621","Nk52_evm4s293","Nk52_evm20s2474","Nk52_evm62s1020","Nk52_evm18s1967","Nk52_evm48s153","Nk52_evm11s256","Nk52_evm45s1569","Nk52_evm29s2367","Nk52_evm50s352","Nk52_evm4s217","Nk52_evm19s307","Nk52_evm18s164","Nk52_evm7s2340","Nk52_evm14s2596","Nk52_evm88s270","Nk52_evm40s292","Nk52_evm1s2634","Nk52_evm50s236","Nk52_evm8s2462","Nk52_evm55s226","Nk52_evm6s217")
Cperk_develOGs<-na.omit(Cperk_OGs2[match(Cperk_develGenes,Cperk_OGs2[,2]),1])
GOgenes<-Cperk_OGs2[match(Cperk_develOGs,Cperk_OGs2[,1]),]
GOogs<-Cperk_develOGs

merged_df <- merge(GOgenes, Gene.x.OG.patternMAP, by.x = "Nk52Cperk_TRX", by.y = "clust3")  
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
mybreaks <- seq(mymin,mymax,by=(abs((mymin-mymax)/20))) #c(-2.5, -2, -1.5,-1,-0.5,0,0.5,1,1.5,2,2.5)

###data tables
genedata.l<-list()
for(i in 1:length(GOogs)){
print(GOogs[i])
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
print(gene_plot)

GObase<-ifelse(grepl("\\.", GOogs[i]), sub("\\..*", "", GOogs[i]), GOogs[i])
#mtitle<-expanded_df[match(GOogs[i],expanded_df$V2),"Animals"]
if(!is.na(OGnameMap[match(GObase,OGnameMap$OG),"NAME"])){
mtitle<-OGnameMap[match(GObase,OGnameMap$OG),"NAME"]
}else{mtitle<-expanded_df[match(GOogs[i],expanded_df$V2),"Animals"]}

genedata.l[[i]]<-gene_plot
names(genedata.l)[i]<-mtitle
}


kw_data<-matrix(nrow=length(genedata.l),ncol=2)
for(i in 1:length(genedata.l)){
kw<-kendall(t(genedata.l[[i]]))
kw_data[i,1]<-names(genedata.l)[i]
kw_data[i,2]<-kw$p.value
}

###will create data table used in Extended Data Table 3, tab Kendall.
write.table(kw_data,file="kendallW.develGenes.txt",sep="\t",quote=F)
