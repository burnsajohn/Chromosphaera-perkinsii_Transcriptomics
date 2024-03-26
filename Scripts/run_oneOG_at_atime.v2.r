
#degreethresh is how many edges are needed to keep a node -- good number is 4, permissive number is 3
#bccs thresh is how many nodes need to be connected to keep a connected components -- good number is 4, permissive number is 3

getOGpat<-function(thisOG, thresh, degreethresh, bccsthresh){
mylines.OG<-mylines[grep(thisOG,mylines$Orthogroup),]
mylines.OG2<-mylines.OG[,1:Cperk_numpatterns]
nan_rows <- rowSums(is.nan(as.matrix(mylines.OG2))) == ncol(as.matrix(mylines.OG2))
# Subset the dataframe to exclude these rows
mylines.OG2 <- mylines.OG2[!nan_rows, ]
OGpatternMAP<-matrix(ncol=2)

#mylines.OG2[is.na(mylines.OG2)] <- 0
# Create a logical index for rows where all numeric columns are NaN
# Exclude the last column (which is a character vector) from this check
mydists<-dist(mylines.OG2)
mydists.mat<-as.matrix(mydists)
#mydists.mat[mydists.mat>=2]<-0
mydists.weight<-exp(-mydists.mat)
mydists.weight.plot<-mydists.weight[mydists.weight>0]
mydists.weight[is.infinite(mydists.weight)]<-0
mydists.weight[is.na(mydists.weight)]<-0
mydists.weight[mydists.weight < thresh]<-0

#Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
g <- graph.adjacency(
  mydists.weight,
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)

# Create a named list where names are your vertex names and values are the indices
names_dict <- setNames(seq_along(V(g)), V(g)$name)

eogs <- as_edgelist(g)

edge_ids <- vector(length = nrow(eogs))

for(i in seq_len(nrow(eogs))) {
  #print(i)
  v1 <- eogs[i, 1]
  v2 <- eogs[i, 2]

  vidx1 <- names_dict[[v1]]
  vidx2 <- names_dict[[v2]]

  eid <- get.edge.ids(g, c(vidx1, vidx2))
  edge_ids[i] <- eid
}


keep_eids<-edge_ids[edge_ids>0]
all_edge_ids <- E(g)

# Delete all edges except those in specific_edges
G2 <- delete_edges(g, all_edge_ids[-keep_eids])

#G2<-set_edge_attr(g, "OGs",0)
#G2<-set_edge_attr(g, "OGs", index = edge_ids[edge_ids>0],1)
#E(G2)$color <- "grey20"
#E(G2)[edge_ids[edge_ids>0]]$color <- "red"

Isolated = which(igraph::degree(G2)<degreethresh)
G3 = delete.vertices(G2, Isolated)
LO = layout_with_fr(G3)
#LO = layout_in_circle(G3) 
#plot(G3, layout = LO,vertex.size=4,vertex.label.cex=0.4)


###cluster to break up OG clusters that have different types of similarity--different patterns!
eb <- edge.betweenness.community(G3)
#plot(eb,G3, layout = LO,vertex.size=4,vertex.label.cex=0.4)

# Get the membership vector
membership <- eb$membership

# Create a new graph with only the edges where both ends belong to the same community
E(G3)$color <- ifelse(membership[get.edges(G3, E(G3))[, 1]] == membership[get.edges(G3, E(G3))[, 2]], "black", "red")

# Create a new graph by removing the red edges
G4 <- delete_edges(G3, E(G3)[color == "red"])

ccs<-components(G4)
bccs<-which(ccs$csize>=bccsthresh)

# Helper function--deal with duplicates--more than one similarity cluster for each OG.
update_string <- function(thisOG, OGsUsed) {
  # If any element in thisOG is already in OGsUsed
  if (any(thisOG %in% OGsUsed)) {
    # Loop over each element in thisOG
    for (i in seq_along(thisOG)) {
      # If this element is in OGsUsed
      if (thisOG[i] %in% OGsUsed) {
        # Find all occurrences of thisOG in OGsUsed
        duplicates <- grepl(paste0("^", thisOG[i], "(\\.\\d+)?$"), OGsUsed)
        # Calculate the next duplicate number
        next_num <- sum(duplicates) + 1
        # Append the duplicate number to this element of thisOG
        thisOG[i] <- paste0(thisOG[i], ".", next_num)
      }
    }
  }
  return(thisOG)
}



###collect genes that fall into a developmental cluster for thisOG
devclust<-vector()
for(i in 1:length(bccs)){
devclust<-c(devclust,names(ccs$membership[which(ccs$membership==bccs[i])]))
}
###all genes matching OG, not just ones in the cluster.
clust2<-OGgene[grep(thisOG, OGgene[,2]), 1]
###remove genes found in one of the developmental clusters
clust3<-setdiff(clust2,devclust)
clust3map<-cbind(clust3,rep(thisOG,length(clust3)))
OGpatternMAP<-rbind(OGpatternMAP,clust3map)

###consolodate patterns for each developmental cluster in each animal
OG_matList<-vector(mode="list",length=length(matrices_list))
OGsUsed<-vector()
for(i in 1:length(bccs)){
clust<-names(ccs$membership[which(ccs$membership==bccs[i])])
#thisOG<-unique(OGgene[match(clust,OGgene[,1]),2])
thisOG.name <- update_string(thisOG, OGsUsed)
#print(thisOG.name)
OGsUsed<-c(OGsUsed,thisOG)
pattern_mat_row<-colMeans(mylines.OG2[clust,])
#pattern_mat_row[length(pattern_mat_row)+1]<-thisOG.name
# Split and extract first element
first_elements <- sapply(strsplit(clust, "_"), `[`, 1)
# Count unique strings
unique_count <- length(unique(first_elements))
# Check if "Cperk" is in the unique strings
contains_cperk <- as.integer("Cperk" %in% first_elements)
contains_sarc <- as.integer("Sarc" %in% first_elements)
# Output vector
pattern_mat_row[length(pattern_mat_row)+1]<-unique_count
pattern_mat_row[length(pattern_mat_row)+1]<-contains_cperk
pattern_mat_row[length(pattern_mat_row)+1]<-contains_sarc
pattern_mat_row<-t(as.data.frame(pattern_mat_row))
rownames(pattern_mat_row)<-thisOG.name
pattern_mat<-rbind(pattern_mat,pattern_mat_row)

clustmap<-cbind(clust,rep(thisOG.name,length(clust)))
OGpatternMAP<-rbind(OGpatternMAP,clustmap)


for(j in 1:length(matrices_list2)){
# Add new row
###if this organism does not have a gene in this cluster, take the overall average for that OG/organism... maybe just fill with NAs instead.
patmap<-vector()
if(length(na.omit(matrices_list2[[j]][clust,])[,1:Cperk_numpatterns][,1])==0){
#clust2--averages all genes for an OG in an organism
#clust3--averages all gense for an OG in an organism that don't match any bccs patterns (shared developmental patterns)
if(length(clust3)>0){
OG_matList[[j]] <- rbind(OG_matList[[j]], colMeans(na.omit(matrices_list2[[j]][clust3,])[,1:Cperk_numpatterns]))
}else{OG_matList[[j]] <- rbind(OG_matList[[j]], colMeans(na.omit(matrices_list2[[j]][clust2,])[,1:Cperk_numpatterns]))}
###if not in similarity cluster, give all genes generic OGname
#OG_matList[[j]] <- rbind(OG_matList[[j]], rep("NA",Cperk_numpatterns))
#print(c("toodiff",j,thisOG))
}else{
OG_matList[[j]] <- rbind(OG_matList[[j]], colMeans(na.omit(matrices_list2[[j]][clust,])[,1:Cperk_numpatterns]))
##if in similarity cluster, give all genes specific OGname associated with cluster.
#print(c("similar",j,thisOG))
}
# Assign rowname
rownames(OG_matList[[j]])[nrow(OG_matList[[j]])] <- thisOG.name
}
}

 return(list(OG_matList = OG_matList, pattern_mat = pattern_mat, OGpatternMAP= OGpatternMAP))

}







getOGpat_shuff<-function(thisOG, thresh){
mylines.OG<-mylines[grep(thisOG,mylines_shuff$Orthogroup),]
mylines.OG2<-mylines.OG[,1:Cperk_numpatterns]
nan_rows <- rowSums(is.nan(as.matrix(mylines.OG2))) == ncol(as.matrix(mylines.OG2))
# Subset the dataframe to exclude these rows
mylines.OG2 <- mylines.OG2[!nan_rows, ]

#mylines.OG2[is.na(mylines.OG2)] <- 0
# Create a logical index for rows where all numeric columns are NaN
# Exclude the last column (which is a character vector) from this check
mydists<-dist(mylines.OG2)
mydists.mat<-as.matrix(mydists)
#mydists.mat[mydists.mat>=2]<-0
mydists.weight<-exp(-mydists.mat)
mydists.weight.plot<-mydists.weight[mydists.weight>0]
mydists.weight[is.infinite(mydists.weight)]<-0
mydists.weight[is.na(mydists.weight)]<-0
mydists.weight[mydists.weight < thresh]<-0

#Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
g <- graph.adjacency(
  mydists.weight,
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)

# Create a named list where names are your vertex names and values are the indices
names_dict <- setNames(seq_along(V(g)), V(g)$name)

eogs <- as_edgelist(g)

edge_ids <- vector(length = nrow(eogs))

for(i in seq_len(nrow(eogs))) {
  #print(i)
  v1 <- eogs[i, 1]
  v2 <- eogs[i, 2]

  vidx1 <- names_dict[[v1]]
  vidx2 <- names_dict[[v2]]

  eid <- get.edge.ids(g, c(vidx1, vidx2))
  edge_ids[i] <- eid
}


keep_eids<-edge_ids[edge_ids>0]
all_edge_ids <- E(g)

# Delete all edges except those in specific_edges
G2 <- delete_edges(g, all_edge_ids[-keep_eids])

#G2<-set_edge_attr(g, "OGs",0)
#G2<-set_edge_attr(g, "OGs", index = edge_ids[edge_ids>0],1)
#E(G2)$color <- "grey20"
#E(G2)[edge_ids[edge_ids>0]]$color <- "red"

Isolated = which(igraph::degree(G2)<4)
G3 = delete.vertices(G2, Isolated)
LO = layout_with_fr(G3)
#LO = layout_in_circle(G3) 
#plot(G3, layout = LO,vertex.size=4,vertex.label.cex=0.4)


###cluster to break up OG clusters that have different types of similarity--different patterns!
eb <- edge.betweenness.community(G3)
#plot(eb,G3, layout = LO,vertex.size=4,vertex.label.cex=0.4)

# Get the membership vector
membership <- eb$membership

# Create a new graph with only the edges where both ends belong to the same community
E(G3)$color <- ifelse(membership[get.edges(G3, E(G3))[, 1]] == membership[get.edges(G3, E(G3))[, 2]], "black", "red")

# Create a new graph by removing the red edges
G4 <- delete_edges(G3, E(G3)[color == "red"])

ccs<-components(G4)
bccs<-which(ccs$csize>=4)

# Helper function--deal with duplicates--more than one similarity cluster for each OG.
update_string <- function(thisOG, OGsUsed) {
  # If any element in thisOG is already in OGsUsed
  if (any(thisOG %in% OGsUsed)) {
    # Loop over each element in thisOG
    for (i in seq_along(thisOG)) {
      # If this element is in OGsUsed
      if (thisOG[i] %in% OGsUsed) {
        # Find all occurrences of thisOG in OGsUsed
        duplicates <- grepl(paste0("^", thisOG[i], "(\\.\\d+)?$"), OGsUsed)
        # Calculate the next duplicate number
        next_num <- sum(duplicates) + 1
        # Append the duplicate number to this element of thisOG
        thisOG[i] <- paste0(thisOG[i], ".", next_num)
      }
    }
  }
  return(thisOG)
}


OG_matList<-vector(mode="list",length=length(matrices_list))
OGsUsed<-vector()
for(i in 1:length(bccs)){
clust<-names(ccs$membership[which(ccs$membership==bccs[i])])
thisOG<-unique(OGgene[match(clust,OGgene[,1]),2])
thisOG.name <- update_string(thisOG, OGsUsed)
#print(thisOG.name)
OGsUsed<-c(OGsUsed,thisOG)
for(j in 1:length(matrices_list)){
# Add new row
###if this organism does not have a gene in this cluster, take the overall average for that OG/organism... maybe just fill with NAs instead.
if(length(na.omit(matrices_list[[j]][clust,])[,1:Cperk_numpatterns][,1])==0){
#clust2<-OGgene[grep(paste(thisOG, collapse = "|"), OGgene[,2]), 1]
#OG_matList[[j]] <- rbind(OG_matList[[j]], colMeans(na.omit(matrices_list[[j]][clust2,])[,1:Cperk_numpatterns]))
OG_matList[[j]] <- rbind(OG_matList[[j]], rep("NA",Cperk_numpatterns))
#print(c("toodiff",j,thisOG))
}else{
OG_matList[[j]] <- rbind(OG_matList[[j]], colMeans(na.omit(matrices_list[[j]][clust,])[,1:Cperk_numpatterns]))
#print(c("similar",j,thisOG))
}
# Assign rowname
rownames(OG_matList[[j]])[nrow(OG_matList[[j]])] <- thisOG.name
}
}

return(OG_matList)

}
