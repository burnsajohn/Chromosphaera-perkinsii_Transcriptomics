###code to generates plots in Figure S2F--example network structures to demonstrate how network structure was used to select genes within an OG with similar expression patterns.

###first run:
###Cperk_TimeCourse_DEGanalyses.R
###Figure_2D.R
###to load and format relevant data structures

og_sample<-allOGs[360]
for(i in 1:length(matrices_list)){
matrices_common[[i]]<-matrices_list[[i]][matrices_list[[i]]$Orthogroup %in% og_sample, ]
}

mylines<-vector()
for(j in 1:length(matrices_common)){
	mylines<-rbind(mylines,matrices_common[[j]])#[grep(allogs[i],matrices_common[[j]]$Orthogroup),])
	
}


mylines<-vector()
for(j in 1:length(matrices_common)){
	mylines<-rbind(mylines,matrices_common[[j]])#[grep(allogs[i],matrices_common[[j]]$Orthogroup),])
	
}

OGgene<-cbind(rownames(mylines),mylines$Orthogroup)
OGedges<-matrix(nrow=length(mylines[,1]),ncol=length(mylines[,1]))
rownames(OGedges)<-rownames(mylines)
colnames(OGedges)<-rownames(mylines)
allOGs<-unique(OGgene[,2])
OGedges[]<-0L
diag(OGedges)<-1

for(i in 1:length(allOGs)){
#print(i)
ogns<-OGgene[grep(allOGs[i],OGgene[,2]),1]
#print(ogns)
if(length(ogns)<2){next}
ognpairs<-t(combn(ogns,2))
for(j in 1:length(ognpairs[,1])){
	OGedges[ognpairs[j,1],ognpairs[j,2]]<-1
	OGedges[ognpairs[j,2],ognpairs[j,1]]<-1
}
}

mythresh<-thresh
OGg <- graph.adjacency(OGedges)
mydists<-dist(mylines[,1:7])
mydists.mat<-as.matrix(mydists)
mydists.weight<-exp(-mydists.mat)
mydists.weight.plot<-mydists.weight[mydists.weight>0]
mydists.weight[is.infinite(mydists.weight)]<-0
mydists.weight[is.na(mydists.weight)]<-0
mydists.weight[mydists.weight<mythresh]<-0

#Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
g <- graph.adjacency(
  mydists.weight,
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)

# Create a named list where names are your vertex names and values are the indices
names_dict <- setNames(seq_along(V(g)), V(g)$name)

eogs <- as_edgelist(OGg)

edge_ids <- vector(length = nrow(eogs))

for(i in seq_len(nrow(eogs))) {
  print(i)
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
###Figure S2F-1: original network
plot(G2)

OGccs<-components(G2)

# Normalize the weights to [0, 1] range for gray scale
normalized_weights <- 1-(E(G2)$weight - min(E(G2)$weight)) / (max(E(G2)$weight) - min(E(G2)$weight))
# Convert the weights to gray levels
gray_levels <- gray(normalized_weights)
E(G2)$color <- gray_levels

###cluster to break up OG clusters that have different types of similarity--different patterns!
eb <- edge.betweenness.community(G2)

###Figure S2F-2: Clustered network
plot(eb,G2,vertex.size=4,vertex.label.cex=0.4)

# Get the membership vector
membership <- eb$membership

# Create a new graph with only the edges where both ends belong to the same community
E(G2)$color <- ifelse(membership[get.edges(G2, E(G2))[, 1]] == membership[get.edges(G2, E(G2))[, 2]], "black", "red")


# Create a new graph by removing the red edges
G3 <- delete_edges(G2, E(G2)[color == "red"])

#plot(G3)

whichdegree<-degreethresh
Isolated = which(igraph::degree(G3)<whichdegree)
G4 = delete.vertices(G3, Isolated)
Isolated = which(igraph::degree(G4)<1)
G4 = delete.vertices(G4, Isolated)

ccs<-components(G4)
bccs<-which(ccs$csize>=bccsthresh)


edge_widths <- rep(2, ecount(G4))  # Let's assume the default width is 1

normalized_weights <- 1-(E(G4)$weight - min(E(G4)$weight)) / (max(E(G4)$weight) - min(E(G4)$weight))
# Convert the weights to gray levels
gray_levels <- gray(E(G4)$weight)
E(G4)$color <- gray_levels

OGccs<-components(G4)

OG1n<-names(OGccs$membership[OGccs$membership==1])
OG2n<-names(OGccs$membership[OGccs$membership==2])

E(G4)[OG1n %--% OG1n]$color="hotpink2"
E(G4)[OG2n %--% OG2n]$color<-"goldenrod2"

LO = layout.fruchterman.reingold(G4) 

#Figure S2F-3 Network connected components
plot(G4, layout = LO,vertex.size=4,vertex.label.cex=0.4,edge.width = edge_widths)

#Figure S2F-4 Unique expression patterns per compnent.
ccst<-components(G4)
matplot(t(mylines[names(ccst$membership[ccst$membership==1]),1:7]),type="l",col="hotpink2",xlim=c(1,6))
matplot(t(mylines[names(ccst$membership[ccst$membership==2]),1:7]),type="l",col="goldenrod2",add=T)
