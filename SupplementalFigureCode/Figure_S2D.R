###First run the code for Figure 2C to load appropriate data into R. Then load additional libraries and run this code to make violin plots of significance values.
###This script uses a correlation test (cor.test()) to test for significant correlation between expression patterns across two matrices.

#Libraries needed
library(eulerr)
library("rstatix")
library(dplyr)

###GOs can be changed to any GOs or list of OGs to run through statistical tests.
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
OGtitles<-c("GO:0007049","GO:0000910","GO:0006260","GO:0007163","GO:0007034","GO:0098609")

OGorder<-list()
p<-list()
plot_list<-list()


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
rownameorder<-rownames(plotcperkA)[row_order]
#vector1 <- as.vector(as.matrix(plotcperkA[row_order,]))
#vector2 <- as.vector(as.matrix(plotAnimalB[row_order,]))
#plot(vector1~vector2,xlim=c(-2,2),ylim=c(-2,2))
#abline(lm(vector1 ~ vector2),col="red")
#test_pearson <- cor.test(vector1, vector2, method = "pearson",alternative="greater")
#text(-1.8,1.7,p_round(test_pearson$p.value,2))

mat1<-plotcperkA[rownameorder,]
mat2<-plotAnimalB[rownameorder,]
##########################################
correlations <- numeric(nrow(mat1))
p_values <- numeric(nrow(mat2))

# Loop through each row
for(k in 1:nrow(mat1)) {
  # Perform correlation test on the ith row of each matrix
  test_result <- cor.test(mat1[k,],mat2[k,],alternative="greater")
  
  # Store the results
  correlations[k] <- test_result$estimate
  p_values[k] <- test_result$p.value
}

adjusted_p_values <- p.adjust(p_values, method = "BH")

# Count significant correlations
significant_correlations <- sum(adjusted_p_values < 0.05)

# Calculate average correlation coefficient
average_correlation <- mean(correlations)


sig_threshold <- 0.05
# Count of significant p-values
sig_count <- sum(adjusted_p_values < sig_threshold)
print(sig_count)
# Total number of p-values
total_count <- length(adjusted_p_values)
# Non-significant count (for illustration purposes)
non_sig_count <- total_count - sig_count
# Sets for the Venn diagram
total_set <- 1:total_count
sig_set <- 1:sig_count

#specify values to use in venn diagram
#fit <- euler(c('Cperk' = non_sig_count, 'Sarc' = non_sig_count, 'Cperk&Sarc' = sig_count))
#create venn diagram with custom colors

#file_name <- paste0("Euler_Diagram_", OGtitles[i], ".pdf")
#CairoPDF(file=file_name) 
#plot_list[[i]]<-grid.grabExpr(print(plot(fit, fill=c(brewer.pal(9,"Pastel1")[7], brewer.pal(9,"Pastel1")[8]),quantities=F) ))
#plot_list[[i]] <- recordPlot()

myps<-as.data.frame(-log(adjusted_p_values,base=10))
#myps<-as.data.frame(p_values)
myps$OG<-rep(OGtitles[i],length(adjusted_p_values))
myps$mycor<-correlations
colnames(myps)<-c("pval","OG","mycor")
rownames(myps)<-rownames(mat1)


p[[i]]<-ggplot(myps,aes(x=OG,y=pval)) + geom_violin(fill="lavenderblush") + geom_hline(yintercept=-log(0.05,base=10),linetype="dashed", color = "red", linewidth=1) + ylim(0,4) + ylab("-log10(adjusted pvalue)") + theme(axis.title.x = element_blank(), axis.title=element_text(size=6,face="bold"))

}

###Figure S2D--violin plots of cor.test significance values.
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],ncol=2)
