#this code will calculate statistics on correlation distributions between C. perkinsii and each animal of genes not involved in development in C. perkinsii and each cluster of genes shown to be important for development in C. perkinsii. It will calculate statistics as reported in the manuscript and indicated in the code and will generate the plot shown in Figure S2i (the plot will be identical except that coloring was added manually after plotting for each cluster for the manuscript).


###First you must run: 
#source("Cperk_TimeCourse_DEGanalyses.R")
#source("Figure_2D.R")

#Library needed:
library(ggridges)

#look at clusters vs. nottimeDEs
Cperk_timeOGs<-unique(na.omit(Cperk_OGs2[match(timeDEs,Cperk_OGs2[,2]),1]))
Cperk_notTimeOGs<-setdiff(Cperk_OGs2[,1],Cperk_timeOGs)

# Define your desired cluster order
desired_cluster_order <- c(4,3,1,2,5)

GOogs.nt<-unique(Cperk_notTimeOGs) ###not DE expressed in development
Cperk.k1<-names(kClusters[kClusters==4])
GOogs.k1<-unique(na.omit(Cperk_OGs2[match(Cperk.k1,Cperk_OGs2[,2]),1]))
Cperk.k2<-names(kClusters[kClusters==3])
GOogs.k2<-unique(na.omit(Cperk_OGs2[match(Cperk.k2,Cperk_OGs2[,2]),1]))
Cperk.k3<-names(kClusters[kClusters==1])
GOogs.k3<-unique(na.omit(Cperk_OGs2[match(Cperk.k3,Cperk_OGs2[,2]),1]))
Cperk.k4<-names(kClusters[kClusters==2])
GOogs.k4<-unique(na.omit(Cperk_OGs2[match(Cperk.k4,Cperk_OGs2[,2]),1]))
Cperk.k5<-names(kClusters[kClusters==5])
GOogs.k5<-unique(na.omit(Cperk_OGs2[match(Cperk.k5,Cperk_OGs2[,2]),1]))

testcors<-list(GOogs.nt,GOogs.k1,GOogs.k2,GOogs.k3,GOogs.k4,GOogs.k5)

cors.list<-list()
for(tc in 1:length(testcors)){
GOogs<-testcors[[tc]]

# Assuming OG_matList is a list of lists of matrices, and we're focusing on the first list of matrices
matrices <- OG_matList[[1]]  # Extract the first list of matrices

###lets make sure to look at all patterns per OG!
base_ids <- ifelse(grepl("\\.", GOogs), sub("\\..*", "", GOogs), GOogs)
patterns <- paste0("^", base_ids, "(\\.\\d+)?$")
OGpats <- unlist(lapply(patterns, function(pat) grep(pat, rownames(matrices[[1]]), value = TRUE)))

if(length(OGpats)<2){next}

# Subset all matrices by the selected rownames
basals <- lapply(matrices, function(x) x[OGpats, ])


# Initialize a vector to store the average correlations
correlations <- c()

####compute pairwise correlations across every row for Cperk vs. each animal in "basals":
nrow(basals[[1]])
# Loop over each row
for (i in 1:nrow(basals[[1]])) {
  # Compute pairwise correlations
  for (j in 2:11) {
      row1 <- basals[[1]][i, ]
      row2 <- basals[[j]][i, ]
	  if(all(is.na(row1)) || all(is.na(row2))){
		correlations <- c(correlations, NA)
	  }else{
      correlations <- c(correlations, cor(row1, row2,method="pearson"))
	  }
  }
}
cors.list[[tc]]<-correlations
}

names(cors.list)<-c("notDE","Clust1","Clust2","Clust3","Clust4","Clust5")


# Combine the list into a data frame
combined_data <- do.call(rbind, lapply(names(cors.list), function(group_name) {
  data.frame(value = cors.list[[group_name]], group = factor(group_name))
}))

# Perform the Kruskal-Wallis test
kruskal_result <- kruskal.test(value ~ group, data = combined_data)

# Print the summary of the Kruskal-Wallis test--this is the overall significance value indicating whether there are any significant differences in the data.
print(kruskal_result)

# Perform pairwise one-sided Wilcoxon rank-sum tests
pairwise_results <- lapply(names(cors.list)[-1], function(test_group) {
  # Subset the data for control and the current test group
  subset_data <- combined_data[combined_data$group %in% c("notDE", test_group), ]
  
  # Ensure 'control' is the first group and test_group is the second group
  subset_data$group <- factor(subset_data$group, levels = c(test_group,"notDE"))
  
  # Perform one-sided Wilcoxon test (greater)
  wilcox_result <- wilcox.test(value ~ group, data = subset_data, alternative = "greater")
  
  # Return the result
  list(group = test_group, p.value = wilcox_result$p.value)
})

# Adjust p-values for multiple comparisons
adjusted_p.values <- p.adjust(sapply(pairwise_results, `[[`, "p.value"), method = "bonferroni")

# Combine results into a data frame for easier viewing
final_results <- data.frame(
  group = sapply(pairwise_results, `[[`, "group"),
  p.value = sapply(pairwise_results, `[[`, "p.value"),
  adjusted.p.value = adjusted_p.values
)

# Print the final results--this shows p-values for pairwise comparisons of notDE to each cluster.
print(final_results)

##set the order of the levels for plotting
combined_data$group <- factor(combined_data$group, levels = c("Clust5", "Clust4", "Clust3","Clust2","Clust1","notDE"))

##Figure S2i--individual colors for each cluster were added manually after plotting, but could also be coded into this function.
ggplot(combined_data, aes(x = value, y = group)) +
  stat_density_ridges(fill = "lightblue", quantile_lines = TRUE, alpha = 0.5, quantiles = 2)  +
  geom_vline(xintercept = median(combined_data[combined_data$group=="notDE","value"],na.rm=T), linetype = "dashed", color = "red") + # Change xintercept to your desired value
  theme_ridges() + 
  labs(title = "", x = "Correlation between C. perkinsii and animals", y = "")


