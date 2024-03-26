##First run analyses from: "Cperk_TimeCourse_DEGanalyses.R"
##Then run this script to make the raw plots for Figure S2B. Boxplots were recolored manully outside of R for the manuscript. 

#######facet plot all clusters#########################
i=1
#Subset the cores molten dataframe so we can plot the core
core2 <- Kmolten[Kmolten$cluster==i,]
core2$CLUSTER<-i

#get cluster 1
K2 <- (scaledata[names(kClusters[kClusters==i]),])
#calculate the correlation with the core
corscore <- function(x){cor(x,core2$value)}
score <- apply(K2, 1, corscore)
#get the data frame into long format for plotting
K2molten <- melt(K2)
colnames(K2molten) <- c('gene','sample','value')
#add the score
K2molten <- merge(K2molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K2molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K2molten$order_factor <- 1:length(K2molten$gene)
#order the dataframe by score
K2molten <- K2molten[order(K2molten$score),]
#set the order by setting the factors
K2molten$order_factor <- factor(K2molten$order_factor , levels = K2molten$order_factor)
K2molten$CLUSTER<-i
plotclusters<-K2molten
coreclusters<-core2

for(i in 2:length(unique(Kmolten$cluster))){
#Subset the cores molten dataframe so we can plot the core
core2 <- Kmolten[Kmolten$cluster==i,]
core2$CLUSTER<-i
#get cluster 2
K2 <- (scaledata[names(kClusters[kClusters==i]),])
#calculate the correlation with the core
corscore <- function(x){cor(x,core2$value)}
score <- apply(K2, 1, corscore)
#get the data frame into long format for plotting
K2molten <- melt(K2)
colnames(K2molten) <- c('gene','sample','value')
#add the score
K2molten <- merge(K2molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K2molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K2molten$order_factor <- 1:length(K2molten$gene)
#order the dataframe by score
K2molten <- K2molten[order(K2molten$score),]
#set the order by setting the factors
K2molten$order_factor <- factor(K2molten$order_factor , levels = K2molten$order_factor)
K2molten$CLUSTER<-i

plotclusters<-rbind(plotclusters,K2molten)
coreclusters<-rbind(coreclusters,core2)

}


# Everything on the same plot
plotclusters2<-plotclusters

plotclusters2$time <- as.numeric(sapply(strsplit(as.character(plotclusters2$sample), "_"), function(x) as.numeric(x[2])))

# Assuming your data is in a dataframe called plotclusters2
# Calculate medians
medians <- plotclusters2 %>% dplyr::group_by(CLUSTER, time) %>% dplyr::summarise(median_value = median(value))


# Merge medians back with original data for plotting
plot_data <- merge(plotclusters2, medians, by = c("CLUSTER", "time"))

###Figure S2B--clusters were renamed/reordered for the manuscript
# Create the plot
p2 <- ggplot(plot_data, aes(x = as.factor(time), y = value)) +
  geom_boxplot() +
  geom_line(data = medians, aes(x = as.factor(time), y = median_value, group = CLUSTER), color = "blue") +
  facet_grid(cols = vars(CLUSTER)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Print the plot
print(p2)

