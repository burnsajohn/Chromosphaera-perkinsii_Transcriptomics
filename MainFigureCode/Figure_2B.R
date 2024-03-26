###This code will make a plot similar to Figure 2B from the manuscript. 

###Can be run with provided GO data tables in the /GOtables directory

#Or: all required data can be generated using the two following scripts
#First run: "Cperk_TimeCourse_DEGanalyses.R" to create scaled gene matrix and clustering objects
#Next run: "Cperk_TimeCourse_GOanalyses.R" to create GO tables for plotting.

#libraries needed
library(ggplot2)
library(DBI)
library(gridExtra)
library(Cairo)

###clusters renumbered based on expression patterns after initial random assignment (reproducible seed=20).
#Old Cluster 5 becomes Cluster 1
#Old Cluster 2 remains Cluster 2
#Old Cluster 1 becomes Cluster 3
#Old Cluster 3 becomes Cluster 4
#Old Cluster 4 becomes Cluster 5

###read in GO tables for plotting
c1<-read.table("GOtables/5_GOPBtbls.txt",sep="\t",header=T)
c2<-read.table("GOtables/2_GOPBtbls.txt",sep="\t",header=T)
c3<-read.table("GOtables/1_GOPBtbls.txt",sep="\t",header=T)
c4<-read.table("GOtables/3_GOPBtbls.txt",sep="\t",header=T)
c5<-read.table("GOtables/4_GOPBtbls.txt",sep="\t",header=T)

all_labels<-rbind(c1[1:5,],c2[1:5,],c3[1:5,],c4[1:5,],c5[1:5,])$Term

###define a function to plot the GO data
plotOGs<-function(x){
y = "Ratio"
df<-x[1:5,]
df$Ratio<-df$Significant/df$Annotated
colnames(df)[4]<-"Count"
colnames(df)[6]<-"Scores"
df$Scores<-as.numeric(gsub("< ","",df$Scores))

# Find the length of the longest label
max_length <- max(nchar(as.character(all_labels)))
  
# Pad all labels to have the same length, right-justified, using a specific character (e.g., '.')
df$Term <- factor(sapply(as.character(df$Term), function(t) {
  padding_length <- max_length - nchar(t)
  padding <- paste(rep(".", padding_length), collapse = "")
  paste0(padding, t)
}))

# Define variable mapping
map <- aes_string(x = "Term", y = y, fill = "Scores")
df$Term <- factor(df$Term, levels = df$Term[order(-df$Scores)])
library(scales)
my_breaks = c(1e-1, 1e-3, 1e-5, 1e-7, 1e-9, 1e-11)

# Make the ggplot
p <- ggplot(df, map) + geom_bar(stat = "identity") + coord_flip() + theme_bw() + ylim(0,1) +
scale_fill_gradient2(low = "blue",mid = "lightblue", high = "red",trans = "log",midpoint=-10,breaks=my_breaks,limits=c(1e-11, 1e-1)) +
guides(fill = guide_colorbar(title = "", reverse = TRUE))+
theme(legend.key.size = unit(3, 'mm'))

# Adjust theme components
p <- p + theme(axis.text.x = element_text(colour = "black", vjust = 1), 
text = element_text(family = "Courier New"), # Set font family to Courier
axis.text.y = element_text(colour = "black", hjust = 1),
axis.title = element_text(color = "black", margin = margin(10, 5, 0, 0)),
axis.title.y = element_text(angle = 90))

print(p)
}

p1<-plotOGs(c1)
p2<-plotOGs(c2)
p3<-plotOGs(c3)
p4<-plotOGs(c4)
p5<-plotOGs(c5)

#plot to screen using grid.arrange (from gridExtra)
grid.arrange(p1, p2, p3, p4, p5, nrow = 5)

###make a pdf--aligns plots.
CairoPDF("GOplots.pdf", width = 11, height = 8.5) # Adjust size as needed

# Arrange and plot using grid.arrange
grid.arrange(p1, p2, p3, p4, p5, nrow = 5)

# Close the PDF device
dev.off()
