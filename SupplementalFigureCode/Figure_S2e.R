###Script to produce Figure S2E--correlation plot of specific genes identified as important for S. arctica development in https://elifesciences.org/articles/49801 (Figure 5)

###First run: Cperk_TimeCourse_DEGanalyses.R to read in C. perkinsii data and format data structures for this analysis.

source("Scripts/animal_analysis2.rscript")
###elifeFig5 analyses
Cperk_noreps<-average_columns_pairs(logcounts)
CPerc_scaledataFULL <- t(scale(t(Cperk_noreps),scale=T))
Cperk_sd_cut<-CPerc_scaledataFULL[,c(1:6)]

genes_df<-read.table("Annotations/elife_fig5_genemap.txt")
colnames(genes_df)<-c("GeneName","Sarc","Cperk")

# Initialize a vector to store the correlation for each gene
correlations <- numeric(nrow(genes_df))

###############getSARCdata######################################
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

Sarc_logcounts<-cbind(rowMeans(Sarc_logcounts[,c(1,11)]),rowMeans(Sarc_logcounts[,c(2,12)]),rowMeans(Sarc_logcounts[,c(3,13)]),rowMeans(Sarc_logcounts[,c(4,14)]),rowMeans(Sarc_logcounts[,c(5,15)]),rowMeans(Sarc_logcounts[,c(6,16)]),Sarc_logcounts[,7],rowMeans(Sarc_logcounts[,c(8,17)]),rowMeans(Sarc_logcounts[,c(9,18)]),rowMeans(Sarc_logcounts[,c(10,19)]),Sarc_logcounts[,20])
colnames(Sarc_logcounts)<-c("T12","T18","T24","T30","T36","T42","T48","T54","T60","T66","T72")


ntimes="all"
if(ntimes=="all"){ntimes=length(Sarc_logcounts[1,])}

Cperk_numpatterns<-length(Cperk_sd_cut[1,])
Sarc_logcounts_red <- adaptive_rolling_average_across_columns(Sarc_logcounts[,1:ntimes], Cperk_numpatterns,1)
Sarc_scale <- t(scale(t(Sarc_logcounts_red),scale=T))


# Loop through each row in genes_df to calculate correlation
for(i in 1:nrow(genes_df)) {
	# Extract the identifiers for the current row
	identifier1 <- genes_df[i, "Sarc"]
	identifier2 <- genes_df[i, "Cperk"]
	if(identifier1 %in% rownames(Sarc_scale) && identifier2 %in% rownames(Cperk_sd_cut)) {  
		# Fetch the corresponding vectors from df1 and df2
		vector1 <- Sarc_scale[identifier1,]
		vector2 <- Cperk_sd_cut[identifier2,]
  
		# Calculate the correlation between the vectors
		correlations[i] <- cor(vector1, vector2, use = "complete.obs") # 'use' parameter handles missing values
	}else {
    # If the identifiers do not exist, set the correlation to NA or some other indicator
    correlations[i] <- NA
	}
}
# Add the correlation values to the genes_df
genes_df$Correlation <- correlations

# Filter out rows where the correlation is NA
genes_df <- genes_df[!is.na(genes_df$Correlation), ]

# Now you can plot these correlations with the gene names as labels
print(ggplot(genes_df, aes(x = reorder(GeneName, Correlation), y = Correlation, label = GeneName)) +
  geom_point(stat = "identity") +
  #geom_text(nudge_y = 0.05, check_overlap = TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Gene Name", y = "Correlation between vectors"))
