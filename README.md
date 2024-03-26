# Chromosphaera-perkinsii_Transcriptomics
Data and R scripts for direct and comparative transcriptomic analyses of the holozoan Chromosphaera perkinsii during its multicellular developmental period.

# R Libraries needed
All of these libraries need to be available in your system to run all of the code and generate all plots:
viridis,superheat,ggplot2,DBI,gridExtra,Cairo,scales,dplyr,tibble,zoo,matrixStats,tximport,edgeR,igraph,topGO,pheatmap,dendextend,eulerr,rstatix,plyr,gplots,purrr,reshape2,cowplot,RColorBrewer,grid,tidyr,data.table,progress

# This code and data is sufficient to make all base plots used in the manuscript. 
For some plots, colors and axis labels were added after plotting, so the colors and labels may not exactly match what is seen in the manuscript, however the data is the same.

# To use the code it should be sufficient to navigate to the base directory and source each script. 
However, it may also be useful to run the code line by line or in chunks.

# Most plots will require information from the main analysis script: Cperk_TimeCourse_DEGanalyses.R
setwd("Chromosphaera-perkinsii_Transcriptomics")
source("Cperk_TimeCourse_DEGanalyses.R")

