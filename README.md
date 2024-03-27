# Chromosphaera-perkinsii_Transcriptomics
Data and R scripts for direct and comparative transcriptomic analyses of the holozoan Chromosphaera perkinsii during its multicellular developmental period.

# R Libraries needed
All of these libraries need to be available in your system to run all of the code and generate all plots:
viridis,superheat,ggplot2,DBI,gridExtra,Cairo,scales,dplyr,tibble,zoo,matrixStats,tximport,edgeR,igraph,topGO,pheatmap,dendextend,eulerr,rstatix,plyr,gplots,purrr,reshape2,cowplot,RColorBrewer,grid,tidyr,data.table,progress

# This code and data is sufficient to make all base plots used in the manuscript. 
For some plots, colors and axis labels were added after plotting, so the colors and labels may not exactly match what is seen in the manuscript, however the data is the same.

# Before running any code, first expand the directory quants.zip 
quants.zip contains read mapping information from the RNAseq libraries from Salmon
Make sure that quants.zip is extracted such that the top directory is quants/ and inside that are the count directories. If it extracts with a nested "quants" directory (e.g. quants/quants/) it will not operate in downstream analyses.

# It may also be useful to run the code line by line or in chunks.

# Most plots will require information from one or both main analyses scripts: 
Cperk_TimeCourse_DEGanalyses.R, Cperk_TimeCourse_GOanalyses.R

# To make Figure 2A you can run:

setwd("Chromosphaera-perkinsii_Transcriptomics-main")

source("Cperk_TimeCourse_DEGanalyses.R")

source("MainFigureCode/Figure_2A.R")

# To make Figure 2B you need the output of the GO analyses. 
*The script Cperk_TimeCourse_GOanalyses.R is included for transparency, to produce GO enrichment data starting from the gene clusters. To save time for plotting, the output files from GO enrichment are already provided in directory "GOtables", so the running the GO enrichment script is not strictly necessary to generate the plots in Figure 2B.
**We normalized the scale for the enrichment plots in Figure 2B and manually colored bars with p-values smaller than 1x10^-11 dark blue. The raw plots color those p-values a dark gray.

source("Cperk_TimeCourse_DEGanalyses.R")

#uncomment the next line if you would like to run GO analyses--not strictly necessary for plotting Figure 2B since the output data files are already provided.
#source("Cperk_TimeCourse_GOanalyses.R")

source("MainFigureCode/Figure_2B.R")

# To make Figure 2C you can run:
*does not need any other data pre-loaded
**takes time to run through pattern clustering within orthogroups.

source("MainFigureCode/Figure_2C.R")

# To make Figure 2D you can run:

source("MainFigureCode/Figure_2D.R")

#The above code will print select transcription factors by default.
#To plot additional genes, after running the default, go to the following lines in the code and alter the functional class as described, or provide custom OGs for plotting:

###Here can select which functional class of gene to plot: TF, Adhesion, or Signaling. They are read from file DevelGenes_MAP.txt. Other OGs can be selected from annotation files for comparison. A subset is shown in Figure 2D in the manuscript.
GOgenes<-dgenes[dgenes[,3]=="TF",]


