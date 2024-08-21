# DOI: 10.5281/zenodo.13352464

# Chromosphaera-perkinsii_Transcriptomics
Data and R scripts for direct and comparative transcriptomic analyses of the holozoan Chromosphaera perkinsii during its multicellular developmental period.
*Note: Most scripts produce plots to the screen, not directly to pdfs. Those plots can be savedexported as PDFs, or the code can be altered to export pdfs directly, as desired.

# R Libraries needed to run these scripts
viridis, superheat, ggplot2, DBI, gridExtra, Cairo, scales, dplyr, tibble, zoo, matrixStats, tximport, edgeR, igraph, topGO, pheatmap, dendextend, eulerr, rstatix, plyr, gplots, purrr, reshape2, cowplot, RColorBrewer, grid, tidyr, data.table, progress, ggridges, corrplot, irr

#dependencies for individual scripts are listed near the top of the code files. Version info at time of publication can be found at the bottom of this README. Installation time for all dependencies should be around 2 hours.

# This code and data is sufficient to make all base plots used in the manuscript. 
For some plots, colors and axis labels were added after plotting, so the colors and labels may not exactly match what is seen in the manuscript, however the data is the same.

Run time for each plot varies from seconds to approximately 30 minutes.

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

# Supplementary figures can be generated using a similar code snippets: either by sourcing the plotting files, or by running the code by line or in chunks.

source("SupplementalFigureCode/Figure_S2a.R")

source("SupplementalFigureCode/Figure_S2b.R")

# For Figure S2C, source the script for Figure 2C:
source("MainFigureCode/Figure_2C.R")
#Then, alter the GO terms that you would like to plot for comparison: In the section labeled:

########################GOs for comparison########################################################################

#Change the GO terms to those for Figure S2C, or any other categories you would like to compare between C. perkinsii and S. arctica as developmental expression patterns.

#Figure S2C used: GO:0000226 (MT Organization); GO:0030154 (Cell differentiation); GO:0051301 (Cell Division)

source("SupplementalFigureCode/Figure_S2d.R")

source("SupplementalFigureCode/Figure_S2e.R")

source("SupplementalFigureCode/Figure_S2f.R")

#To vary Figure S2F you can vary the GO term genes are selected from and can vary the number of genes plotted to visualize different sets of genes. 

source("SupplementalFigureCode/Figure_S2g.R")

source("SupplementalFigureCode/Figure_S2h.R")

source("SupplementalFigureCode/Figure_S2i.R")

source("SupplementalFigureCode/Figure_S2j.R")

source("SupplementalFigureCode/Figure_S2k.R")

source("SupplementalFigureCode/Figure_S2l.R")

#code for calculating Kendall's p-value for correspondence among animals for individual genes can be found in the file "KendallsPvalueCalcs.R" in the SupplementalFigureCode directory.

# Package Version Info at the time of publication:
 [1]irr_0.84.1           lpSolve_5.6.20       corrplot_0.92       
 [4] ggridges_0.5.6       progress_1.2.3       data.table_1.15.4   
 [7] tidyr_1.3.0          RColorBrewer_1.1-3   cowplot_1.1.3       
[10] reshape2_1.4.4       purrr_1.0.2          gplots_3.1.3.1      
[13] plyr_1.8.9           rstatix_0.7.2        eulerr_7.0.2        
[16] dendextend_1.17.1    pheatmap_1.0.12      topGO_2.54.0        
[19] SparseM_1.81         GO.db_3.18.0         AnnotationDbi_1.64.1
[22] IRanges_2.36.0       S4Vectors_0.40.2     Biobase_2.62.0      
[25] graph_1.80.0         BiocGenerics_0.48.1  igraph_2.0.3        
[28] edgeR_4.0.3          limma_3.58.1         tximport_1.30.0     
[31] matrixStats_1.3.0    zoo_1.8-12           tibble_3.2.1        
[34] dplyr_1.1.4          scales_1.3.0         Cairo_1.6-2         
[37] gridExtra_2.3        DBI_1.2.2            ggplot2_3.5.1       
[40] superheat_1.0.0      viridis_0.6.5        viridisLite_0.4.2 





