###This code will the plots for Figure 2C from the manuscript. 
###It uses count data for C. perkinsii generated for this project and incoporates count data generated for Sphaeroforma arctica previously: https://elifesciences.org/articles/49801
###Figure 3—source data 2: https://cdn.elifesciences.org/articles/49801/elife-49801-fig3-data2-v4.txt


#libraries needed:
library(dplyr)
library(tibble)
library(zoo)
library(matrixStats)
library(scales)
library(tximport)
library(edgeR)
library(matrixStats)

####read in the C. perkinsii data:
sampletable <- read.table("Cperk_samples.txt", header = TRUE)
sampletable

files <- file.path("quants", sampletable$quant, "quant.sf")
filenames<-matrix(unlist(strsplit(files,"/")),ncol=21)[2,]

namearr<-vector(length=length(filenames))
for(i in 1:length(filenames)){
namearr[i]<-sampletable[match(filenames[i],sampletable[,3]),2]
}

names(files) <- namearr
txi.salmon <- tximport(files, type = "salmon", txOut=T)
head(txi.salmon$abundance)

###use TPMs to be compatible with 2012 paper.
Cperk.tpm<-txi.salmon$abundance

###remove "B" replicates
Cperk.tpm<-Cperk.tpm[,c(1,3,4,6,7,9,10,12,13,15,16,18,19,21)]

minsum<-min(Cperk.tpm[Cperk.tpm>0])/2
Cperk_logcounts<-log(Cperk.tpm+minsum, base=2)

###open map of C.perkinsii gene names to orthogroups defined across animals:
Cperk_OGs<-read.table("Orthogroups/Nk52Cperk_TRX_OGdf.txt",header=T)
Cperk_OGs2<-Cperk_OGs

#convert logcount data to data frame and split out A and C replicates
Cperk_OGplot<-as.data.frame(Cperk_logcounts)
CperkA_OGplot<-Cperk_OGplot[,c(1,3,5,7,9,11,13)]
CperkC_OGplot<-Cperk_OGplot[,c(2,4,6,8,10,12,14)]

###format the dataframes for mapping OGs to count data
CperkA_OGplot <- CperkA_OGplot %>% rownames_to_column("Nk52Cperk_TRX")
CperkA_OGplot <- CperkA_OGplot %>% left_join(Cperk_OGs2, by = "Nk52Cperk_TRX")
CperkC_OGplot <- CperkC_OGplot %>% rownames_to_column("Nk52Cperk_TRX")
CperkC_OGplot <- CperkC_OGplot %>% left_join(Cperk_OGs2, by = "Nk52Cperk_TRX")

rownames(CperkA_OGplot)<-CperkA_OGplot$Nk52Cperk_TRX
CperkA_OGplot_OG<-CperkA_OGplot[, !(colnames(CperkA_OGplot) %in% c("Nk52Cperk_TRX", "vars"))]
rownames(CperkC_OGplot)<-CperkC_OGplot$Nk52Cperk_TRX
CperkC_OGplot_OG<-CperkC_OGplot[, !(colnames(CperkC_OGplot) %in% c("Nk52Cperk_TRX", "vars"))]
   
CperkA_OGplot<-t(scale(t(CperkA_OGplot_OG[, !(colnames(CperkA_OGplot_OG) %in% c("Orthogroup"))])))
CperkC_OGplot<-t(scale(t(CperkC_OGplot_OG[, !(colnames(CperkC_OGplot_OG) %in% c("Orthogroup"))])))

CperkA_OGplot<-as.data.frame(CperkA_OGplot)
CperkC_OGplot<-as.data.frame(CperkC_OGplot)

