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
Cperk_OGs<-read.table("Nk52Cperk_TRX_OGdf.txt",header=T)


