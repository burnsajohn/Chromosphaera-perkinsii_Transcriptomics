
library(dplyr)
library(tibble)
library(zoo)
library(matrixStats)
library(scales)


###Blastula numbers
Cperk_numpatterns<-length(Cperk_OGplot[1,])-1
PT_OGplot <- animal_analysis("PT",14,Cperk_numpatterns)
FW_OGplot <- animal_analysis("FW",9,Cperk_numpatterns)
CE_OGplot <- animal_analysis("CE",20,Cperk_numpatterns)
HD_OGplot <- animal_analysis("HD",15,Cperk_numpatterns)
DM_OGplot <- animal_analysis("DM",9,Cperk_numpatterns)
SU_OGplot <- animal_analysis("SU",22,Cperk_numpatterns)
ZF_OGplot <- animal_analysis("ZF",6,Cperk_numpatterns)
NV_OGplot <- animal_analysis("NV",10,Cperk_numpatterns)
AQ_OGplot <- animal_analysis("AQ",11,Cperk_numpatterns)
ML_OGplot <- animal_analysis("ML",11,Cperk_numpatterns)

