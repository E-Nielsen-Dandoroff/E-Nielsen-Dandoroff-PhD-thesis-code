# 24/06/2024

# EntireCohort_pathvar_CentCytgenes.R

# this script is for the Centrosomal and cytoskeletal associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)


# CentCytgene
# POC1A
# CDK6
# ARHGEF2
# ASPM
# BUB1
# CDK5RAP2
# CENPE
# CENPF
# CENPJ
# CENPT
# CEP135
# CEP152
# CEP63
# CIT
# CKAP2L
# KATNB1
# KIF11
# KIF22
# KIFBP
# KNL1
# NDE1
# NIN
# PCNT
# PLK4
# RTTN
# SASS6
# STIL
# TRAPPC14
# TUBGCP6
# WDR62
# NEK1
# CFAP410
# LMNA
# LMNB1
# LMNB2
# ZMPSTE24
# FLNA
# FLNB
# TUBA1A
# TUBB2B
# TUBGCP4
# ACTB
# ACTG1
# MACF1
# VPS4A












# CentCytgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')


# go thr all the genes to check that they're present
grep("VPS4A", EntireCohort_pathvar$GENE)
# all seem to be present, except for:
# KIF22 = 0
# LMNA = 0
# LMNB1 = 0
# LMNB2 = 0
# FLNA = 0
# TUBA1A = 0
# ACTG1 = 0
# MACF1 = 0
# VPS4A = 0


# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$CentCytgene <- NA
EntireCohort_pathvar$CentCytgene[EntireCohort_pathvar$GENE1 %in% c('POC1A', 'CDK6', 'ARHGEF2', 'ASPM', 'BUB1', 'CDK5RAP2', 'CENPE', 'CENPF', 'CENPJ', 'CENPT', 'CEP135', 'CEP152', 'CEP63', 'CIT', 'CKAP2L', 'KATNB1', 'KIF11', 'KIFBP', 'KNL1', 'NDE1', 'NIN', 'PCNT', 'PLK4', 'RTTN', 'SASS6', 'STIL', 'TRAPPC14', 'TUBGCP6', 'WDR62', 'NEK1', 'CFAP410', 'ZMPSTE24', 'FLNB', 'TUBB2B', 'TUBGCP4', 'ACTB')] <- 1
EntireCohort_pathvar$CentCytgene[EntireCohort_pathvar$GENE2 %in% c('POC1A', 'CDK6', 'ARHGEF2', 'ASPM', 'BUB1', 'CDK5RAP2', 'CENPE', 'CENPF', 'CENPJ', 'CENPT', 'CEP135', 'CEP152', 'CEP63', 'CIT', 'CKAP2L', 'KATNB1', 'KIF11', 'KIFBP', 'KNL1', 'NDE1', 'NIN', 'PCNT', 'PLK4', 'RTTN', 'SASS6', 'STIL', 'TRAPPC14', 'TUBGCP6', 'WDR62', 'NEK1', 'CFAP410', 'ZMPSTE24', 'FLNB', 'TUBB2B', 'TUBGCP4', 'ACTB')] <- 1
EntireCohort_pathvar$CentCytgene[EntireCohort_pathvar$GENE3 %in% c('POC1A', 'CDK6', 'ARHGEF2', 'ASPM', 'BUB1', 'CDK5RAP2', 'CENPE', 'CENPF', 'CENPJ', 'CENPT', 'CEP135', 'CEP152', 'CEP63', 'CIT', 'CKAP2L', 'KATNB1', 'KIF11', 'KIFBP', 'KNL1', 'NDE1', 'NIN', 'PCNT', 'PLK4', 'RTTN', 'SASS6', 'STIL', 'TRAPPC14', 'TUBGCP6', 'WDR62', 'NEK1', 'CFAP410', 'ZMPSTE24', 'FLNB', 'TUBB2B', 'TUBGCP4', 'ACTB')] <- 1
EntireCohort_pathvar$CentCytgene[EntireCohort_pathvar$GENE4 %in% c('POC1A', 'CDK6', 'ARHGEF2', 'ASPM', 'BUB1', 'CDK5RAP2', 'CENPE', 'CENPF', 'CENPJ', 'CENPT', 'CEP135', 'CEP152', 'CEP63', 'CIT', 'CKAP2L', 'KATNB1', 'KIF11', 'KIFBP', 'KNL1', 'NDE1', 'NIN', 'PCNT', 'PLK4', 'RTTN', 'SASS6', 'STIL', 'TRAPPC14', 'TUBGCP6', 'WDR62', 'NEK1', 'CFAP410', 'ZMPSTE24', 'FLNB', 'TUBB2B', 'TUBGCP4', 'ACTB')] <- 1
EntireCohort_pathvar$CentCytgene[EntireCohort_pathvar$GENE5 %in% c('POC1A', 'CDK6', 'ARHGEF2', 'ASPM', 'BUB1', 'CDK5RAP2', 'CENPE', 'CENPF', 'CENPJ', 'CENPT', 'CEP135', 'CEP152', 'CEP63', 'CIT', 'CKAP2L', 'KATNB1', 'KIF11', 'KIFBP', 'KNL1', 'NDE1', 'NIN', 'PCNT', 'PLK4', 'RTTN', 'SASS6', 'STIL', 'TRAPPC14', 'TUBGCP6', 'WDR62', 'NEK1', 'CFAP410', 'ZMPSTE24', 'FLNB', 'TUBB2B', 'TUBGCP4', 'ACTB')] <- 1
EntireCohort_pathvar$CentCytgene[EntireCohort_pathvar$GENE6 %in% c('POC1A', 'CDK6', 'ARHGEF2', 'ASPM', 'BUB1', 'CDK5RAP2', 'CENPE', 'CENPF', 'CENPJ', 'CENPT', 'CEP135', 'CEP152', 'CEP63', 'CIT', 'CKAP2L', 'KATNB1', 'KIF11', 'KIFBP', 'KNL1', 'NDE1', 'NIN', 'PCNT', 'PLK4', 'RTTN', 'SASS6', 'STIL', 'TRAPPC14', 'TUBGCP6', 'WDR62', 'NEK1', 'CFAP410', 'ZMPSTE24', 'FLNB', 'TUBB2B', 'TUBGCP4', 'ACTB')] <- 1
EntireCohort_pathvar$CentCytgene[EntireCohort_pathvar$GENE7 %in% c('POC1A', 'CDK6', 'ARHGEF2', 'ASPM', 'BUB1', 'CDK5RAP2', 'CENPE', 'CENPF', 'CENPJ', 'CENPT', 'CEP135', 'CEP152', 'CEP63', 'CIT', 'CKAP2L', 'KATNB1', 'KIF11', 'KIFBP', 'KNL1', 'NDE1', 'NIN', 'PCNT', 'PLK4', 'RTTN', 'SASS6', 'STIL', 'TRAPPC14', 'TUBGCP6', 'WDR62', 'NEK1', 'CFAP410', 'ZMPSTE24', 'FLNB', 'TUBB2B', 'TUBGCP4', 'ACTB')] <- 1
EntireCohort_pathvar$CentCytgene[EntireCohort_pathvar$GENE8 %in% c('POC1A', 'CDK6', 'ARHGEF2', 'ASPM', 'BUB1', 'CDK5RAP2', 'CENPE', 'CENPF', 'CENPJ', 'CENPT', 'CEP135', 'CEP152', 'CEP63', 'CIT', 'CKAP2L', 'KATNB1', 'KIF11', 'KIFBP', 'KNL1', 'NDE1', 'NIN', 'PCNT', 'PLK4', 'RTTN', 'SASS6', 'STIL', 'TRAPPC14', 'TUBGCP6', 'WDR62', 'NEK1', 'CFAP410', 'ZMPSTE24', 'FLNB', 'TUBB2B', 'TUBGCP4', 'ACTB')] <- 1
EntireCohort_pathvar$CentCytgene[EntireCohort_pathvar$GENE9 %in% c('POC1A', 'CDK6', 'ARHGEF2', 'ASPM', 'BUB1', 'CDK5RAP2', 'CENPE', 'CENPF', 'CENPJ', 'CENPT', 'CEP135', 'CEP152', 'CEP63', 'CIT', 'CKAP2L', 'KATNB1', 'KIF11', 'KIFBP', 'KNL1', 'NDE1', 'NIN', 'PCNT', 'PLK4', 'RTTN', 'SASS6', 'STIL', 'TRAPPC14', 'TUBGCP6', 'WDR62', 'NEK1', 'CFAP410', 'ZMPSTE24', 'FLNB', 'TUBB2B', 'TUBGCP4', 'ACTB')] <- 1


EntireCohort_pathvar$CentCytgene[is.na(EntireCohort_pathvar$CentCytgene)] <- 0
EntireCohort_pathvar$CentCytgene <- as.numeric(EntireCohort_pathvar$CentCytgene)
unique(EntireCohort_pathvar$CentCytgene)
# 0 1

table(EntireCohort_pathvar$CentCytgene)
#  0      1 
# 404209  65590 




# now to make a column for each of the CentCytgenes
### 'POC1A', 
EntireCohort_pathvar$POC1A <- NA
EntireCohort_pathvar$POC1A[EntireCohort_pathvar$GENE1 %in% c("POC1A")] <- 1
EntireCohort_pathvar$POC1A[EntireCohort_pathvar$GENE2 %in% c("POC1A")] <- 1
EntireCohort_pathvar$POC1A[EntireCohort_pathvar$GENE3 %in% c("POC1A")] <- 1
EntireCohort_pathvar$POC1A[EntireCohort_pathvar$GENE4 %in% c("POC1A")] <- 1
EntireCohort_pathvar$POC1A[EntireCohort_pathvar$GENE5 %in% c("POC1A")] <- 1
EntireCohort_pathvar$POC1A[EntireCohort_pathvar$GENE6 %in% c("POC1A")] <- 1
EntireCohort_pathvar$POC1A[EntireCohort_pathvar$GENE7 %in% c("POC1A")] <- 1
EntireCohort_pathvar$POC1A[EntireCohort_pathvar$GENE8 %in% c("POC1A")] <- 1
EntireCohort_pathvar$POC1A[EntireCohort_pathvar$GENE9 %in% c("POC1A")] <- 1

EntireCohort_pathvar$POC1A[is.na(EntireCohort_pathvar$POC1A)] <- 0
EntireCohort_pathvar$POC1A <- as.numeric(EntireCohort_pathvar$POC1A)
unique(EntireCohort_pathvar$POC1A)
# 0 1

table(EntireCohort_pathvar$POC1A)
# 0      1 
# 466825   2974





### 'CDK6', 
EntireCohort_pathvar$CDK6 <- NA
EntireCohort_pathvar$CDK6[EntireCohort_pathvar$GENE1 %in% c("CDK6")] <- 1
EntireCohort_pathvar$CDK6[EntireCohort_pathvar$GENE2 %in% c("CDK6")] <- 1
EntireCohort_pathvar$CDK6[EntireCohort_pathvar$GENE3 %in% c("CDK6")] <- 1
EntireCohort_pathvar$CDK6[EntireCohort_pathvar$GENE4 %in% c("CDK6")] <- 1
EntireCohort_pathvar$CDK6[EntireCohort_pathvar$GENE5 %in% c("CDK6")] <- 1
EntireCohort_pathvar$CDK6[EntireCohort_pathvar$GENE6 %in% c("CDK6")] <- 1
EntireCohort_pathvar$CDK6[EntireCohort_pathvar$GENE7 %in% c("CDK6")] <- 1
EntireCohort_pathvar$CDK6[EntireCohort_pathvar$GENE8 %in% c("CDK6")] <- 1
EntireCohort_pathvar$CDK6[EntireCohort_pathvar$GENE9 %in% c("CDK6")] <- 1

EntireCohort_pathvar$CDK6[is.na(EntireCohort_pathvar$CDK6)] <- 0
EntireCohort_pathvar$CDK6 <- as.numeric(EntireCohort_pathvar$CDK6)
unique(EntireCohort_pathvar$CDK6)
# 0 1

table(EntireCohort_pathvar$CDK6)
# 0      1 
# 469671    128




### 'ARHGEF2', 
EntireCohort_pathvar$ARHGEF2 <- NA
EntireCohort_pathvar$ARHGEF2[EntireCohort_pathvar$GENE1 %in% c("ARHGEF2")] <- 1
EntireCohort_pathvar$ARHGEF2[EntireCohort_pathvar$GENE2 %in% c("ARHGEF2")] <- 1
EntireCohort_pathvar$ARHGEF2[EntireCohort_pathvar$GENE3 %in% c("ARHGEF2")] <- 1
EntireCohort_pathvar$ARHGEF2[EntireCohort_pathvar$GENE4 %in% c("ARHGEF2")] <- 1
EntireCohort_pathvar$ARHGEF2[EntireCohort_pathvar$GENE5 %in% c("ARHGEF2")] <- 1
EntireCohort_pathvar$ARHGEF2[EntireCohort_pathvar$GENE6 %in% c("ARHGEF2")] <- 1
EntireCohort_pathvar$ARHGEF2[EntireCohort_pathvar$GENE7 %in% c("ARHGEF2")] <- 1
EntireCohort_pathvar$ARHGEF2[EntireCohort_pathvar$GENE8 %in% c("ARHGEF2")] <- 1
EntireCohort_pathvar$ARHGEF2[EntireCohort_pathvar$GENE9 %in% c("ARHGEF2")] <- 1

EntireCohort_pathvar$ARHGEF2[is.na(EntireCohort_pathvar$ARHGEF2)] <- 0
EntireCohort_pathvar$ARHGEF2 <- as.numeric(EntireCohort_pathvar$ARHGEF2)
unique(EntireCohort_pathvar$ARHGEF2)
# 0 1

table(EntireCohort_pathvar$ARHGEF2)
# 0      1 
# 468776   1023




### 'ASPM',  
EntireCohort_pathvar$ASPM <- NA
EntireCohort_pathvar$ASPM[EntireCohort_pathvar$GENE1 %in% c("ASPM")] <- 1
EntireCohort_pathvar$ASPM[EntireCohort_pathvar$GENE2 %in% c("ASPM")] <- 1
EntireCohort_pathvar$ASPM[EntireCohort_pathvar$GENE3 %in% c("ASPM")] <- 1
EntireCohort_pathvar$ASPM[EntireCohort_pathvar$GENE4 %in% c("ASPM")] <- 1
EntireCohort_pathvar$ASPM[EntireCohort_pathvar$GENE5 %in% c("ASPM")] <- 1
EntireCohort_pathvar$ASPM[EntireCohort_pathvar$GENE6 %in% c("ASPM")] <- 1
EntireCohort_pathvar$ASPM[EntireCohort_pathvar$GENE7 %in% c("ASPM")] <- 1
EntireCohort_pathvar$ASPM[EntireCohort_pathvar$GENE8 %in% c("ASPM")] <- 1
EntireCohort_pathvar$ASPM[EntireCohort_pathvar$GENE9 %in% c("ASPM")] <- 1

EntireCohort_pathvar$ASPM[is.na(EntireCohort_pathvar$ASPM)] <- 0
EntireCohort_pathvar$ASPM <- as.numeric(EntireCohort_pathvar$ASPM)
unique(EntireCohort_pathvar$ASPM)
# 0 1

table(EntireCohort_pathvar$ASPM)
#   0      1 
# 460203   9596




### 'BUB1', 
EntireCohort_pathvar$BUB1 <- NA
EntireCohort_pathvar$BUB1[EntireCohort_pathvar$GENE1 %in% c("BUB1")] <- 1
EntireCohort_pathvar$BUB1[EntireCohort_pathvar$GENE2 %in% c("BUB1")] <- 1
EntireCohort_pathvar$BUB1[EntireCohort_pathvar$GENE3 %in% c("BUB1")] <- 1
EntireCohort_pathvar$BUB1[EntireCohort_pathvar$GENE4 %in% c("BUB1")] <- 1
EntireCohort_pathvar$BUB1[EntireCohort_pathvar$GENE5 %in% c("BUB1")] <- 1
EntireCohort_pathvar$BUB1[EntireCohort_pathvar$GENE6 %in% c("BUB1")] <- 1
EntireCohort_pathvar$BUB1[EntireCohort_pathvar$GENE7 %in% c("BUB1")] <- 1
EntireCohort_pathvar$BUB1[EntireCohort_pathvar$GENE8 %in% c("BUB1")] <- 1
EntireCohort_pathvar$BUB1[EntireCohort_pathvar$GENE9 %in% c("BUB1")] <- 1

EntireCohort_pathvar$BUB1[is.na(EntireCohort_pathvar$BUB1)] <- 0
EntireCohort_pathvar$BUB1 <- as.numeric(EntireCohort_pathvar$BUB1)
unique(EntireCohort_pathvar$BUB1)
# 0 1

table(EntireCohort_pathvar$BUB1)
#   0      1 
# 468857    942 




### 'CDK5RAP2', 
EntireCohort_pathvar$CDK5RAP2 <- NA
EntireCohort_pathvar$CDK5RAP2[EntireCohort_pathvar$GENE1 %in% c("CDK5RAP2")] <- 1
EntireCohort_pathvar$CDK5RAP2[EntireCohort_pathvar$GENE2 %in% c("CDK5RAP2")] <- 1
EntireCohort_pathvar$CDK5RAP2[EntireCohort_pathvar$GENE3 %in% c("CDK5RAP2")] <- 1
EntireCohort_pathvar$CDK5RAP2[EntireCohort_pathvar$GENE4 %in% c("CDK5RAP2")] <- 1
EntireCohort_pathvar$CDK5RAP2[EntireCohort_pathvar$GENE5 %in% c("CDK5RAP2")] <- 1
EntireCohort_pathvar$CDK5RAP2[EntireCohort_pathvar$GENE6 %in% c("CDK5RAP2")] <- 1
EntireCohort_pathvar$CDK5RAP2[EntireCohort_pathvar$GENE7 %in% c("CDK5RAP2")] <- 1
EntireCohort_pathvar$CDK5RAP2[EntireCohort_pathvar$GENE8 %in% c("CDK5RAP2")] <- 1
EntireCohort_pathvar$CDK5RAP2[EntireCohort_pathvar$GENE9 %in% c("CDK5RAP2")] <- 1

EntireCohort_pathvar$CDK5RAP2[is.na(EntireCohort_pathvar$CDK5RAP2)] <- 0
EntireCohort_pathvar$CDK5RAP2 <- as.numeric(EntireCohort_pathvar$CDK5RAP2)
unique(EntireCohort_pathvar$CDK5RAP2)
# 0 1

table(EntireCohort_pathvar$CDK5RAP2)
#  0      1 
# 468404   1395




### 'CENPE', 
EntireCohort_pathvar$CENPE <- NA
EntireCohort_pathvar$CENPE[EntireCohort_pathvar$GENE1 %in% c("CENPE")] <- 1
EntireCohort_pathvar$CENPE[EntireCohort_pathvar$GENE2 %in% c("CENPE")] <- 1
EntireCohort_pathvar$CENPE[EntireCohort_pathvar$GENE3 %in% c("CENPE")] <- 1
EntireCohort_pathvar$CENPE[EntireCohort_pathvar$GENE4 %in% c("CENPE")] <- 1
EntireCohort_pathvar$CENPE[EntireCohort_pathvar$GENE5 %in% c("CENPE")] <- 1
EntireCohort_pathvar$CENPE[EntireCohort_pathvar$GENE6 %in% c("CENPE")] <- 1
EntireCohort_pathvar$CENPE[EntireCohort_pathvar$GENE7 %in% c("CENPE")] <- 1
EntireCohort_pathvar$CENPE[EntireCohort_pathvar$GENE8 %in% c("CENPE")] <- 1
EntireCohort_pathvar$CENPE[EntireCohort_pathvar$GENE9 %in% c("CENPE")] <- 1

EntireCohort_pathvar$CENPE[is.na(EntireCohort_pathvar$CENPE)] <- 0
EntireCohort_pathvar$CENPE <- as.numeric(EntireCohort_pathvar$CENPE)
unique(EntireCohort_pathvar$CENPE)
# 0 1

table(EntireCohort_pathvar$CENPE)
#  0      1 
# 469004    795




### 'CENPF', 
EntireCohort_pathvar$CENPF <- NA
EntireCohort_pathvar$CENPF[EntireCohort_pathvar$GENE1 %in% c("CENPF")] <- 1
EntireCohort_pathvar$CENPF[EntireCohort_pathvar$GENE2 %in% c("CENPF")] <- 1
EntireCohort_pathvar$CENPF[EntireCohort_pathvar$GENE3 %in% c("CENPF")] <- 1
EntireCohort_pathvar$CENPF[EntireCohort_pathvar$GENE4 %in% c("CENPF")] <- 1
EntireCohort_pathvar$CENPF[EntireCohort_pathvar$GENE5 %in% c("CENPF")] <- 1
EntireCohort_pathvar$CENPF[EntireCohort_pathvar$GENE6 %in% c("CENPF")] <- 1
EntireCohort_pathvar$CENPF[EntireCohort_pathvar$GENE7 %in% c("CENPF")] <- 1
EntireCohort_pathvar$CENPF[EntireCohort_pathvar$GENE8 %in% c("CENPF")] <- 1
EntireCohort_pathvar$CENPF[EntireCohort_pathvar$GENE9 %in% c("CENPF")] <- 1

EntireCohort_pathvar$CENPF[is.na(EntireCohort_pathvar$CENPF)] <- 0
EntireCohort_pathvar$CENPF <- as.numeric(EntireCohort_pathvar$CENPF)
unique(EntireCohort_pathvar$CENPF)
# 0 1

table(EntireCohort_pathvar$CENPF)
#  0      1 
# 467852   1947




### 'CENPJ', 
EntireCohort_pathvar$CENPJ <- NA
EntireCohort_pathvar$CENPJ[EntireCohort_pathvar$GENE1 %in% c("CENPJ")] <- 1
EntireCohort_pathvar$CENPJ[EntireCohort_pathvar$GENE2 %in% c("CENPJ")] <- 1
EntireCohort_pathvar$CENPJ[EntireCohort_pathvar$GENE3 %in% c("CENPJ")] <- 1
EntireCohort_pathvar$CENPJ[EntireCohort_pathvar$GENE4 %in% c("CENPJ")] <- 1
EntireCohort_pathvar$CENPJ[EntireCohort_pathvar$GENE5 %in% c("CENPJ")] <- 1
EntireCohort_pathvar$CENPJ[EntireCohort_pathvar$GENE6 %in% c("CENPJ")] <- 1
EntireCohort_pathvar$CENPJ[EntireCohort_pathvar$GENE7 %in% c("CENPJ")] <- 1
EntireCohort_pathvar$CENPJ[EntireCohort_pathvar$GENE8 %in% c("CENPJ")] <- 1
EntireCohort_pathvar$CENPJ[EntireCohort_pathvar$GENE9 %in% c("CENPJ")] <- 1

EntireCohort_pathvar$CENPJ[is.na(EntireCohort_pathvar$CENPJ)] <- 0
EntireCohort_pathvar$CENPJ <- as.numeric(EntireCohort_pathvar$CENPJ)
unique(EntireCohort_pathvar$CENPJ)
# 0 1

table(EntireCohort_pathvar$CENPJ)
#  0      1 
# 467484   2315





### 'CENPT', 
EntireCohort_pathvar$CENPT <- NA
EntireCohort_pathvar$CENPT[EntireCohort_pathvar$GENE1 %in% c("CENPT")] <- 1
EntireCohort_pathvar$CENPT[EntireCohort_pathvar$GENE2 %in% c("CENPT")] <- 1
EntireCohort_pathvar$CENPT[EntireCohort_pathvar$GENE3 %in% c("CENPT")] <- 1
EntireCohort_pathvar$CENPT[EntireCohort_pathvar$GENE4 %in% c("CENPT")] <- 1
EntireCohort_pathvar$CENPT[EntireCohort_pathvar$GENE5 %in% c("CENPT")] <- 1
EntireCohort_pathvar$CENPT[EntireCohort_pathvar$GENE6 %in% c("CENPT")] <- 1
EntireCohort_pathvar$CENPT[EntireCohort_pathvar$GENE7 %in% c("CENPT")] <- 1
EntireCohort_pathvar$CENPT[EntireCohort_pathvar$GENE8 %in% c("CENPT")] <- 1
EntireCohort_pathvar$CENPT[EntireCohort_pathvar$GENE9 %in% c("CENPT")] <- 1

EntireCohort_pathvar$CENPT[is.na(EntireCohort_pathvar$CENPT)] <- 0
EntireCohort_pathvar$CENPT <- as.numeric(EntireCohort_pathvar$CENPT)
unique(EntireCohort_pathvar$CENPT)
# 0 1

table(EntireCohort_pathvar$CENPT)
#  0      1 
# 466837   2962




### 'CEP135', 
EntireCohort_pathvar$CEP135 <- NA
EntireCohort_pathvar$CEP135[EntireCohort_pathvar$GENE1 %in% c("CEP135")] <- 1
EntireCohort_pathvar$CEP135[EntireCohort_pathvar$GENE2 %in% c("CEP135")] <- 1
EntireCohort_pathvar$CEP135[EntireCohort_pathvar$GENE3 %in% c("CEP135")] <- 1
EntireCohort_pathvar$CEP135[EntireCohort_pathvar$GENE4 %in% c("CEP135")] <- 1
EntireCohort_pathvar$CEP135[EntireCohort_pathvar$GENE5 %in% c("CEP135")] <- 1
EntireCohort_pathvar$CEP135[EntireCohort_pathvar$GENE6 %in% c("CEP135")] <- 1
EntireCohort_pathvar$CEP135[EntireCohort_pathvar$GENE7 %in% c("CEP135")] <- 1
EntireCohort_pathvar$CEP135[EntireCohort_pathvar$GENE8 %in% c("CEP135")] <- 1
EntireCohort_pathvar$CEP135[EntireCohort_pathvar$GENE9 %in% c("CEP135")] <- 1

EntireCohort_pathvar$CEP135[is.na(EntireCohort_pathvar$CEP135)] <- 0
EntireCohort_pathvar$CEP135 <- as.numeric(EntireCohort_pathvar$CEP135)
unique(EntireCohort_pathvar$CEP135)
# 0 1

table(EntireCohort_pathvar$CEP135)
#    0      1 
# 467023   2776




### 'CEP152', 
EntireCohort_pathvar$CEP152 <- NA
EntireCohort_pathvar$CEP152[EntireCohort_pathvar$GENE1 %in% c("CEP152")] <- 1
EntireCohort_pathvar$CEP152[EntireCohort_pathvar$GENE2 %in% c("CEP152")] <- 1
EntireCohort_pathvar$CEP152[EntireCohort_pathvar$GENE3 %in% c("CEP152")] <- 1
EntireCohort_pathvar$CEP152[EntireCohort_pathvar$GENE4 %in% c("CEP152")] <- 1
EntireCohort_pathvar$CEP152[EntireCohort_pathvar$GENE5 %in% c("CEP152")] <- 1
EntireCohort_pathvar$CEP152[EntireCohort_pathvar$GENE6 %in% c("CEP152")] <- 1
EntireCohort_pathvar$CEP152[EntireCohort_pathvar$GENE7 %in% c("CEP152")] <- 1
EntireCohort_pathvar$CEP152[EntireCohort_pathvar$GENE8 %in% c("CEP152")] <- 1
EntireCohort_pathvar$CEP152[EntireCohort_pathvar$GENE9 %in% c("CEP152")] <- 1

EntireCohort_pathvar$CEP152[is.na(EntireCohort_pathvar$CEP152)] <- 0
EntireCohort_pathvar$CEP152 <- as.numeric(EntireCohort_pathvar$CEP152)
unique(EntireCohort_pathvar$CEP152)
# 0 1

table(EntireCohort_pathvar$CEP152)
#  0      1 
# 467630   2169




### 'CEP63', 
EntireCohort_pathvar$CEP63 <- NA
EntireCohort_pathvar$CEP63[EntireCohort_pathvar$GENE1 %in% c("CEP63")] <- 1
EntireCohort_pathvar$CEP63[EntireCohort_pathvar$GENE2 %in% c("CEP63")] <- 1
EntireCohort_pathvar$CEP63[EntireCohort_pathvar$GENE3 %in% c("CEP63")] <- 1
EntireCohort_pathvar$CEP63[EntireCohort_pathvar$GENE4 %in% c("CEP63")] <- 1
EntireCohort_pathvar$CEP63[EntireCohort_pathvar$GENE5 %in% c("CEP63")] <- 1
EntireCohort_pathvar$CEP63[EntireCohort_pathvar$GENE6 %in% c("CEP63")] <- 1
EntireCohort_pathvar$CEP63[EntireCohort_pathvar$GENE7 %in% c("CEP63")] <- 1
EntireCohort_pathvar$CEP63[EntireCohort_pathvar$GENE8 %in% c("CEP63")] <- 1
EntireCohort_pathvar$CEP63[EntireCohort_pathvar$GENE9 %in% c("CEP63")] <- 1

EntireCohort_pathvar$CEP63[is.na(EntireCohort_pathvar$CEP63)] <- 0
EntireCohort_pathvar$CEP63 <- as.numeric(EntireCohort_pathvar$CEP63)
unique(EntireCohort_pathvar$CEP63)
# 0 1

table(EntireCohort_pathvar$CEP63)
#    0      1 
# 468779   1020




### 'CIT', 
EntireCohort_pathvar$CIT <- NA
EntireCohort_pathvar$CIT[EntireCohort_pathvar$GENE1 %in% c("CIT")] <- 1
EntireCohort_pathvar$CIT[EntireCohort_pathvar$GENE2 %in% c("CIT")] <- 1
EntireCohort_pathvar$CIT[EntireCohort_pathvar$GENE3 %in% c("CIT")] <- 1
EntireCohort_pathvar$CIT[EntireCohort_pathvar$GENE4 %in% c("CIT")] <- 1
EntireCohort_pathvar$CIT[EntireCohort_pathvar$GENE5 %in% c("CIT")] <- 1
EntireCohort_pathvar$CIT[EntireCohort_pathvar$GENE6 %in% c("CIT")] <- 1
EntireCohort_pathvar$CIT[EntireCohort_pathvar$GENE7 %in% c("CIT")] <- 1
EntireCohort_pathvar$CIT[EntireCohort_pathvar$GENE8 %in% c("CIT")] <- 1
EntireCohort_pathvar$CIT[EntireCohort_pathvar$GENE9 %in% c("CIT")] <- 1

EntireCohort_pathvar$CIT[is.na(EntireCohort_pathvar$CIT)] <- 0
EntireCohort_pathvar$CIT <- as.numeric(EntireCohort_pathvar$CIT)
unique(EntireCohort_pathvar$CIT)
# 0 1

table(EntireCohort_pathvar$CIT)
#  0      1 
# 466010   3789




### 'CKAP2L', 
EntireCohort_pathvar$CKAP2L <- NA
EntireCohort_pathvar$CKAP2L[EntireCohort_pathvar$GENE1 %in% c("CKAP2L")] <- 1
EntireCohort_pathvar$CKAP2L[EntireCohort_pathvar$GENE2 %in% c("CKAP2L")] <- 1
EntireCohort_pathvar$CKAP2L[EntireCohort_pathvar$GENE3 %in% c("CKAP2L")] <- 1
EntireCohort_pathvar$CKAP2L[EntireCohort_pathvar$GENE4 %in% c("CKAP2L")] <- 1
EntireCohort_pathvar$CKAP2L[EntireCohort_pathvar$GENE5 %in% c("CKAP2L")] <- 1
EntireCohort_pathvar$CKAP2L[EntireCohort_pathvar$GENE6 %in% c("CKAP2L")] <- 1
EntireCohort_pathvar$CKAP2L[EntireCohort_pathvar$GENE7 %in% c("CKAP2L")] <- 1
EntireCohort_pathvar$CKAP2L[EntireCohort_pathvar$GENE8 %in% c("CKAP2L")] <- 1
EntireCohort_pathvar$CKAP2L[EntireCohort_pathvar$GENE9 %in% c("CKAP2L")] <- 1

EntireCohort_pathvar$CKAP2L[is.na(EntireCohort_pathvar$CKAP2L)] <- 0
EntireCohort_pathvar$CKAP2L <- as.numeric(EntireCohort_pathvar$CKAP2L)
unique(EntireCohort_pathvar$CKAP2L)
# 0 1

table(EntireCohort_pathvar$CKAP2L)
#  0      1 
# 469225    574 





### 'KATNB1', 
EntireCohort_pathvar$KATNB1 <- NA
EntireCohort_pathvar$KATNB1[EntireCohort_pathvar$GENE1 %in% c("KATNB1")] <- 1
EntireCohort_pathvar$KATNB1[EntireCohort_pathvar$GENE2 %in% c("KATNB1")] <- 1
EntireCohort_pathvar$KATNB1[EntireCohort_pathvar$GENE3 %in% c("KATNB1")] <- 1
EntireCohort_pathvar$KATNB1[EntireCohort_pathvar$GENE4 %in% c("KATNB1")] <- 1
EntireCohort_pathvar$KATNB1[EntireCohort_pathvar$GENE5 %in% c("KATNB1")] <- 1
EntireCohort_pathvar$KATNB1[EntireCohort_pathvar$GENE6 %in% c("KATNB1")] <- 1
EntireCohort_pathvar$KATNB1[EntireCohort_pathvar$GENE7 %in% c("KATNB1")] <- 1
EntireCohort_pathvar$KATNB1[EntireCohort_pathvar$GENE8 %in% c("KATNB1")] <- 1
EntireCohort_pathvar$KATNB1[EntireCohort_pathvar$GENE9 %in% c("KATNB1")] <- 1

EntireCohort_pathvar$KATNB1[is.na(EntireCohort_pathvar$KATNB1)] <- 0
EntireCohort_pathvar$KATNB1 <- as.numeric(EntireCohort_pathvar$KATNB1)
unique(EntireCohort_pathvar$KATNB1)
# 0 1

table(EntireCohort_pathvar$KATNB1)
#   0      1 
# 468689   1110





### 'KIF11', 
EntireCohort_pathvar$KIF11 <- NA
EntireCohort_pathvar$KIF11[EntireCohort_pathvar$GENE1 %in% c("KIF11")] <- 1
EntireCohort_pathvar$KIF11[EntireCohort_pathvar$GENE2 %in% c("KIF11")] <- 1
EntireCohort_pathvar$KIF11[EntireCohort_pathvar$GENE3 %in% c("KIF11")] <- 1
EntireCohort_pathvar$KIF11[EntireCohort_pathvar$GENE4 %in% c("KIF11")] <- 1
EntireCohort_pathvar$KIF11[EntireCohort_pathvar$GENE5 %in% c("KIF11")] <- 1
EntireCohort_pathvar$KIF11[EntireCohort_pathvar$GENE6 %in% c("KIF11")] <- 1
EntireCohort_pathvar$KIF11[EntireCohort_pathvar$GENE7 %in% c("KIF11")] <- 1
EntireCohort_pathvar$KIF11[EntireCohort_pathvar$GENE8 %in% c("KIF11")] <- 1
EntireCohort_pathvar$KIF11[EntireCohort_pathvar$GENE9 %in% c("KIF11")] <- 1

EntireCohort_pathvar$KIF11[is.na(EntireCohort_pathvar$KIF11)] <- 0
EntireCohort_pathvar$KIF11 <- as.numeric(EntireCohort_pathvar$KIF11)
unique(EntireCohort_pathvar$KIF11)
# 0 1

table(EntireCohort_pathvar$KIF11)
#  0      1 
# 469287    512 




### 'KIF22', 
# no variants present




### 'KIFBP', 
EntireCohort_pathvar$KIFBP <- NA
EntireCohort_pathvar$KIFBP[EntireCohort_pathvar$GENE1 %in% c("KIFBP")] <- 1
EntireCohort_pathvar$KIFBP[EntireCohort_pathvar$GENE2 %in% c("KIFBP")] <- 1
EntireCohort_pathvar$KIFBP[EntireCohort_pathvar$GENE3 %in% c("KIFBP")] <- 1
EntireCohort_pathvar$KIFBP[EntireCohort_pathvar$GENE4 %in% c("KIFBP")] <- 1
EntireCohort_pathvar$KIFBP[EntireCohort_pathvar$GENE5 %in% c("KIFBP")] <- 1
EntireCohort_pathvar$KIFBP[EntireCohort_pathvar$GENE6 %in% c("KIFBP")] <- 1
EntireCohort_pathvar$KIFBP[EntireCohort_pathvar$GENE7 %in% c("KIFBP")] <- 1
EntireCohort_pathvar$KIFBP[EntireCohort_pathvar$GENE8 %in% c("KIFBP")] <- 1
EntireCohort_pathvar$KIFBP[EntireCohort_pathvar$GENE9 %in% c("KIFBP")] <- 1

EntireCohort_pathvar$KIFBP[is.na(EntireCohort_pathvar$KIFBP)] <- 0
EntireCohort_pathvar$KIFBP <- as.numeric(EntireCohort_pathvar$KIFBP)
unique(EntireCohort_pathvar$KIFBP)
# 0 1

table(EntireCohort_pathvar$KIFBP)
# 0      1 
# 468843    956




### 'KNL1', 
EntireCohort_pathvar$KNL1 <- NA
EntireCohort_pathvar$KNL1[EntireCohort_pathvar$GENE1 %in% c("KNL1")] <- 1
EntireCohort_pathvar$KNL1[EntireCohort_pathvar$GENE2 %in% c("KNL1")] <- 1
EntireCohort_pathvar$KNL1[EntireCohort_pathvar$GENE3 %in% c("KNL1")] <- 1
EntireCohort_pathvar$KNL1[EntireCohort_pathvar$GENE4 %in% c("KNL1")] <- 1
EntireCohort_pathvar$KNL1[EntireCohort_pathvar$GENE5 %in% c("KNL1")] <- 1
EntireCohort_pathvar$KNL1[EntireCohort_pathvar$GENE6 %in% c("KNL1")] <- 1
EntireCohort_pathvar$KNL1[EntireCohort_pathvar$GENE7 %in% c("KNL1")] <- 1
EntireCohort_pathvar$KNL1[EntireCohort_pathvar$GENE8 %in% c("KNL1")] <- 1
EntireCohort_pathvar$KNL1[EntireCohort_pathvar$GENE9 %in% c("KNL1")] <- 1

EntireCohort_pathvar$KNL1[is.na(EntireCohort_pathvar$KNL1)] <- 0
EntireCohort_pathvar$KNL1 <- as.numeric(EntireCohort_pathvar$KNL1)
unique(EntireCohort_pathvar$KNL1)
# 0 1

table(EntireCohort_pathvar$KNL1)
#  0      1 
# 469353    446




### 'NDE1', 
EntireCohort_pathvar$NDE1 <- NA
EntireCohort_pathvar$NDE1[EntireCohort_pathvar$GENE1 %in% c("NDE1")] <- 1
EntireCohort_pathvar$NDE1[EntireCohort_pathvar$GENE2 %in% c("NDE1")] <- 1
EntireCohort_pathvar$NDE1[EntireCohort_pathvar$GENE3 %in% c("NDE1")] <- 1
EntireCohort_pathvar$NDE1[EntireCohort_pathvar$GENE4 %in% c("NDE1")] <- 1
EntireCohort_pathvar$NDE1[EntireCohort_pathvar$GENE5 %in% c("NDE1")] <- 1
EntireCohort_pathvar$NDE1[EntireCohort_pathvar$GENE6 %in% c("NDE1")] <- 1
EntireCohort_pathvar$NDE1[EntireCohort_pathvar$GENE7 %in% c("NDE1")] <- 1
EntireCohort_pathvar$NDE1[EntireCohort_pathvar$GENE8 %in% c("NDE1")] <- 1
EntireCohort_pathvar$NDE1[EntireCohort_pathvar$GENE9 %in% c("NDE1")] <- 1

EntireCohort_pathvar$NDE1[is.na(EntireCohort_pathvar$NDE1)] <- 0
EntireCohort_pathvar$NDE1 <- as.numeric(EntireCohort_pathvar$NDE1)
unique(EntireCohort_pathvar$NDE1)
# 0 1

table(EntireCohort_pathvar$NDE1)
#  0      1 
# 469363    436




### 'NIN', 
EntireCohort_pathvar$NIN <- NA
EntireCohort_pathvar$NIN[EntireCohort_pathvar$GENE1 %in% c("NIN")] <- 1
EntireCohort_pathvar$NIN[EntireCohort_pathvar$GENE2 %in% c("NIN")] <- 1
EntireCohort_pathvar$NIN[EntireCohort_pathvar$GENE3 %in% c("NIN")] <- 1
EntireCohort_pathvar$NIN[EntireCohort_pathvar$GENE4 %in% c("NIN")] <- 1
EntireCohort_pathvar$NIN[EntireCohort_pathvar$GENE5 %in% c("NIN")] <- 1
EntireCohort_pathvar$NIN[EntireCohort_pathvar$GENE6 %in% c("NIN")] <- 1
EntireCohort_pathvar$NIN[EntireCohort_pathvar$GENE7 %in% c("NIN")] <- 1
EntireCohort_pathvar$NIN[EntireCohort_pathvar$GENE8 %in% c("NIN")] <- 1
EntireCohort_pathvar$NIN[EntireCohort_pathvar$GENE9 %in% c("NIN")] <- 1

EntireCohort_pathvar$NIN[is.na(EntireCohort_pathvar$NIN)] <- 0
EntireCohort_pathvar$NIN <- as.numeric(EntireCohort_pathvar$NIN)
unique(EntireCohort_pathvar$NIN)
# 0 1

table(EntireCohort_pathvar$NIN)
#  0      1 
# 468112   1687






### 'PCNT', 
EntireCohort_pathvar$PCNT <- NA
EntireCohort_pathvar$PCNT[EntireCohort_pathvar$GENE1 %in% c("PCNT")] <- 1
EntireCohort_pathvar$PCNT[EntireCohort_pathvar$GENE2 %in% c("PCNT")] <- 1
EntireCohort_pathvar$PCNT[EntireCohort_pathvar$GENE3 %in% c("PCNT")] <- 1
EntireCohort_pathvar$PCNT[EntireCohort_pathvar$GENE4 %in% c("PCNT")] <- 1
EntireCohort_pathvar$PCNT[EntireCohort_pathvar$GENE5 %in% c("PCNT")] <- 1
EntireCohort_pathvar$PCNT[EntireCohort_pathvar$GENE6 %in% c("PCNT")] <- 1
EntireCohort_pathvar$PCNT[EntireCohort_pathvar$GENE7 %in% c("PCNT")] <- 1
EntireCohort_pathvar$PCNT[EntireCohort_pathvar$GENE8 %in% c("PCNT")] <- 1
EntireCohort_pathvar$PCNT[EntireCohort_pathvar$GENE9 %in% c("PCNT")] <- 1

EntireCohort_pathvar$PCNT[is.na(EntireCohort_pathvar$PCNT)] <- 0
EntireCohort_pathvar$PCNT <- as.numeric(EntireCohort_pathvar$PCNT)
unique(EntireCohort_pathvar$PCNT)
# 0 1

table(EntireCohort_pathvar$PCNT)
#  0      1 
# 468225   1574





### 'PLK4', 
EntireCohort_pathvar$PLK4 <- NA
EntireCohort_pathvar$PLK4[EntireCohort_pathvar$GENE1 %in% c("PLK4")] <- 1
EntireCohort_pathvar$PLK4[EntireCohort_pathvar$GENE2 %in% c("PLK4")] <- 1
EntireCohort_pathvar$PLK4[EntireCohort_pathvar$GENE3 %in% c("PLK4")] <- 1
EntireCohort_pathvar$PLK4[EntireCohort_pathvar$GENE4 %in% c("PLK4")] <- 1
EntireCohort_pathvar$PLK4[EntireCohort_pathvar$GENE5 %in% c("PLK4")] <- 1
EntireCohort_pathvar$PLK4[EntireCohort_pathvar$GENE6 %in% c("PLK4")] <- 1
EntireCohort_pathvar$PLK4[EntireCohort_pathvar$GENE7 %in% c("PLK4")] <- 1
EntireCohort_pathvar$PLK4[EntireCohort_pathvar$GENE8 %in% c("PLK4")] <- 1
EntireCohort_pathvar$PLK4[EntireCohort_pathvar$GENE9 %in% c("PLK4")] <- 1

EntireCohort_pathvar$PLK4[is.na(EntireCohort_pathvar$PLK4)] <- 0
EntireCohort_pathvar$PLK4 <- as.numeric(EntireCohort_pathvar$PLK4)
unique(EntireCohort_pathvar$PLK4)
# 0 1

table(EntireCohort_pathvar$PLK4)
#   0      1 
# 468907    892




### 'RTTN', 
EntireCohort_pathvar$RTTN <- NA
EntireCohort_pathvar$RTTN[EntireCohort_pathvar$GENE1 %in% c("RTTN")] <- 1
EntireCohort_pathvar$RTTN[EntireCohort_pathvar$GENE2 %in% c("RTTN")] <- 1
EntireCohort_pathvar$RTTN[EntireCohort_pathvar$GENE3 %in% c("RTTN")] <- 1
EntireCohort_pathvar$RTTN[EntireCohort_pathvar$GENE4 %in% c("RTTN")] <- 1
EntireCohort_pathvar$RTTN[EntireCohort_pathvar$GENE5 %in% c("RTTN")] <- 1
EntireCohort_pathvar$RTTN[EntireCohort_pathvar$GENE6 %in% c("RTTN")] <- 1
EntireCohort_pathvar$RTTN[EntireCohort_pathvar$GENE7 %in% c("RTTN")] <- 1
EntireCohort_pathvar$RTTN[EntireCohort_pathvar$GENE8 %in% c("RTTN")] <- 1
EntireCohort_pathvar$RTTN[EntireCohort_pathvar$GENE9 %in% c("RTTN")] <- 1

EntireCohort_pathvar$RTTN[is.na(EntireCohort_pathvar$RTTN)] <- 0
EntireCohort_pathvar$RTTN <- as.numeric(EntireCohort_pathvar$RTTN)
unique(EntireCohort_pathvar$RTTN)
# 0 1

table(EntireCohort_pathvar$RTTN)
#   0      1 
# 467184   2615




### 'SASS6',        
EntireCohort_pathvar$SASS6 <- NA
EntireCohort_pathvar$SASS6[EntireCohort_pathvar$GENE1 %in% c("SASS6")] <- 1
EntireCohort_pathvar$SASS6[EntireCohort_pathvar$GENE2 %in% c("SASS6")] <- 1
EntireCohort_pathvar$SASS6[EntireCohort_pathvar$GENE3 %in% c("SASS6")] <- 1
EntireCohort_pathvar$SASS6[EntireCohort_pathvar$GENE4 %in% c("SASS6")] <- 1
EntireCohort_pathvar$SASS6[EntireCohort_pathvar$GENE5 %in% c("SASS6")] <- 1
EntireCohort_pathvar$SASS6[EntireCohort_pathvar$GENE6 %in% c("SASS6")] <- 1
EntireCohort_pathvar$SASS6[EntireCohort_pathvar$GENE7 %in% c("SASS6")] <- 1
EntireCohort_pathvar$SASS6[EntireCohort_pathvar$GENE8 %in% c("SASS6")] <- 1
EntireCohort_pathvar$SASS6[EntireCohort_pathvar$GENE9 %in% c("SASS6")] <- 1

EntireCohort_pathvar$SASS6[is.na(EntireCohort_pathvar$SASS6)] <- 0
EntireCohort_pathvar$SASS6 <- as.numeric(EntireCohort_pathvar$SASS6)
unique(EntireCohort_pathvar$SASS6)
# 0 1

table(EntireCohort_pathvar$SASS6)
#  0      1 
# 469060    739




### 'STIL',  
EntireCohort_pathvar$STIL <- NA
EntireCohort_pathvar$STIL[EntireCohort_pathvar$GENE1 %in% c("STIL")] <- 1
EntireCohort_pathvar$STIL[EntireCohort_pathvar$GENE2 %in% c("STIL")] <- 1
EntireCohort_pathvar$STIL[EntireCohort_pathvar$GENE3 %in% c("STIL")] <- 1
EntireCohort_pathvar$STIL[EntireCohort_pathvar$GENE4 %in% c("STIL")] <- 1
EntireCohort_pathvar$STIL[EntireCohort_pathvar$GENE5 %in% c("STIL")] <- 1
EntireCohort_pathvar$STIL[EntireCohort_pathvar$GENE6 %in% c("STIL")] <- 1
EntireCohort_pathvar$STIL[EntireCohort_pathvar$GENE7 %in% c("STIL")] <- 1
EntireCohort_pathvar$STIL[EntireCohort_pathvar$GENE8 %in% c("STIL")] <- 1
EntireCohort_pathvar$STIL[EntireCohort_pathvar$GENE9 %in% c("STIL")] <- 1

EntireCohort_pathvar$STIL[is.na(EntireCohort_pathvar$STIL)] <- 0
EntireCohort_pathvar$STIL <- as.numeric(EntireCohort_pathvar$STIL)
unique(EntireCohort_pathvar$STIL)
# 0 1

table(EntireCohort_pathvar$STIL)
#  0      1 
# 467747   2052




### 'TRAPPC14', 
EntireCohort_pathvar$TRAPPC14 <- NA
EntireCohort_pathvar$TRAPPC14[EntireCohort_pathvar$GENE1 %in% c("TRAPPC14")] <- 1
EntireCohort_pathvar$TRAPPC14[EntireCohort_pathvar$GENE2 %in% c("TRAPPC14")] <- 1
EntireCohort_pathvar$TRAPPC14[EntireCohort_pathvar$GENE3 %in% c("TRAPPC14")] <- 1
EntireCohort_pathvar$TRAPPC14[EntireCohort_pathvar$GENE4 %in% c("TRAPPC14")] <- 1
EntireCohort_pathvar$TRAPPC14[EntireCohort_pathvar$GENE5 %in% c("TRAPPC14")] <- 1
EntireCohort_pathvar$TRAPPC14[EntireCohort_pathvar$GENE6 %in% c("TRAPPC14")] <- 1
EntireCohort_pathvar$TRAPPC14[EntireCohort_pathvar$GENE7 %in% c("TRAPPC14")] <- 1
EntireCohort_pathvar$TRAPPC14[EntireCohort_pathvar$GENE8 %in% c("TRAPPC14")] <- 1
EntireCohort_pathvar$TRAPPC14[EntireCohort_pathvar$GENE9 %in% c("TRAPPC14")] <- 1

EntireCohort_pathvar$TRAPPC14[is.na(EntireCohort_pathvar$TRAPPC14)] <- 0
EntireCohort_pathvar$TRAPPC14 <- as.numeric(EntireCohort_pathvar$TRAPPC14)
unique(EntireCohort_pathvar$TRAPPC14)
# 0 1

table(EntireCohort_pathvar$TRAPPC14)
#  0      1 
# 469239    560




### 'TUBGCP6',  
EntireCohort_pathvar$TUBGCP6 <- NA
EntireCohort_pathvar$TUBGCP6[EntireCohort_pathvar$GENE1 %in% c("TUBGCP6")] <- 1
EntireCohort_pathvar$TUBGCP6[EntireCohort_pathvar$GENE2 %in% c("TUBGCP6")] <- 1
EntireCohort_pathvar$TUBGCP6[EntireCohort_pathvar$GENE3 %in% c("TUBGCP6")] <- 1
EntireCohort_pathvar$TUBGCP6[EntireCohort_pathvar$GENE4 %in% c("TUBGCP6")] <- 1
EntireCohort_pathvar$TUBGCP6[EntireCohort_pathvar$GENE5 %in% c("TUBGCP6")] <- 1
EntireCohort_pathvar$TUBGCP6[EntireCohort_pathvar$GENE6 %in% c("TUBGCP6")] <- 1
EntireCohort_pathvar$TUBGCP6[EntireCohort_pathvar$GENE7 %in% c("TUBGCP6")] <- 1
EntireCohort_pathvar$TUBGCP6[EntireCohort_pathvar$GENE8 %in% c("TUBGCP6")] <- 1
EntireCohort_pathvar$TUBGCP6[EntireCohort_pathvar$GENE9 %in% c("TUBGCP6")] <- 1

EntireCohort_pathvar$TUBGCP6[is.na(EntireCohort_pathvar$TUBGCP6)] <- 0
EntireCohort_pathvar$TUBGCP6 <- as.numeric(EntireCohort_pathvar$TUBGCP6)
unique(EntireCohort_pathvar$TUBGCP6)
# 0 1

table(EntireCohort_pathvar$TUBGCP6)
#   0      1 
# 468262   1537 




### 'WDR62', 
EntireCohort_pathvar$WDR62 <- NA
EntireCohort_pathvar$WDR62[EntireCohort_pathvar$GENE1 %in% c("WDR62")] <- 1
EntireCohort_pathvar$WDR62[EntireCohort_pathvar$GENE2 %in% c("WDR62")] <- 1
EntireCohort_pathvar$WDR62[EntireCohort_pathvar$GENE3 %in% c("WDR62")] <- 1
EntireCohort_pathvar$WDR62[EntireCohort_pathvar$GENE4 %in% c("WDR62")] <- 1
EntireCohort_pathvar$WDR62[EntireCohort_pathvar$GENE5 %in% c("WDR62")] <- 1
EntireCohort_pathvar$WDR62[EntireCohort_pathvar$GENE6 %in% c("WDR62")] <- 1
EntireCohort_pathvar$WDR62[EntireCohort_pathvar$GENE7 %in% c("WDR62")] <- 1
EntireCohort_pathvar$WDR62[EntireCohort_pathvar$GENE8 %in% c("WDR62")] <- 1
EntireCohort_pathvar$WDR62[EntireCohort_pathvar$GENE9 %in% c("WDR62")] <- 1

EntireCohort_pathvar$WDR62[is.na(EntireCohort_pathvar$WDR62)] <- 0
EntireCohort_pathvar$WDR62 <- as.numeric(EntireCohort_pathvar$WDR62)
unique(EntireCohort_pathvar$WDR62)
# 0 1

table(EntireCohort_pathvar$WDR62)
#  0      1 
# 467948   1851




### 'NEK1', 
EntireCohort_pathvar$NEK1 <- NA
EntireCohort_pathvar$NEK1[EntireCohort_pathvar$GENE1 %in% c("NEK1")] <- 1
EntireCohort_pathvar$NEK1[EntireCohort_pathvar$GENE2 %in% c("NEK1")] <- 1
EntireCohort_pathvar$NEK1[EntireCohort_pathvar$GENE3 %in% c("NEK1")] <- 1
EntireCohort_pathvar$NEK1[EntireCohort_pathvar$GENE4 %in% c("NEK1")] <- 1
EntireCohort_pathvar$NEK1[EntireCohort_pathvar$GENE5 %in% c("NEK1")] <- 1
EntireCohort_pathvar$NEK1[EntireCohort_pathvar$GENE6 %in% c("NEK1")] <- 1
EntireCohort_pathvar$NEK1[EntireCohort_pathvar$GENE7 %in% c("NEK1")] <- 1
EntireCohort_pathvar$NEK1[EntireCohort_pathvar$GENE8 %in% c("NEK1")] <- 1
EntireCohort_pathvar$NEK1[EntireCohort_pathvar$GENE9 %in% c("NEK1")] <- 1

EntireCohort_pathvar$NEK1[is.na(EntireCohort_pathvar$NEK1)] <- 0
EntireCohort_pathvar$NEK1 <- as.numeric(EntireCohort_pathvar$NEK1)
unique(EntireCohort_pathvar$NEK1)
# 0 1

table(EntireCohort_pathvar$NEK1)
#  0      1 
# 467572   2227





### 'CFAP410', 
EntireCohort_pathvar$CFAP410 <- NA
EntireCohort_pathvar$CFAP410[EntireCohort_pathvar$GENE1 %in% c("CFAP410")] <- 1
EntireCohort_pathvar$CFAP410[EntireCohort_pathvar$GENE2 %in% c("CFAP410")] <- 1
EntireCohort_pathvar$CFAP410[EntireCohort_pathvar$GENE3 %in% c("CFAP410")] <- 1
EntireCohort_pathvar$CFAP410[EntireCohort_pathvar$GENE4 %in% c("CFAP410")] <- 1
EntireCohort_pathvar$CFAP410[EntireCohort_pathvar$GENE5 %in% c("CFAP410")] <- 1
EntireCohort_pathvar$CFAP410[EntireCohort_pathvar$GENE6 %in% c("CFAP410")] <- 1
EntireCohort_pathvar$CFAP410[EntireCohort_pathvar$GENE7 %in% c("CFAP410")] <- 1
EntireCohort_pathvar$CFAP410[EntireCohort_pathvar$GENE8 %in% c("CFAP410")] <- 1
EntireCohort_pathvar$CFAP410[EntireCohort_pathvar$GENE9 %in% c("CFAP410")] <- 1

EntireCohort_pathvar$CFAP410[is.na(EntireCohort_pathvar$CFAP410)] <- 0
EntireCohort_pathvar$CFAP410 <- as.numeric(EntireCohort_pathvar$CFAP410)
unique(EntireCohort_pathvar$CFAP410)
# 0 1

table(EntireCohort_pathvar$CFAP410)
#  0      1 
# 469297    502




### 'LMNA', 
# no variants present




### 'LMNB1', 
# no variants present




### 'LMNB2', 
# no variants present




### 'ZMPSTE24',  
EntireCohort_pathvar$ZMPSTE24 <- NA
EntireCohort_pathvar$ZMPSTE24[EntireCohort_pathvar$GENE1 %in% c("ZMPSTE24")] <- 1
EntireCohort_pathvar$ZMPSTE24[EntireCohort_pathvar$GENE2 %in% c("ZMPSTE24")] <- 1
EntireCohort_pathvar$ZMPSTE24[EntireCohort_pathvar$GENE3 %in% c("ZMPSTE24")] <- 1
EntireCohort_pathvar$ZMPSTE24[EntireCohort_pathvar$GENE4 %in% c("ZMPSTE24")] <- 1
EntireCohort_pathvar$ZMPSTE24[EntireCohort_pathvar$GENE5 %in% c("ZMPSTE24")] <- 1
EntireCohort_pathvar$ZMPSTE24[EntireCohort_pathvar$GENE6 %in% c("ZMPSTE24")] <- 1
EntireCohort_pathvar$ZMPSTE24[EntireCohort_pathvar$GENE7 %in% c("ZMPSTE24")] <- 1
EntireCohort_pathvar$ZMPSTE24[EntireCohort_pathvar$GENE8 %in% c("ZMPSTE24")] <- 1
EntireCohort_pathvar$ZMPSTE24[EntireCohort_pathvar$GENE9 %in% c("ZMPSTE24")] <- 1

EntireCohort_pathvar$ZMPSTE24[is.na(EntireCohort_pathvar$ZMPSTE24)] <- 0
EntireCohort_pathvar$ZMPSTE24 <- as.numeric(EntireCohort_pathvar$ZMPSTE24)
unique(EntireCohort_pathvar$ZMPSTE24)
# 0 1

table(EntireCohort_pathvar$ZMPSTE24)
#  0      1 
# 469152    647




### 'FLNA', 
# no varians present 




### 'FLNB', 
EntireCohort_pathvar$FLNB <- NA
EntireCohort_pathvar$FLNB[EntireCohort_pathvar$GENE1 %in% c("FLNB")] <- 1
EntireCohort_pathvar$FLNB[EntireCohort_pathvar$GENE2 %in% c("FLNB")] <- 1
EntireCohort_pathvar$FLNB[EntireCohort_pathvar$GENE3 %in% c("FLNB")] <- 1
EntireCohort_pathvar$FLNB[EntireCohort_pathvar$GENE4 %in% c("FLNB")] <- 1
EntireCohort_pathvar$FLNB[EntireCohort_pathvar$GENE5 %in% c("FLNB")] <- 1
EntireCohort_pathvar$FLNB[EntireCohort_pathvar$GENE6 %in% c("FLNB")] <- 1
EntireCohort_pathvar$FLNB[EntireCohort_pathvar$GENE7 %in% c("FLNB")] <- 1
EntireCohort_pathvar$FLNB[EntireCohort_pathvar$GENE8 %in% c("FLNB")] <- 1
EntireCohort_pathvar$FLNB[EntireCohort_pathvar$GENE9 %in% c("FLNB")] <- 1

EntireCohort_pathvar$FLNB[is.na(EntireCohort_pathvar$FLNB)] <- 0
EntireCohort_pathvar$FLNB <- as.numeric(EntireCohort_pathvar$FLNB)
unique(EntireCohort_pathvar$FLNB)
# 0 1

table(EntireCohort_pathvar$FLNB)
#   0      1 
# 455054  14745




### 'TUBA1A', 
# no variants present





### 'TUBB2B', 
EntireCohort_pathvar$TUBB2B <- NA
EntireCohort_pathvar$TUBB2B[EntireCohort_pathvar$GENE1 %in% c("TUBB2B")] <- 1
EntireCohort_pathvar$TUBB2B[EntireCohort_pathvar$GENE2 %in% c("TUBB2B")] <- 1
EntireCohort_pathvar$TUBB2B[EntireCohort_pathvar$GENE3 %in% c("TUBB2B")] <- 1
EntireCohort_pathvar$TUBB2B[EntireCohort_pathvar$GENE4 %in% c("TUBB2B")] <- 1
EntireCohort_pathvar$TUBB2B[EntireCohort_pathvar$GENE5 %in% c("TUBB2B")] <- 1
EntireCohort_pathvar$TUBB2B[EntireCohort_pathvar$GENE6 %in% c("TUBB2B")] <- 1
EntireCohort_pathvar$TUBB2B[EntireCohort_pathvar$GENE7 %in% c("TUBB2B")] <- 1
EntireCohort_pathvar$TUBB2B[EntireCohort_pathvar$GENE8 %in% c("TUBB2B")] <- 1
EntireCohort_pathvar$TUBB2B[EntireCohort_pathvar$GENE9 %in% c("TUBB2B")] <- 1

EntireCohort_pathvar$TUBB2B[is.na(EntireCohort_pathvar$TUBB2B)] <- 0
EntireCohort_pathvar$TUBB2B <- as.numeric(EntireCohort_pathvar$TUBB2B)
unique(EntireCohort_pathvar$TUBB2B)
# 0 1

table(EntireCohort_pathvar$TUBB2B)
#   0      1 
# 469745     54




### 'TUBGCP4', 
EntireCohort_pathvar$TUBGCP4 <- NA
EntireCohort_pathvar$TUBGCP4[EntireCohort_pathvar$GENE1 %in% c("TUBGCP4")] <- 1
EntireCohort_pathvar$TUBGCP4[EntireCohort_pathvar$GENE2 %in% c("TUBGCP4")] <- 1
EntireCohort_pathvar$TUBGCP4[EntireCohort_pathvar$GENE3 %in% c("TUBGCP4")] <- 1
EntireCohort_pathvar$TUBGCP4[EntireCohort_pathvar$GENE4 %in% c("TUBGCP4")] <- 1
EntireCohort_pathvar$TUBGCP4[EntireCohort_pathvar$GENE5 %in% c("TUBGCP4")] <- 1
EntireCohort_pathvar$TUBGCP4[EntireCohort_pathvar$GENE6 %in% c("TUBGCP4")] <- 1
EntireCohort_pathvar$TUBGCP4[EntireCohort_pathvar$GENE7 %in% c("TUBGCP4")] <- 1
EntireCohort_pathvar$TUBGCP4[EntireCohort_pathvar$GENE8 %in% c("TUBGCP4")] <- 1
EntireCohort_pathvar$TUBGCP4[EntireCohort_pathvar$GENE9 %in% c("TUBGCP4")] <- 1

EntireCohort_pathvar$TUBGCP4[is.na(EntireCohort_pathvar$TUBGCP4)] <- 0
EntireCohort_pathvar$TUBGCP4 <- as.numeric(EntireCohort_pathvar$TUBGCP4)
unique(EntireCohort_pathvar$TUBGCP4)
# 0 1

table(EntireCohort_pathvar$TUBGCP4)
#  0      1 
# 469089    710




### 'ACTB',  
EntireCohort_pathvar$ACTB <- NA
EntireCohort_pathvar$ACTB[EntireCohort_pathvar$GENE1 %in% c("ACTB")] <- 1
EntireCohort_pathvar$ACTB[EntireCohort_pathvar$GENE2 %in% c("ACTB")] <- 1
EntireCohort_pathvar$ACTB[EntireCohort_pathvar$GENE3 %in% c("ACTB")] <- 1
EntireCohort_pathvar$ACTB[EntireCohort_pathvar$GENE4 %in% c("ACTB")] <- 1
EntireCohort_pathvar$ACTB[EntireCohort_pathvar$GENE5 %in% c("ACTB")] <- 1
EntireCohort_pathvar$ACTB[EntireCohort_pathvar$GENE6 %in% c("ACTB")] <- 1
EntireCohort_pathvar$ACTB[EntireCohort_pathvar$GENE7 %in% c("ACTB")] <- 1
EntireCohort_pathvar$ACTB[EntireCohort_pathvar$GENE8 %in% c("ACTB")] <- 1
EntireCohort_pathvar$ACTB[EntireCohort_pathvar$GENE9 %in% c("ACTB")] <- 1

EntireCohort_pathvar$ACTB[is.na(EntireCohort_pathvar$ACTB)] <- 0
EntireCohort_pathvar$ACTB <- as.numeric(EntireCohort_pathvar$ACTB)
unique(EntireCohort_pathvar$ACTB)
# 0 1

table(EntireCohort_pathvar$ACTB)
# 0      1 
# 469787     12




### 'ACTG1', 
# no variants present




### 'MACF1', 
# no variants present




### 'VPS4A'  
# no variants present




# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)

















# CentCytgene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_CentCytgene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_CentCytgene_, "pathvar_CentCytgene_.csv", row.names=FALSE)
# system('dx upload pathvar_CentCytgene_.csv --path Emily-folder/results_2/pathvar_CentCytgene_.csv')


# CentCytgene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ CentCytgene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_CentCytgene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_incidence_1, "pathvar_CentCytgene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_incidence_1.csv --path Emily-folder/results_2/pathvar_CentCytgene_incidence_1.csv')


# GOIs
LR_2 <- glm(formula = cancer_ROOC_dummy ~ POC1A+CDK6+ARHGEF2+ASPM+BUB1+CDK5RAP2+CENPE+CENPF+CENPJ+CENPT+CEP135+CEP152+CEP63+CIT+CKAP2L+KATNB1+KIF11+KIFBP+KNL1+NDE1+NIN+PCNT+PLK4+RTTN+SASS6+STIL+TRAPPC14+TUBGCP6+WDR62+NEK1+CFAP410+ZMPSTE24+FLNB+TUBB2B+TUBGCP4+ACTB, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_CentCytgene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_incidence_2, "pathvar_CentCytgene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_incidence_2.csv --path Emily-folder/results_2/pathvar_CentCytgene_incidence_2.csv')







### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "CentCytgene", 'POC1A', 'CDK6', 'ARHGEF2', 'ASPM', 
                                     'BUB1', 'CDK5RAP2', 'CENPE', 'CENPF', 'CENPJ', 'CENPT', 
                                     'CEP135', 'CEP152', 'CEP63', 'CIT', 'CKAP2L', 'KATNB1', 
                                     'KIF11', 'KIFBP', 'KNL1', 'NDE1', 'NIN', 
                                     'PCNT', 'PLK4', 'RTTN', 'SASS6', 'STIL', 'TRAPPC14', 
                                     'TUBGCP6', 'WDR62', 'NEK1', 'CFAP410',  
                                     'ZMPSTE24', 'FLNB', 'TUBB2B', 
                                     'TUBGCP4', 'ACTB')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
#465881     43

with(data_setA, table(ACTB))
# TUBB2B = 52
# ACTB = 12





# CentCytgene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ CentCytgene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_CentCytgene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_incidence_3, "pathvar_CentCytgene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_incidence_3.csv --path Emily-folder/results_2/pathvar_CentCytgene_incidence_3.csv')


# age + ethnicity + smoking + sex + CentCytgene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_CentCytgene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_incidence_4, "pathvar_CentCytgene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_incidence_4.csv --path Emily-folder/results_2/pathvar_CentCytgene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ POC1A+CDK6+ARHGEF2+ASPM+BUB1+CDK5RAP2+CENPE+CENPF+CENPJ+CENPT+CEP135+CEP152+CEP63+CIT+CKAP2L+KATNB1+KIF11+KIFBP+KNL1+NDE1+NIN+PCNT+PLK4+RTTN+SASS6+STIL+TRAPPC14+TUBGCP6+WDR62+NEK1+CFAP410+ZMPSTE24+FLNB+TUBB2B+TUBGCP4+ACTB, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_CentCytgene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_incidence_5, "pathvar_CentCytgene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_incidence_5.csv --path Emily-folder/results_2/pathvar_CentCytgene_incidence_5.csv')


# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POC1A+CDK6+ARHGEF2+ASPM+BUB1+CDK5RAP2+CENPE+CENPF+CENPJ+CENPT+CEP135+CEP152+CEP63+CIT+CKAP2L+KATNB1+KIF11+KIFBP+KNL1+NDE1+NIN+PCNT+PLK4+RTTN+SASS6+STIL+TRAPPC14+TUBGCP6+WDR62+NEK1+CFAP410+ZMPSTE24+FLNB+TUBB2B+TUBGCP4+ACTB, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_CentCytgene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_incidence_6, "pathvar_CentCytgene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_incidence_6.csv --path Emily-folder/results_2/pathvar_CentCytgene_incidence_6.csv')




### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "CentCytgene", 'POC1A', 'CDK6', 'ARHGEF2', 'ASPM', 
                                     'BUB1', 'CDK5RAP2', 'CENPE', 'CENPF', 'CENPJ', 'CENPT', 
                                     'CEP135', 'CEP152', 'CEP63', 'CIT', 'CKAP2L', 'KATNB1', 
                                     'KIF11', 'KIFBP', 'KNL1', 'NDE1', 'NIN', 
                                     'PCNT', 'PLK4', 'RTTN', 'SASS6', 'STIL', 'TRAPPC14', 
                                     'TUBGCP6', 'WDR62', 'NEK1', 'CFAP410', 
                                     'ZMPSTE24', 'FLNB', 'TUBB2B', 
                                     'TUBGCP4', 'ACTB')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     48

with(data_setB, table(ACTB))
# CDK6 = 99
# TUBB2B = 41
# ACTB = 10





# CentCytgene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ CentCytgene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_CentCytgene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_incidence_7, "pathvar_CentCytgene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_incidence_7.csv --path Emily-folder/results_2/pathvar_CentCytgene_incidence_7.csv')


# all variables + CentCytgene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+CentCytgene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_CentCytgene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_incidence_8, "pathvar_CentCytgene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_incidence_8.csv --path Emily-folder/results_2/pathvar_CentCytgene_incidence_8.csv')


# GOIs alone
LR_9 <- glm(formula = cancer_ROOC_dummy ~ POC1A+CDK6+ARHGEF2+ASPM+BUB1+CDK5RAP2+CENPE+CENPF+CENPJ+CENPT+CEP135+CEP152+CEP63+CIT+CKAP2L+KATNB1+KIF11+KIFBP+KNL1+NDE1+NIN+PCNT+PLK4+RTTN+SASS6+STIL+TRAPPC14+TUBGCP6+WDR62+NEK1+CFAP410+ZMPSTE24+FLNB+TUBB2B+TUBGCP4+ACTB, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_CentCytgene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_incidence_9, "pathvar_CentCytgene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_incidence_9.csv --path Emily-folder/results_2/pathvar_CentCytgene_incidence_9.csv')


# all variables + GOIs
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+POC1A+CDK6+ARHGEF2+ASPM+BUB1+CDK5RAP2+CENPE+CENPF+CENPJ+CENPT+CEP135+CEP152+CEP63+CIT+CKAP2L+KATNB1+KIF11+KIFBP+KNL1+NDE1+NIN+PCNT+PLK4+RTTN+SASS6+STIL+TRAPPC14+TUBGCP6+WDR62+NEK1+CFAP410+ZMPSTE24+FLNB+TUBB2B+TUBGCP4+ACTB, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_CentCytgene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_incidence_10, "pathvar_CentCytgene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_incidence_10.csv --path Emily-folder/results_2/pathvar_CentCytgene_incidence_10.csv')










### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(CentCytgene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(CentCytgene, cancer_ROOC_dummy)))
# p-value = 1.129e-06
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9339476 0.9714315
# sample estimates:
#   odds ratio 
# 0.9525361

with(EntireCohort_pathvar, table(CentCytgene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(CentCytgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5562106 0.5926229
# sample estimates:
#   odds ratio 
# 0.5740907

with(EntireCohort_pathvar, table(CentCytgene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(CentCytgene, Ever_smoked_dummy)))
# p-value = 1.272e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.947059 0.979554
# sample estimates:
#   odds ratio 
# 0.9631353

with(EntireCohort_pathvar, table(CentCytgene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(CentCytgene, Sex_Male)))
# p-value = 0.366

with(EntireCohort_pathvar, table(CentCytgene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(CentCytgene, University_education)))
# p-value = 0.1616

with(EntireCohort_pathvar, table(CentCytgene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(CentCytgene, Annual_income))	)
# p-value = 0.3089

with(EntireCohort_pathvar, table(CentCytgene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(CentCytgene, Smoker_previous)))
# p-value = 1.961e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9460265 0.9796993
# sample estimates:
#   odds ratio 
# 0.962705

with(EntireCohort_pathvar, table(CentCytgene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(CentCytgene, Smoker_current))	)
# p-value = 0.8962



### data_setA
with(data_setA, table(CentCytgene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(CentCytgene, cancer_ROOC_dummy)))
# p-value = 9.598e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9333783 0.9709947
# sample estimates:
#   odds ratio 
# 0.9520222

with(data_setA, table(CentCytgene, Ethnicity_White))			
fisher.test(with(data_setA, table(CentCytgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5561442 0.5926992
# sample estimates:
#   odds ratio 
# 0.5740933

with(data_setA, table(CentCytgene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(CentCytgene, Ever_smoked_dummy)))
# p-value = 1.961e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9477320 0.9802821
# sample estimates:
#   odds ratio 
# 0.9638545

with(data_setA, table(CentCytgene, Sex_Male))			
fisher.test(with(data_setA, table(CentCytgene, Sex_Male)))
# p-value = 0.4945



### data_setB
with(data_setB, table(CentCytgene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(CentCytgene, cancer_ROOC_dummy)))
# p-value = 8.507e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9361682 0.9782911
# sample estimates:
#   odds ratio 
# 0.9570424

with(data_setB, table(CentCytgene, Ethnicity_White))			
fisher.test(with(data_setB, table(CentCytgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5666326 0.6119772
# sample estimates:
#   odds ratio 
# 0.5888041

with(data_setB, table(CentCytgene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(CentCytgene, Ever_smoked_dummy)))
# p-value = 0.002456
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9533526 0.9898671
# sample estimates:
#   odds ratio 
# 0.9714117

with(data_setB, table(CentCytgene, Sex_Male))			
fisher.test(with(data_setB, table(CentCytgene, Sex_Male)))
# p-value = 0.6525

with(data_setB, table(CentCytgene, University_education))			
fisher.test(with(data_setB, table(CentCytgene, University_education)))
# p-value = 0.11

with(data_setB, table(CentCytgene, Annual_income))			
fisher.test(with(data_setB, table(CentCytgene, Annual_income)))
# p-value = 0.4933






### can plot the continuous variables: 

# make CentCytgene into factor
df_for_plot <- EntireCohort_pathvar
df_for_plot$CentCytgene <- as.factor(df_for_plot$CentCytgene)

# try a violin plot

# Age_at_recruitment  
df_for_plot %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ CentCytgene, data = df_for_plot)
# W = 1.343e+10, p-value = 6.257e-08


# BMI  
df_for_plot %>% 
  ggplot(mapping = aes(x = CentCytgene, y = BMI)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = CentCytgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ CentCytgene, data = df_for_plot)
# W = 1.3058e+10, p-value = 0.009021

# Standing_height 
df_for_plot %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ CentCytgene, data = df_for_plot)
# W = 1.3323e+10, p-value = 1.648e-06

# make this a bar plot
# Weekly_exercise
df_for_plot %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(CentCytgene), scales = "free_y") +
  theme_light()
# 





# save plots --------------------------------------------------------------


### data_setA

# make CentCytgene into factor
data_setA$CentCytgene <- as.factor(data_setA$CentCytgene)

# Age_at_recruitment  
jpeg(file="CentCytgene_Age_v_A.jpeg", width = 300, height = 300)
data_setA %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
dev.off()
system('dx upload CentCytgene_Age_v_A.jpeg --path Emily-folder/pathvar/plots/CentCytgene_Age_v_A.jpeg')

jpeg(file="CentCytgene_Age_b_A.jpeg", width = 300, height = 300)
data_setA %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
dev.off()
system('dx upload CentCytgene_Age_b_A.jpeg --path Emily-folder/pathvar/plots/CentCytgene_Age_b_A.jpeg')

wilcox.test(Age_at_recruitment ~ CentCytgene, data = data_setA)
# W = 1.3199e+10, p-value = 7.8e-08





### data_setB

# make CentCytgene into factor
data_setB$CentCytgene <- as.factor(data_setB$CentCytgene)

# Age_at_recruitment  

jpeg(file="CentCytgene_Age_v_B.jpeg", width = 300, height = 300)
data_setB %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
dev.off()
system('dx upload CentCytgene_Age_v_B.jpeg --path Emily-folder/pathvar/plots/CentCytgene_Age_v_B.jpeg')


jpeg(file="CentCytgene_Age_b_B.jpeg", width = 300, height = 300)
data_setB %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
dev.off()
system('dx upload CentCytgene_Age_b_B.jpeg --path Emily-folder/pathvar/plots/CentCytgene_Age_b_B.jpeg')

wilcox.test(Age_at_recruitment ~ CentCytgene, data = data_setB)
# W = 8756296456, p-value = 0.0008004



# BMI 

jpeg(file="CentCytgene_BMI_v_B.jpeg", width = 300, height = 300)
data_setB %>% 
  ggplot(mapping = aes(x = CentCytgene, y = BMI)) +
  geom_violin() +
  theme_light()
dev.off()
system('dx upload CentCytgene_BMI_v_B.jpeg --path Emily-folder/pathvar/plots/CentCytgene_BMI_v_B.jpeg')

jpeg(file="CentCytgene_BMI_b_B.jpeg", width = 300, height = 300)
data_setB %>% 
  ggplot(mapping = aes(x = CentCytgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
dev.off()
system('dx upload CentCytgene_BMI_b_B.jpeg --path Emily-folder/pathvar/plots/CentCytgene_BMI_b_B.jpeg')

wilcox.test(BMI ~ CentCytgene, data = data_setB)
# W = 8625186722, p-value = 0.0255




# Standing_height 

jpeg(file="CentCytgene_Height_v_B.jpeg", width = 300, height = 300)
data_setB %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
dev.off()
system('dx upload CentCytgene_Height_v_B.jpeg --path Emily-folder/pathvar/plots/CentCytgene_Height_v_B.jpeg')

jpeg(file="CentCytgene_Height_b_B.jpeg", width = 300, height = 300)
data_setB %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
dev.off()
system('dx upload CentCytgene_Height_b_B.jpeg --path Emily-folder/pathvar/plots/CentCytgene_Height_b_B.jpeg')

wilcox.test(Standing_height ~ CentCytgene, data = data_setB)
# W = 8778870334, p-value = 1.605e-05





# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(CentCytgene), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# CentCytgene Type ------------------------------------------------------------

# make sure the execute the code in the 'CentCytgenes' section


## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_CentCytgene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_type_1, "pathvar_CentCytgene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_type_1.csv --path Emily-folder/results_2/pathvar_CentCytgene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_CentCytgene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_type_2, "pathvar_CentCytgene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_type_2.csv --path Emily-folder/results_2/pathvar_CentCytgene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_CentCytgene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_type_3, "pathvar_CentCytgene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_type_3.csv --path Emily-folder/results_2/pathvar_CentCytgene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_CentCytgene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_type_4, "pathvar_CentCytgene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_type_4.csv --path Emily-folder/results_2/pathvar_CentCytgene_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_CentCytgene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_type_5, "pathvar_CentCytgene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_type_5.csv --path Emily-folder/results_2/pathvar_CentCytgene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_CentCytgene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_type_6, "pathvar_CentCytgene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_type_6.csv --path Emily-folder/results_2/pathvar_CentCytgene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_CentCytgene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_type_7, "pathvar_CentCytgene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_type_7.csv --path Emily-folder/results_2/pathvar_CentCytgene_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_CentCytgene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_type_8, "pathvar_CentCytgene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_type_8.csv --path Emily-folder/results_2/pathvar_CentCytgene_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_CentCytgene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_type_9, "pathvar_CentCytgene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_type_9.csv --path Emily-folder/results_2/pathvar_CentCytgene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_CentCytgene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_type_10, "pathvar_CentCytgene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_type_10.csv --path Emily-folder/results_2/pathvar_CentCytgene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_CentCytgene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_type_11, "pathvar_CentCytgene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_type_11.csv --path Emily-folder/results_2/pathvar_CentCytgene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_CentCytgene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_type_12, "pathvar_CentCytgene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_type_12.csv --path Emily-folder/results_2/pathvar_CentCytgene_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_CentCytgene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_type_13, "pathvar_CentCytgene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_type_13.csv --path Emily-folder/results_2/pathvar_CentCytgene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_CentCytgene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_type_14, "pathvar_CentCytgene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_type_14.csv --path Emily-folder/results_2/pathvar_CentCytgene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_CentCytgene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CentCytgene_type_15, "pathvar_CentCytgene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_type_15.csv --path Emily-folder/results_2/pathvar_CentCytgene_type_15.csv')







# CKAP2L Type -------------------------------------------------------------
# repeat the type analysis for just CKAP2L




LR_oral_CKAP2L <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CKAP2L, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_CKAP2L)
# exp(cbind(OR = coef(LR_oral_CKAP2L), confint(LR_oral_CKAP2L)))
pathvar_CKAP2L_type_1 <- tidy(LR_oral_CKAP2L, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CKAP2L_type_1, "pathvar_CKAP2L_type_1.csv", row.names=FALSE)
system('dx upload pathvar_CKAP2L_type_1.csv --path Emily-folder/results_2/pathvar_CKAP2L_type_1.csv')



LR_digestive_CKAP2L <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CKAP2L, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_CKAP2L)
# exp(cbind(OR = coef(LR_digestive_CKAP2L), confint(LR_digestive_CKAP2L)))
pathvar_CKAP2L_type_2 <- tidy(LR_digestive_CKAP2L, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CKAP2L_type_2, "pathvar_CKAP2L_type_2.csv", row.names=FALSE)
system('dx upload pathvar_CKAP2L_type_2.csv --path Emily-folder/results_2/pathvar_CKAP2L_type_2.csv')



LR_respiratory_CKAP2L <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CKAP2L, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_CKAP2L)
# exp(cbind(OR = coef(LR_respiratory_CKAP2L), confint(LR_respiratory_CKAP2L)))
pathvar_CKAP2L_type_3 <- tidy(LR_respiratory_CKAP2L, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CKAP2L_type_3, "pathvar_CKAP2L_type_3.csv", row.names=FALSE)
system('dx upload pathvar_CKAP2L_type_3.csv --path Emily-folder/results_2/pathvar_CKAP2L_type_3.csv')



LR_bone_CKAP2L <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CKAP2L, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_CKAP2L)
# exp(cbind(OR = coef(LR_bone_CKAP2L), confint(LR_bone_CKAP2L)))
pathvar_CKAP2L_type_4 <- tidy(LR_bone_CKAP2L, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CKAP2L_type_4, "pathvar_CKAP2L_type_4.csv", row.names=FALSE)
system('dx upload pathvar_CKAP2L_type_4.csv --path Emily-folder/results_2/pathvar_CKAP2L_type_4.csv')



LR_skin_CKAP2L <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CKAP2L, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_CKAP2L)
# exp(cbind(OR = coef(LR_skin_CKAP2L), confint(LR_skin_CKAP2L)))
pathvar_CKAP2L_type_5 <- tidy(LR_skin_CKAP2L, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CKAP2L_type_5, "pathvar_CKAP2L_type_5.csv", row.names=FALSE)
system('dx upload pathvar_CKAP2L_type_5.csv --path Emily-folder/results_2/pathvar_CKAP2L_type_5.csv')



LR_mesothelium_CKAP2L <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CKAP2L, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_CKAP2L)
# exp(cbind(OR = coef(LR_mesothelium_CKAP2L), confint(LR_mesothelium_CKAP2L)))
pathvar_CKAP2L_type_6 <- tidy(LR_mesothelium_CKAP2L, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CKAP2L_type_6, "pathvar_CKAP2L_type_6.csv", row.names=FALSE)
system('dx upload pathvar_CKAP2L_type_6.csv --path Emily-folder/results_2/pathvar_CKAP2L_type_6.csv')



LR_breast_CKAP2L <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CKAP2L, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_CKAP2L)
# exp(cbind(OR = coef(LR_breast_CKAP2L), confint(LR_breast_CKAP2L)))
pathvar_CKAP2L_type_7 <- tidy(LR_breast_CKAP2L, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CKAP2L_type_7, "pathvar_CKAP2L_type_7.csv", row.names=FALSE)
system('dx upload pathvar_CKAP2L_type_7.csv --path Emily-folder/results_2/pathvar_CKAP2L_type_7.csv')



LR_femalegenital_CKAP2L <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CKAP2L, 
                               data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_CKAP2L)
# exp(cbind(OR = coef(LR_femalegenital_CKAP2L), confint(LR_femalegenital_CKAP2L)))
pathvar_CKAP2L_type_8 <- tidy(LR_femalegenital_CKAP2L, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CKAP2L_type_8, "pathvar_CKAP2L_type_8.csv", row.names=FALSE)
system('dx upload pathvar_CKAP2L_type_8.csv --path Emily-folder/results_2/pathvar_CKAP2L_type_8.csv')



LR_malegenital_CKAP2L <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CKAP2L, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_CKAP2L)
# exp(cbind(OR = coef(LR_malegenital_CKAP2L), confint(LR_malegenital_CKAP2L)))
pathvar_CKAP2L_type_9 <- tidy(LR_malegenital_CKAP2L, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CKAP2L_type_9, "pathvar_CKAP2L_type_9.csv", row.names=FALSE)
system('dx upload pathvar_CKAP2L_type_9.csv --path Emily-folder/results_2/pathvar_CKAP2L_type_9.csv')



LR_urinary_CKAP2L <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CKAP2L, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_CKAP2L)
# exp(cbind(OR = coef(LR_urinary_CKAP2L), confint(LR_urinary_CKAP2L)))
pathvar_CKAP2L_type_10 <- tidy(LR_urinary_CKAP2L, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CKAP2L_type_10, "pathvar_CKAP2L_type_10.csv", row.names=FALSE)
system('dx upload pathvar_CKAP2L_type_10.csv --path Emily-folder/results_2/pathvar_CKAP2L_type_10.csv')



LR_cns_CKAP2L <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CKAP2L, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_CKAP2L)
# exp(cbind(OR = coef(LR_cns_CKAP2L), confint(LR_cns_CKAP2L)))
pathvar_CKAP2L_type_11 <- tidy(LR_cns_CKAP2L, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CKAP2L_type_11, "pathvar_CKAP2L_type_11.csv", row.names=FALSE)
system('dx upload pathvar_CKAP2L_type_11.csv --path Emily-folder/results_2/pathvar_CKAP2L_type_11.csv')



LR_endocrine_CKAP2L <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CKAP2L, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_CKAP2L)
# exp(cbind(OR = coef(LR_endocrine_CKAP2L), confint(LR_endocrine_CKAP2L)))
pathvar_CKAP2L_type_12 <- tidy(LR_endocrine_CKAP2L, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CKAP2L_type_12, "pathvar_CKAP2L_type_12.csv", row.names=FALSE)
system('dx upload pathvar_CKAP2L_type_12.csv --path Emily-folder/results_2/pathvar_CKAP2L_type_12.csv')



LR_lymphatic_CKAP2L <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CKAP2L, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_CKAP2L)
# exp(cbind(OR = coef(LR_lymphatic_CKAP2L), confint(LR_lymphatic_CKAP2L)))
pathvar_CKAP2L_type_13 <- tidy(LR_lymphatic_CKAP2L, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CKAP2L_type_13, "pathvar_CKAP2L_type_13.csv", row.names=FALSE)
system('dx upload pathvar_CKAP2L_type_13.csv --path Emily-folder/results_2/pathvar_CKAP2L_type_13.csv')



LR_secondary_CKAP2L <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CKAP2L, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_CKAP2L)
# exp(cbind(OR = coef(LR_secondary_CKAP2L), confint(LR_secondary_CKAP2L)))
pathvar_CKAP2L_type_14 <- tidy(LR_secondary_CKAP2L, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CKAP2L_type_14, "pathvar_CKAP2L_type_14.csv", row.names=FALSE)
system('dx upload pathvar_CKAP2L_type_14.csv --path Emily-folder/results_2/pathvar_CKAP2L_type_14.csv')



LR_other_CKAP2L <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CKAP2L, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_CKAP2L)
# exp(cbind(OR = coef(LR_other_CKAP2L), confint(LR_other_CKAP2L)))
pathvar_CKAP2L_type_15 <- tidy(LR_other_CKAP2L, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CKAP2L_type_15, "pathvar_CKAP2L_type_15.csv", row.names=FALSE)
system('dx upload pathvar_CKAP2L_type_15.csv --path Emily-folder/results_2/pathvar_CKAP2L_type_15.csv')











# CentCytgene Age -------------------------------------------------------------

# make sure the execute the code in the 'CentCytgenes' section

### linear regression
# "CentCytgene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ CentCytgene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_CentCytgene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_CentCytgene_age_1, "pathvar_CentCytgene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_age_1.csv --path Emily-folder/results_2/pathvar_CentCytgene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ POC1A+CDK6+ARHGEF2+ASPM+BUB1+CDK5RAP2+CENPE+CENPF+CENPJ+CENPT+CEP135+CEP152+CEP63+CIT+CKAP2L+KATNB1+KIF11+KIFBP+KNL1+NDE1+NIN+PCNT+PLK4+RTTN+SASS6+STIL+TRAPPC14+TUBGCP6+WDR62+NEK1+CFAP410+ZMPSTE24+FLNB+TUBB2B+TUBGCP4+ACTB, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_CentCytgene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_CentCytgene_age_2, "pathvar_CentCytgene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_age_2.csv --path Emily-folder/results_2/pathvar_CentCytgene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "CentCytgene", 'POC1A', 'CDK6', 'ARHGEF2', 'ASPM', 'BUB1', 
                                     'CDK5RAP2', 'CENPE', 'CENPF', 'CENPJ', 'CENPT', 'CEP135', 
                                     'CEP152', 'CEP63', 'CIT', 'CKAP2L', 'KATNB1', 'KIF11', 
                                     'KIFBP', 'KNL1', 'NDE1', 'NIN', 'PCNT', 'PLK4', 
                                     'RTTN', 'SASS6', 'STIL', 'TRAPPC14', 'TUBGCP6', 'WDR62', 
                                     'NEK1', 'CFAP410', 'ZMPSTE24', 
                                     'FLNB', 'TUBB2B', 'TUBGCP4', 'ACTB')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     43


with(data_setD, table(ACTB))
# CDK6 = 31
# TUBB2B = 15
# ACTB = 3




data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "CentCytgene", 'POC1A', 'CDK6', 'ARHGEF2', 
                                     'ASPM', 'BUB1', 'CDK5RAP2', 'CENPE', 'CENPF', 'CENPJ', 'CENPT', 
                                     'CEP135', 'CEP152', 'CEP63', 'CIT', 'CKAP2L', 'KATNB1', 'KIF11', 
                                     'KIFBP', 'KNL1', 'NDE1', 'NIN', 'PCNT', 'PLK4', 'RTTN', 
                                     'SASS6', 'STIL', 'TRAPPC14', 'TUBGCP6', 'WDR62', 'NEK1', 
                                     'CFAP410', 'ZMPSTE24', 
                                     'FLNB', 'TUBB2B', 'TUBGCP4', 'ACTB')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    48


with(data_setE, table(ACTB))
# CDK6 = 25
# CKAP2L = 88
# KIF11 = 90
# KNL1 = 85
# NDE1 = 84
# TRAPPC14 = 93
# TUBB2B = 11
# ACTB = 2



### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ CentCytgene, data = data_setD)
# summary(LM_15)
pathvar_CentCytgene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_CentCytgene_age_3, "pathvar_CentCytgene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_age_3.csv --path Emily-folder/results_2/pathvar_CentCytgene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CentCytgene, data = data_setD)
# summary(LM_16)
pathvar_CentCytgene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_CentCytgene_age_4, "pathvar_CentCytgene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_age_4.csv --path Emily-folder/results_2/pathvar_CentCytgene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ POC1A+CDK6+ARHGEF2+ASPM+BUB1+CDK5RAP2+CENPE+CENPF+CENPJ+CENPT+CEP135+CEP152+CEP63+CIT+CKAP2L+KATNB1+KIF11+KIFBP+KNL1+NDE1+NIN+PCNT+PLK4+RTTN+SASS6+STIL+TRAPPC14+TUBGCP6+WDR62+NEK1+CFAP410+ZMPSTE24+FLNB+TUBB2B+TUBGCP4+ACTB, data = data_setD)
# summary(LM_17)
pathvar_CentCytgene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_CentCytgene_age_5, "pathvar_CentCytgene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_age_5.csv --path Emily-folder/results_2/pathvar_CentCytgene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POC1A+CDK6+ARHGEF2+ASPM+BUB1+CDK5RAP2+CENPE+CENPF+CENPJ+CENPT+CEP135+CEP152+CEP63+CIT+CKAP2L+KATNB1+KIF11+KIFBP+KNL1+NDE1+NIN+PCNT+PLK4+RTTN+SASS6+STIL+TRAPPC14+TUBGCP6+WDR62+NEK1+CFAP410+ZMPSTE24+FLNB+TUBB2B+TUBGCP4+ACTB, data = data_setD)
# summary(LM_18)
pathvar_CentCytgene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_CentCytgene_age_6, "pathvar_CentCytgene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_age_6.csv --path Emily-folder/results_2/pathvar_CentCytgene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ CentCytgene, data = data_setE)
# summary(LM_15)
pathvar_CentCytgene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_CentCytgene_age_7, "pathvar_CentCytgene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_age_7.csv --path Emily-folder/results_2/pathvar_CentCytgene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+CentCytgene, data = data_setE)
# summary(LM_16)
pathvar_CentCytgene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_CentCytgene_age_8, "pathvar_CentCytgene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_age_8.csv --path Emily-folder/results_2/pathvar_CentCytgene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ POC1A+CDK6+ARHGEF2+ASPM+BUB1+CDK5RAP2+CENPE+CENPF+CENPJ+CENPT+CEP135+CEP152+CEP63+CIT+CKAP2L+KATNB1+KIF11+KIFBP+KNL1+NDE1+NIN+PCNT+PLK4+RTTN+SASS6+STIL+TRAPPC14+TUBGCP6+WDR62+NEK1+CFAP410+ZMPSTE24+FLNB+TUBB2B+TUBGCP4+ACTB, data = data_setE)
# summary(LM_17)
pathvar_CentCytgene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_CentCytgene_age_9, "pathvar_CentCytgene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_age_9.csv --path Emily-folder/results_2/pathvar_CentCytgene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+POC1A+CDK6+ARHGEF2+ASPM+BUB1+CDK5RAP2+CENPE+CENPF+CENPJ+CENPT+CEP135+CEP152+CEP63+CIT+CKAP2L+KATNB1+KIF11+KIFBP+KNL1+NDE1+NIN+PCNT+PLK4+RTTN+SASS6+STIL+TRAPPC14+TUBGCP6+WDR62+NEK1+CFAP410+ZMPSTE24+FLNB+TUBB2B+TUBGCP4+ACTB, data = data_setE)
# summary(LM_18)
pathvar_CentCytgene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_CentCytgene_age_10, "pathvar_CentCytgene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_CentCytgene_age_10.csv --path Emily-folder/results_2/pathvar_CentCytgene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(df_cancer_diagnosis, table(CentCytgene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(CentCytgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5348349 0.6421838
# sample estimates:
#   odds ratio 
# 0.5857165 

with(df_cancer_diagnosis, table(CentCytgene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(CentCytgene, Ever_smoked_dummy)))
# p-value = 0.1706

with(df_cancer_diagnosis, table(CentCytgene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(CentCytgene, Sex_Male)))
# p-value = 0.2306

with(df_cancer_diagnosis, table(CentCytgene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(CentCytgene, University_education)))
# p-value = 0.009509
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.012200 1.091441
# sample estimates:
#   odds ratio 
# 1.051118

with(df_cancer_diagnosis, table(CentCytgene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(CentCytgene, Annual_income)))
# p-value = 0.3701


### data_setD
with(data_setD, table(CentCytgene, Ethnicity_White))			
fisher.test(with(data_setD, table(CentCytgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5348931 0.6428253
# sample estimates:
#   odds ratio 
# 0.5860343

with(data_setD, table(CentCytgene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(CentCytgene, Ever_smoked_dummy)))
# p-value = 0.1508

with(data_setD, table(CentCytgene, Sex_Male))			
fisher.test(with(data_setD, table(CentCytgene, Sex_Male)))
# p-value = 0.3212


### data_setE
with(data_setE, table(CentCytgene, Ethnicity_White))			
fisher.test(with(data_setE, table(CentCytgene, Ethnicity_White)))
# p-value = 2.617e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5482282 0.6877785
# sample estimates:
#   odds ratio 
# 0.6135075

with(data_setE, table(CentCytgene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(CentCytgene, Ever_smoked_dummy)))
# p-value = 0.3309

with(data_setE, table(CentCytgene, Sex_Male))			
fisher.test(with(data_setE, table(CentCytgene, Sex_Male)))
# p-value = 0.153

with(data_setE, table(CentCytgene, University_education))			
fisher.test(with(data_setE, table(CentCytgene, University_education)))
# p-value = 0.02399
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.006083 1.091940
# sample estimates:
#   odds ratio 
# 1.048182

with(data_setE, table(CentCytgene, Annual_income))			
fisher.test(with(data_setE, table(CentCytgene, Annual_income)))
# p-value = 0.2855



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 

# make CentCytgene into factor
df_cancer_diagnosis$CentCytgene <- as.factor(df_cancer_diagnosis$CentCytgene)

# try a violin plot



# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ CentCytgene, data = df_cancer_diagnosis)
# 



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ CentCytgene, data = df_cancer_diagnosis)
# W = 723844474, p-value = 0.1492



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CentCytgene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CentCytgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ CentCytgene, data = df_cancer_diagnosis)
# W = 708287134, p-value = 0.1722



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ CentCytgene, data = df_cancer_diagnosis)
# W = 719368324, p-value = 0.2116



# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(CentCytgene), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$CentCytgene <- as.factor(data_setD$CentCytgene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ CentCytgene, data = data_setD)
# 



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ CentCytgene, data = data_setD)
# W = 711115546, p-value = 0.146




### data_setE

data_setE$CentCytgene <- as.factor(data_setE$CentCytgene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ CentCytgene, data = data_setE)
# 



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ CentCytgene, data = data_setE)
# W = 459043695, p-value = 0.2625



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = CentCytgene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = CentCytgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ CentCytgene, data = data_setE)
# W = 452133781, p-value = 0.1223



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = CentCytgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ CentCytgene, data = data_setE)
# W = 458183100, p-value = 0.4308


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(CentCytgene), scales = "free_y") +
  theme_light()
# 














# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_CentCytgenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_CentCytgenes.R')
