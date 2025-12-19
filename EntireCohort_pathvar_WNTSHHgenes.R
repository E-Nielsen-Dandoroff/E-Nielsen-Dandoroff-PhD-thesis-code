# 20/06/2024

# EntireCohort_pathvar_WNTSHHgenes.R

# this script is for the WNT and SHH associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)


# WNTSHHgene
# IHH
# GPC6
# WDR19
# CHMP1A
# ROR2
# DVL3
# EVC
# EVC2
# ICK
# IFT122
# WDR35
# GLI1
# GLI2
# GLI3
# FZD2
# HHAT
# SHH
# MYCN
# WNT1
# WNT5A
# NXN
# WDFY3







# WNTSHHgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')


# go thr all the genes to check that they're present
grep("WDFY3", EntireCohort_pathvar$GENE)
# all seem to be present except for:
# DVL3 = 0 present
# FZD2 = 0 present
# WDFY3 = 0 present


# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$WNTSHHgene <- NA
EntireCohort_pathvar$WNTSHHgene[EntireCohort_pathvar$GENE1 %in% c('IHH', 'GPC6', 'WDR19', 'CHMP1A', 'ROR2', 'EVC', 'EVC2', 'ICK', 'IFT122', 'WDR35', 'GLI1', 'GLI2', 'GLI3', 'HHAT', 'SHH', 'MYCN', 'WNT1', 'WNT5A', 'NXN')] <- 1
EntireCohort_pathvar$WNTSHHgene[EntireCohort_pathvar$GENE2 %in% c('IHH', 'GPC6', 'WDR19', 'CHMP1A', 'ROR2', 'EVC', 'EVC2', 'ICK', 'IFT122', 'WDR35', 'GLI1', 'GLI2', 'GLI3', 'HHAT', 'SHH', 'MYCN', 'WNT1', 'WNT5A', 'NXN')] <- 1
EntireCohort_pathvar$WNTSHHgene[EntireCohort_pathvar$GENE3 %in% c('IHH', 'GPC6', 'WDR19', 'CHMP1A', 'ROR2', 'EVC', 'EVC2', 'ICK', 'IFT122', 'WDR35', 'GLI1', 'GLI2', 'GLI3', 'HHAT', 'SHH', 'MYCN', 'WNT1', 'WNT5A', 'NXN')] <- 1
EntireCohort_pathvar$WNTSHHgene[EntireCohort_pathvar$GENE4 %in% c('IHH', 'GPC6', 'WDR19', 'CHMP1A', 'ROR2', 'EVC', 'EVC2', 'ICK', 'IFT122', 'WDR35', 'GLI1', 'GLI2', 'GLI3', 'HHAT', 'SHH', 'MYCN', 'WNT1', 'WNT5A', 'NXN')] <- 1
EntireCohort_pathvar$WNTSHHgene[EntireCohort_pathvar$GENE5 %in% c('IHH', 'GPC6', 'WDR19', 'CHMP1A', 'ROR2', 'EVC', 'EVC2', 'ICK', 'IFT122', 'WDR35', 'GLI1', 'GLI2', 'GLI3', 'HHAT', 'SHH', 'MYCN', 'WNT1', 'WNT5A', 'NXN')] <- 1
EntireCohort_pathvar$WNTSHHgene[EntireCohort_pathvar$GENE6 %in% c('IHH', 'GPC6', 'WDR19', 'CHMP1A', 'ROR2', 'EVC', 'EVC2', 'ICK', 'IFT122', 'WDR35', 'GLI1', 'GLI2', 'GLI3', 'HHAT', 'SHH', 'MYCN', 'WNT1', 'WNT5A', 'NXN')] <- 1
EntireCohort_pathvar$WNTSHHgene[EntireCohort_pathvar$GENE7 %in% c('IHH', 'GPC6', 'WDR19', 'CHMP1A', 'ROR2', 'EVC', 'EVC2', 'ICK', 'IFT122', 'WDR35', 'GLI1', 'GLI2', 'GLI3', 'HHAT', 'SHH', 'MYCN', 'WNT1', 'WNT5A', 'NXN')] <- 1
EntireCohort_pathvar$WNTSHHgene[EntireCohort_pathvar$GENE8 %in% c('IHH', 'GPC6', 'WDR19', 'CHMP1A', 'ROR2', 'EVC', 'EVC2', 'ICK', 'IFT122', 'WDR35', 'GLI1', 'GLI2', 'GLI3', 'HHAT', 'SHH', 'MYCN', 'WNT1', 'WNT5A', 'NXN')] <- 1
EntireCohort_pathvar$WNTSHHgene[EntireCohort_pathvar$GENE9 %in% c('IHH', 'GPC6', 'WDR19', 'CHMP1A', 'ROR2', 'EVC', 'EVC2', 'ICK', 'IFT122', 'WDR35', 'GLI1', 'GLI2', 'GLI3', 'HHAT', 'SHH', 'MYCN', 'WNT1', 'WNT5A', 'NXN')] <- 1


EntireCohort_pathvar$WNTSHHgene[is.na(EntireCohort_pathvar$WNTSHHgene)] <- 0
EntireCohort_pathvar$WNTSHHgene <- as.numeric(EntireCohort_pathvar$WNTSHHgene)
unique(EntireCohort_pathvar$WNTSHHgene)
# 0 1 

table(EntireCohort_pathvar$WNTSHHgene)
#      0      1 
# 437000  32799




# now to make a column for each of the WNTSHHgenes
### 'IHH'
EntireCohort_pathvar$IHH <- NA
EntireCohort_pathvar$IHH[EntireCohort_pathvar$GENE1 %in% c("IHH")] <- 1
EntireCohort_pathvar$IHH[EntireCohort_pathvar$GENE2 %in% c("IHH")] <- 1
EntireCohort_pathvar$IHH[EntireCohort_pathvar$GENE3 %in% c("IHH")] <- 1
EntireCohort_pathvar$IHH[EntireCohort_pathvar$GENE4 %in% c("IHH")] <- 1
EntireCohort_pathvar$IHH[EntireCohort_pathvar$GENE5 %in% c("IHH")] <- 1
EntireCohort_pathvar$IHH[EntireCohort_pathvar$GENE6 %in% c("IHH")] <- 1
EntireCohort_pathvar$IHH[EntireCohort_pathvar$GENE7 %in% c("IHH")] <- 1
EntireCohort_pathvar$IHH[EntireCohort_pathvar$GENE8 %in% c("IHH")] <- 1
EntireCohort_pathvar$IHH[EntireCohort_pathvar$GENE9 %in% c("IHH")] <- 1

EntireCohort_pathvar$IHH[is.na(EntireCohort_pathvar$IHH)] <- 0
EntireCohort_pathvar$IHH <- as.numeric(EntireCohort_pathvar$IHH)
unique(EntireCohort_pathvar$IHH)
# 0 1

table(EntireCohort_pathvar$IHH)
#      0      1 
# 469006    793 





### 'GPC6'
EntireCohort_pathvar$GPC6 <- NA
EntireCohort_pathvar$GPC6[EntireCohort_pathvar$GENE1 %in% c("GPC6")] <- 1
EntireCohort_pathvar$GPC6[EntireCohort_pathvar$GENE2 %in% c("GPC6")] <- 1
EntireCohort_pathvar$GPC6[EntireCohort_pathvar$GENE3 %in% c("GPC6")] <- 1
EntireCohort_pathvar$GPC6[EntireCohort_pathvar$GENE4 %in% c("GPC6")] <- 1
EntireCohort_pathvar$GPC6[EntireCohort_pathvar$GENE5 %in% c("GPC6")] <- 1
EntireCohort_pathvar$GPC6[EntireCohort_pathvar$GENE6 %in% c("GPC6")] <- 1
EntireCohort_pathvar$GPC6[EntireCohort_pathvar$GENE7 %in% c("GPC6")] <- 1
EntireCohort_pathvar$GPC6[EntireCohort_pathvar$GENE8 %in% c("GPC6")] <- 1
EntireCohort_pathvar$GPC6[EntireCohort_pathvar$GENE9 %in% c("GPC6")] <- 1

EntireCohort_pathvar$GPC6[is.na(EntireCohort_pathvar$GPC6)] <- 0
EntireCohort_pathvar$GPC6 <- as.numeric(EntireCohort_pathvar$GPC6)
unique(EntireCohort_pathvar$GPC6)
# 0 1

table(EntireCohort_pathvar$GPC6)
#      0      1 
# 468845    954 




### 'WDR19'
EntireCohort_pathvar$WDR19 <- NA
EntireCohort_pathvar$WDR19[EntireCohort_pathvar$GENE1 %in% c("WDR19")] <- 1
EntireCohort_pathvar$WDR19[EntireCohort_pathvar$GENE2 %in% c("WDR19")] <- 1
EntireCohort_pathvar$WDR19[EntireCohort_pathvar$GENE3 %in% c("WDR19")] <- 1
EntireCohort_pathvar$WDR19[EntireCohort_pathvar$GENE4 %in% c("WDR19")] <- 1
EntireCohort_pathvar$WDR19[EntireCohort_pathvar$GENE5 %in% c("WDR19")] <- 1
EntireCohort_pathvar$WDR19[EntireCohort_pathvar$GENE6 %in% c("WDR19")] <- 1
EntireCohort_pathvar$WDR19[EntireCohort_pathvar$GENE7 %in% c("WDR19")] <- 1
EntireCohort_pathvar$WDR19[EntireCohort_pathvar$GENE8 %in% c("WDR19")] <- 1
EntireCohort_pathvar$WDR19[EntireCohort_pathvar$GENE9 %in% c("WDR19")] <- 1

EntireCohort_pathvar$WDR19[is.na(EntireCohort_pathvar$WDR19)] <- 0
EntireCohort_pathvar$WDR19 <- as.numeric(EntireCohort_pathvar$WDR19)
unique(EntireCohort_pathvar$WDR19)
# 0 1

table(EntireCohort_pathvar$WDR19)
# 0      1 
# 465411   4388




### 'CHMP1A'
EntireCohort_pathvar$CHMP1A <- NA
EntireCohort_pathvar$CHMP1A[EntireCohort_pathvar$GENE1 %in% c("CHMP1A")] <- 1
EntireCohort_pathvar$CHMP1A[EntireCohort_pathvar$GENE2 %in% c("CHMP1A")] <- 1
EntireCohort_pathvar$CHMP1A[EntireCohort_pathvar$GENE3 %in% c("CHMP1A")] <- 1
EntireCohort_pathvar$CHMP1A[EntireCohort_pathvar$GENE4 %in% c("CHMP1A")] <- 1
EntireCohort_pathvar$CHMP1A[EntireCohort_pathvar$GENE5 %in% c("CHMP1A")] <- 1
EntireCohort_pathvar$CHMP1A[EntireCohort_pathvar$GENE6 %in% c("CHMP1A")] <- 1
EntireCohort_pathvar$CHMP1A[EntireCohort_pathvar$GENE7 %in% c("CHMP1A")] <- 1
EntireCohort_pathvar$CHMP1A[EntireCohort_pathvar$GENE8 %in% c("CHMP1A")] <- 1
EntireCohort_pathvar$CHMP1A[EntireCohort_pathvar$GENE9 %in% c("CHMP1A")] <- 1

EntireCohort_pathvar$CHMP1A[is.na(EntireCohort_pathvar$CHMP1A)] <- 0
EntireCohort_pathvar$CHMP1A <- as.numeric(EntireCohort_pathvar$CHMP1A)
unique(EntireCohort_pathvar$CHMP1A)
# 0 1

table(EntireCohort_pathvar$CHMP1A)
# 0      1 
# 469227    572




### 'ROR2'
EntireCohort_pathvar$ROR2 <- NA
EntireCohort_pathvar$ROR2[EntireCohort_pathvar$GENE1 %in% c("ROR2")] <- 1
EntireCohort_pathvar$ROR2[EntireCohort_pathvar$GENE2 %in% c("ROR2")] <- 1
EntireCohort_pathvar$ROR2[EntireCohort_pathvar$GENE3 %in% c("ROR2")] <- 1
EntireCohort_pathvar$ROR2[EntireCohort_pathvar$GENE4 %in% c("ROR2")] <- 1
EntireCohort_pathvar$ROR2[EntireCohort_pathvar$GENE5 %in% c("ROR2")] <- 1
EntireCohort_pathvar$ROR2[EntireCohort_pathvar$GENE6 %in% c("ROR2")] <- 1
EntireCohort_pathvar$ROR2[EntireCohort_pathvar$GENE7 %in% c("ROR2")] <- 1
EntireCohort_pathvar$ROR2[EntireCohort_pathvar$GENE8 %in% c("ROR2")] <- 1
EntireCohort_pathvar$ROR2[EntireCohort_pathvar$GENE9 %in% c("ROR2")] <- 1

EntireCohort_pathvar$ROR2[is.na(EntireCohort_pathvar$ROR2)] <- 0
EntireCohort_pathvar$ROR2 <- as.numeric(EntireCohort_pathvar$ROR2)
unique(EntireCohort_pathvar$ROR2)
# 0 1

table(EntireCohort_pathvar$ROR2)
# 0      1 
# 468028   1771




### 'DVL3'
# no variants in this gene




### 'EVC'
EntireCohort_pathvar$EVC <- NA
EntireCohort_pathvar$EVC[EntireCohort_pathvar$GENE1 %in% c("EVC")] <- 1
EntireCohort_pathvar$EVC[EntireCohort_pathvar$GENE2 %in% c("EVC")] <- 1
EntireCohort_pathvar$EVC[EntireCohort_pathvar$GENE3 %in% c("EVC")] <- 1
EntireCohort_pathvar$EVC[EntireCohort_pathvar$GENE4 %in% c("EVC")] <- 1
EntireCohort_pathvar$EVC[EntireCohort_pathvar$GENE5 %in% c("EVC")] <- 1
EntireCohort_pathvar$EVC[EntireCohort_pathvar$GENE6 %in% c("EVC")] <- 1
EntireCohort_pathvar$EVC[EntireCohort_pathvar$GENE7 %in% c("EVC")] <- 1
EntireCohort_pathvar$EVC[EntireCohort_pathvar$GENE8 %in% c("EVC")] <- 1
EntireCohort_pathvar$EVC[EntireCohort_pathvar$GENE9 %in% c("EVC")] <- 1

EntireCohort_pathvar$EVC[is.na(EntireCohort_pathvar$EVC)] <- 0
EntireCohort_pathvar$EVC <- as.numeric(EntireCohort_pathvar$EVC)
unique(EntireCohort_pathvar$EVC)
# 0 1

table(EntireCohort_pathvar$EVC)
#  0      1 
# 468424   1375




### 'EVC2'
EntireCohort_pathvar$EVC2 <- NA
EntireCohort_pathvar$EVC2[EntireCohort_pathvar$GENE1 %in% c("EVC2")] <- 1
EntireCohort_pathvar$EVC2[EntireCohort_pathvar$GENE2 %in% c("EVC2")] <- 1
EntireCohort_pathvar$EVC2[EntireCohort_pathvar$GENE3 %in% c("EVC2")] <- 1
EntireCohort_pathvar$EVC2[EntireCohort_pathvar$GENE4 %in% c("EVC2")] <- 1
EntireCohort_pathvar$EVC2[EntireCohort_pathvar$GENE5 %in% c("EVC2")] <- 1
EntireCohort_pathvar$EVC2[EntireCohort_pathvar$GENE6 %in% c("EVC2")] <- 1
EntireCohort_pathvar$EVC2[EntireCohort_pathvar$GENE7 %in% c("EVC2")] <- 1
EntireCohort_pathvar$EVC2[EntireCohort_pathvar$GENE8 %in% c("EVC2")] <- 1
EntireCohort_pathvar$EVC2[EntireCohort_pathvar$GENE9 %in% c("EVC2")] <- 1

EntireCohort_pathvar$EVC2[is.na(EntireCohort_pathvar$EVC2)] <- 0
EntireCohort_pathvar$EVC2 <- as.numeric(EntireCohort_pathvar$EVC2)
unique(EntireCohort_pathvar$EVC2)
# 0 1

table(EntireCohort_pathvar$EVC2)
# 0      1 
# 467384   2415




### 'ICK'
EntireCohort_pathvar$ICK <- NA
EntireCohort_pathvar$ICK[EntireCohort_pathvar$GENE1 %in% c("ICK")] <- 1
EntireCohort_pathvar$ICK[EntireCohort_pathvar$GENE2 %in% c("ICK")] <- 1
EntireCohort_pathvar$ICK[EntireCohort_pathvar$GENE3 %in% c("ICK")] <- 1
EntireCohort_pathvar$ICK[EntireCohort_pathvar$GENE4 %in% c("ICK")] <- 1
EntireCohort_pathvar$ICK[EntireCohort_pathvar$GENE5 %in% c("ICK")] <- 1
EntireCohort_pathvar$ICK[EntireCohort_pathvar$GENE6 %in% c("ICK")] <- 1
EntireCohort_pathvar$ICK[EntireCohort_pathvar$GENE7 %in% c("ICK")] <- 1
EntireCohort_pathvar$ICK[EntireCohort_pathvar$GENE8 %in% c("ICK")] <- 1
EntireCohort_pathvar$ICK[EntireCohort_pathvar$GENE9 %in% c("ICK")] <- 1

EntireCohort_pathvar$ICK[is.na(EntireCohort_pathvar$ICK)] <- 0
EntireCohort_pathvar$ICK <- as.numeric(EntireCohort_pathvar$ICK)
unique(EntireCohort_pathvar$ICK)
# 0 1

table(EntireCohort_pathvar$ICK)
# 0      1 
# 469256    543





### 'IFT122'
EntireCohort_pathvar$IFT122 <- NA
EntireCohort_pathvar$IFT122[EntireCohort_pathvar$GENE1 %in% c("IFT122")] <- 1
EntireCohort_pathvar$IFT122[EntireCohort_pathvar$GENE2 %in% c("IFT122")] <- 1
EntireCohort_pathvar$IFT122[EntireCohort_pathvar$GENE3 %in% c("IFT122")] <- 1
EntireCohort_pathvar$IFT122[EntireCohort_pathvar$GENE4 %in% c("IFT122")] <- 1
EntireCohort_pathvar$IFT122[EntireCohort_pathvar$GENE5 %in% c("IFT122")] <- 1
EntireCohort_pathvar$IFT122[EntireCohort_pathvar$GENE6 %in% c("IFT122")] <- 1
EntireCohort_pathvar$IFT122[EntireCohort_pathvar$GENE7 %in% c("IFT122")] <- 1
EntireCohort_pathvar$IFT122[EntireCohort_pathvar$GENE8 %in% c("IFT122")] <- 1
EntireCohort_pathvar$IFT122[EntireCohort_pathvar$GENE9 %in% c("IFT122")] <- 1

EntireCohort_pathvar$IFT122[is.na(EntireCohort_pathvar$IFT122)] <- 0
EntireCohort_pathvar$IFT122 <- as.numeric(EntireCohort_pathvar$IFT122)
unique(EntireCohort_pathvar$IFT122)
# 0 1

table(EntireCohort_pathvar$IFT122)
#  0      1 
# 464948   4851




### 'WDR35'
EntireCohort_pathvar$WDR35 <- NA
EntireCohort_pathvar$WDR35[EntireCohort_pathvar$GENE1 %in% c("WDR35")] <- 1
EntireCohort_pathvar$WDR35[EntireCohort_pathvar$GENE2 %in% c("WDR35")] <- 1
EntireCohort_pathvar$WDR35[EntireCohort_pathvar$GENE3 %in% c("WDR35")] <- 1
EntireCohort_pathvar$WDR35[EntireCohort_pathvar$GENE4 %in% c("WDR35")] <- 1
EntireCohort_pathvar$WDR35[EntireCohort_pathvar$GENE5 %in% c("WDR35")] <- 1
EntireCohort_pathvar$WDR35[EntireCohort_pathvar$GENE6 %in% c("WDR35")] <- 1
EntireCohort_pathvar$WDR35[EntireCohort_pathvar$GENE7 %in% c("WDR35")] <- 1
EntireCohort_pathvar$WDR35[EntireCohort_pathvar$GENE8 %in% c("WDR35")] <- 1
EntireCohort_pathvar$WDR35[EntireCohort_pathvar$GENE9 %in% c("WDR35")] <- 1

EntireCohort_pathvar$WDR35[is.na(EntireCohort_pathvar$WDR35)] <- 0
EntireCohort_pathvar$WDR35 <- as.numeric(EntireCohort_pathvar$WDR35)
unique(EntireCohort_pathvar$WDR35)
# 0 1

table(EntireCohort_pathvar$WDR35)
# 0      1 
# 464981   4818 




### 'GLI1'
EntireCohort_pathvar$GLI1 <- NA
EntireCohort_pathvar$GLI1[EntireCohort_pathvar$GENE1 %in% c("GLI1")] <- 1
EntireCohort_pathvar$GLI1[EntireCohort_pathvar$GENE2 %in% c("GLI1")] <- 1
EntireCohort_pathvar$GLI1[EntireCohort_pathvar$GENE3 %in% c("GLI1")] <- 1
EntireCohort_pathvar$GLI1[EntireCohort_pathvar$GENE4 %in% c("GLI1")] <- 1
EntireCohort_pathvar$GLI1[EntireCohort_pathvar$GENE5 %in% c("GLI1")] <- 1
EntireCohort_pathvar$GLI1[EntireCohort_pathvar$GENE6 %in% c("GLI1")] <- 1
EntireCohort_pathvar$GLI1[EntireCohort_pathvar$GENE7 %in% c("GLI1")] <- 1
EntireCohort_pathvar$GLI1[EntireCohort_pathvar$GENE8 %in% c("GLI1")] <- 1
EntireCohort_pathvar$GLI1[EntireCohort_pathvar$GENE9 %in% c("GLI1")] <- 1

EntireCohort_pathvar$GLI1[is.na(EntireCohort_pathvar$GLI1)] <- 0
EntireCohort_pathvar$GLI1 <- as.numeric(EntireCohort_pathvar$GLI1)
unique(EntireCohort_pathvar$GLI1)
# 0 1

table(EntireCohort_pathvar$GLI1)
# 0      1 
# 468537   1262




### 'GLI2'
EntireCohort_pathvar$GLI2 <- NA
EntireCohort_pathvar$GLI2[EntireCohort_pathvar$GENE1 %in% c("GLI2")] <- 1
EntireCohort_pathvar$GLI2[EntireCohort_pathvar$GENE2 %in% c("GLI2")] <- 1
EntireCohort_pathvar$GLI2[EntireCohort_pathvar$GENE3 %in% c("GLI2")] <- 1
EntireCohort_pathvar$GLI2[EntireCohort_pathvar$GENE4 %in% c("GLI2")] <- 1
EntireCohort_pathvar$GLI2[EntireCohort_pathvar$GENE5 %in% c("GLI2")] <- 1
EntireCohort_pathvar$GLI2[EntireCohort_pathvar$GENE6 %in% c("GLI2")] <- 1
EntireCohort_pathvar$GLI2[EntireCohort_pathvar$GENE7 %in% c("GLI2")] <- 1
EntireCohort_pathvar$GLI2[EntireCohort_pathvar$GENE8 %in% c("GLI2")] <- 1
EntireCohort_pathvar$GLI2[EntireCohort_pathvar$GENE9 %in% c("GLI2")] <- 1

EntireCohort_pathvar$GLI2[is.na(EntireCohort_pathvar$GLI2)] <- 0
EntireCohort_pathvar$GLI2 <- as.numeric(EntireCohort_pathvar$GLI2)
unique(EntireCohort_pathvar$GLI2)
# 0 1

table(EntireCohort_pathvar$GLI2)
# 0      1 
# 466229   3570 




### 'GLI3'
EntireCohort_pathvar$GLI3 <- NA
EntireCohort_pathvar$GLI3[EntireCohort_pathvar$GENE1 %in% c("GLI3")] <- 1
EntireCohort_pathvar$GLI3[EntireCohort_pathvar$GENE2 %in% c("GLI3")] <- 1
EntireCohort_pathvar$GLI3[EntireCohort_pathvar$GENE3 %in% c("GLI3")] <- 1
EntireCohort_pathvar$GLI3[EntireCohort_pathvar$GENE4 %in% c("GLI3")] <- 1
EntireCohort_pathvar$GLI3[EntireCohort_pathvar$GENE5 %in% c("GLI3")] <- 1
EntireCohort_pathvar$GLI3[EntireCohort_pathvar$GENE6 %in% c("GLI3")] <- 1
EntireCohort_pathvar$GLI3[EntireCohort_pathvar$GENE7 %in% c("GLI3")] <- 1
EntireCohort_pathvar$GLI3[EntireCohort_pathvar$GENE8 %in% c("GLI3")] <- 1
EntireCohort_pathvar$GLI3[EntireCohort_pathvar$GENE9 %in% c("GLI3")] <- 1

EntireCohort_pathvar$GLI3[is.na(EntireCohort_pathvar$GLI3)] <- 0
EntireCohort_pathvar$GLI3 <- as.numeric(EntireCohort_pathvar$GLI3)
unique(EntireCohort_pathvar$GLI3)
# 0 1

table(EntireCohort_pathvar$GLI3)
# 0      1 
# 467986   1813




### 'FZD2'
# no variants in this gene





### 'HHAT'
EntireCohort_pathvar$HHAT <- NA
EntireCohort_pathvar$HHAT[EntireCohort_pathvar$GENE1 %in% c("HHAT")] <- 1
EntireCohort_pathvar$HHAT[EntireCohort_pathvar$GENE2 %in% c("HHAT")] <- 1
EntireCohort_pathvar$HHAT[EntireCohort_pathvar$GENE3 %in% c("HHAT")] <- 1
EntireCohort_pathvar$HHAT[EntireCohort_pathvar$GENE4 %in% c("HHAT")] <- 1
EntireCohort_pathvar$HHAT[EntireCohort_pathvar$GENE5 %in% c("HHAT")] <- 1
EntireCohort_pathvar$HHAT[EntireCohort_pathvar$GENE6 %in% c("HHAT")] <- 1
EntireCohort_pathvar$HHAT[EntireCohort_pathvar$GENE7 %in% c("HHAT")] <- 1
EntireCohort_pathvar$HHAT[EntireCohort_pathvar$GENE8 %in% c("HHAT")] <- 1
EntireCohort_pathvar$HHAT[EntireCohort_pathvar$GENE9 %in% c("HHAT")] <- 1

EntireCohort_pathvar$HHAT[is.na(EntireCohort_pathvar$HHAT)] <- 0
EntireCohort_pathvar$HHAT <- as.numeric(EntireCohort_pathvar$HHAT)
unique(EntireCohort_pathvar$HHAT)
# 0 1

table(EntireCohort_pathvar$HHAT)
#  0      1 
# 468769   1030





### 'SHH'
EntireCohort_pathvar$SHH <- NA
EntireCohort_pathvar$SHH[EntireCohort_pathvar$GENE1 %in% c("SHH")] <- 1
EntireCohort_pathvar$SHH[EntireCohort_pathvar$GENE2 %in% c("SHH")] <- 1
EntireCohort_pathvar$SHH[EntireCohort_pathvar$GENE3 %in% c("SHH")] <- 1
EntireCohort_pathvar$SHH[EntireCohort_pathvar$GENE4 %in% c("SHH")] <- 1
EntireCohort_pathvar$SHH[EntireCohort_pathvar$GENE5 %in% c("SHH")] <- 1
EntireCohort_pathvar$SHH[EntireCohort_pathvar$GENE6 %in% c("SHH")] <- 1
EntireCohort_pathvar$SHH[EntireCohort_pathvar$GENE7 %in% c("SHH")] <- 1
EntireCohort_pathvar$SHH[EntireCohort_pathvar$GENE8 %in% c("SHH")] <- 1
EntireCohort_pathvar$SHH[EntireCohort_pathvar$GENE9 %in% c("SHH")] <- 1

EntireCohort_pathvar$SHH[is.na(EntireCohort_pathvar$SHH)] <- 0
EntireCohort_pathvar$SHH <- as.numeric(EntireCohort_pathvar$SHH)
unique(EntireCohort_pathvar$SHH)
# 0 1

table(EntireCohort_pathvar$SHH)
#  0      1 
# 469497    302




### 'MYCN'
EntireCohort_pathvar$MYCN <- NA
EntireCohort_pathvar$MYCN[EntireCohort_pathvar$GENE1 %in% c("MYCN")] <- 1
EntireCohort_pathvar$MYCN[EntireCohort_pathvar$GENE2 %in% c("MYCN")] <- 1
EntireCohort_pathvar$MYCN[EntireCohort_pathvar$GENE3 %in% c("MYCN")] <- 1
EntireCohort_pathvar$MYCN[EntireCohort_pathvar$GENE4 %in% c("MYCN")] <- 1
EntireCohort_pathvar$MYCN[EntireCohort_pathvar$GENE5 %in% c("MYCN")] <- 1
EntireCohort_pathvar$MYCN[EntireCohort_pathvar$GENE6 %in% c("MYCN")] <- 1
EntireCohort_pathvar$MYCN[EntireCohort_pathvar$GENE7 %in% c("MYCN")] <- 1
EntireCohort_pathvar$MYCN[EntireCohort_pathvar$GENE8 %in% c("MYCN")] <- 1
EntireCohort_pathvar$MYCN[EntireCohort_pathvar$GENE9 %in% c("MYCN")] <- 1

EntireCohort_pathvar$MYCN[is.na(EntireCohort_pathvar$MYCN)] <- 0
EntireCohort_pathvar$MYCN <- as.numeric(EntireCohort_pathvar$MYCN)
unique(EntireCohort_pathvar$MYCN)
# 0 1

table(EntireCohort_pathvar$MYCN)
# 0      1 
# 469359    440




### 'WNT1'
EntireCohort_pathvar$WNT1 <- NA
EntireCohort_pathvar$WNT1[EntireCohort_pathvar$GENE1 %in% c("WNT1")] <- 1
EntireCohort_pathvar$WNT1[EntireCohort_pathvar$GENE2 %in% c("WNT1")] <- 1
EntireCohort_pathvar$WNT1[EntireCohort_pathvar$GENE3 %in% c("WNT1")] <- 1
EntireCohort_pathvar$WNT1[EntireCohort_pathvar$GENE4 %in% c("WNT1")] <- 1
EntireCohort_pathvar$WNT1[EntireCohort_pathvar$GENE5 %in% c("WNT1")] <- 1
EntireCohort_pathvar$WNT1[EntireCohort_pathvar$GENE6 %in% c("WNT1")] <- 1
EntireCohort_pathvar$WNT1[EntireCohort_pathvar$GENE7 %in% c("WNT1")] <- 1
EntireCohort_pathvar$WNT1[EntireCohort_pathvar$GENE8 %in% c("WNT1")] <- 1
EntireCohort_pathvar$WNT1[EntireCohort_pathvar$GENE9 %in% c("WNT1")] <- 1

EntireCohort_pathvar$WNT1[is.na(EntireCohort_pathvar$WNT1)] <- 0
EntireCohort_pathvar$WNT1 <- as.numeric(EntireCohort_pathvar$WNT1)
unique(EntireCohort_pathvar$WNT1)
# 0 1

table(EntireCohort_pathvar$WNT1)
#  0      1 
# 468226   1573




### 'WNT5A'
EntireCohort_pathvar$WNT5A <- NA
EntireCohort_pathvar$WNT5A[EntireCohort_pathvar$GENE1 %in% c("WNT5A")] <- 1
EntireCohort_pathvar$WNT5A[EntireCohort_pathvar$GENE2 %in% c("WNT5A")] <- 1
EntireCohort_pathvar$WNT5A[EntireCohort_pathvar$GENE3 %in% c("WNT5A")] <- 1
EntireCohort_pathvar$WNT5A[EntireCohort_pathvar$GENE4 %in% c("WNT5A")] <- 1
EntireCohort_pathvar$WNT5A[EntireCohort_pathvar$GENE5 %in% c("WNT5A")] <- 1
EntireCohort_pathvar$WNT5A[EntireCohort_pathvar$GENE6 %in% c("WNT5A")] <- 1
EntireCohort_pathvar$WNT5A[EntireCohort_pathvar$GENE7 %in% c("WNT5A")] <- 1
EntireCohort_pathvar$WNT5A[EntireCohort_pathvar$GENE8 %in% c("WNT5A")] <- 1
EntireCohort_pathvar$WNT5A[EntireCohort_pathvar$GENE9 %in% c("WNT5A")] <- 1

EntireCohort_pathvar$WNT5A[is.na(EntireCohort_pathvar$WNT5A)] <- 0
EntireCohort_pathvar$WNT5A <- as.numeric(EntireCohort_pathvar$WNT5A)
unique(EntireCohort_pathvar$WNT5A)
# 0 1

table(EntireCohort_pathvar$WNT5A)
# 0      1 
# 469306    493




### 'NXN'
EntireCohort_pathvar$NXN <- NA
EntireCohort_pathvar$NXN[EntireCohort_pathvar$GENE1 %in% c("NXN")] <- 1
EntireCohort_pathvar$NXN[EntireCohort_pathvar$GENE2 %in% c("NXN")] <- 1
EntireCohort_pathvar$NXN[EntireCohort_pathvar$GENE3 %in% c("NXN")] <- 1
EntireCohort_pathvar$NXN[EntireCohort_pathvar$GENE4 %in% c("NXN")] <- 1
EntireCohort_pathvar$NXN[EntireCohort_pathvar$GENE5 %in% c("NXN")] <- 1
EntireCohort_pathvar$NXN[EntireCohort_pathvar$GENE6 %in% c("NXN")] <- 1
EntireCohort_pathvar$NXN[EntireCohort_pathvar$GENE7 %in% c("NXN")] <- 1
EntireCohort_pathvar$NXN[EntireCohort_pathvar$GENE8 %in% c("NXN")] <- 1
EntireCohort_pathvar$NXN[EntireCohort_pathvar$GENE9 %in% c("NXN")] <- 1

EntireCohort_pathvar$NXN[is.na(EntireCohort_pathvar$NXN)] <- 0
EntireCohort_pathvar$NXN <- as.numeric(EntireCohort_pathvar$NXN)
unique(EntireCohort_pathvar$NXN)
# 0 1

table(EntireCohort_pathvar$NXN)
# 0      1 
# 468775   1024




### 'WDFY3'
# no variants in this gene



# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)

















# WNTSHHgene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_WNTSHHgene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_WNTSHHgene_, "pathvar_WNTSHHgene_.csv", row.names=FALSE)
# system('dx upload pathvar_WNTSHHgene_.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_.csv')


# WNTSHHgene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ WNTSHHgene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_WNTSHHgene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_incidence_1, "pathvar_WNTSHHgene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_incidence_1.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_incidence_1.csv')


# GOIs
LR_2 <- glm(formula = cancer_ROOC_dummy ~ IHH+GPC6+WDR19+CHMP1A+ROR2+EVC+EVC2+ICK+IFT122+WDR35+GLI1+GLI2+GLI3+HHAT+SHH+MYCN+WNT1+WNT5A+NXN, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_WNTSHHgene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_incidence_2, "pathvar_WNTSHHgene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_incidence_2.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_incidence_2.csv')







### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "WNTSHHgene", 'IHH', 'GPC6', 'WDR19', 'CHMP1A', 'ROR2', 
                                     'EVC', 'EVC2', 'ICK', 'IFT122', 'WDR35', 'GLI1', 
                                     'GLI2', 'GLI3', 'HHAT', 'SHH', 'MYCN', 'WNT1', 
                                     'WNT5A', 'NXN')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
#465881     26

with(data_setA, table(NXN))
# all have plenty of participants


# WNTSHHgene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ WNTSHHgene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_WNTSHHgene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_incidence_3, "pathvar_WNTSHHgene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_incidence_3.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_incidence_3.csv')


# age + ethnicity + smoking + sex + WNTSHHgene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_WNTSHHgene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_incidence_4, "pathvar_WNTSHHgene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_incidence_4.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ IHH+GPC6+WDR19+CHMP1A+ROR2+EVC+EVC2+ICK+IFT122+WDR35+GLI1+GLI2+GLI3+HHAT+SHH+MYCN+WNT1+WNT5A+NXN, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_WNTSHHgene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_incidence_5, "pathvar_WNTSHHgene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_incidence_5.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_incidence_5.csv')


# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+IHH+GPC6+WDR19+CHMP1A+ROR2+EVC+EVC2+ICK+IFT122+WDR35+GLI1+GLI2+GLI3+HHAT+SHH+MYCN+WNT1+WNT5A+NXN, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_WNTSHHgene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_incidence_6, "pathvar_WNTSHHgene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_incidence_6.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_incidence_6.csv')




### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "WNTSHHgene", 'IHH', 'GPC6', 'WDR19', 'CHMP1A', 
                                     'ROR2', 'EVC', 'EVC2', 'ICK', 'IFT122', 
                                     'WDR35', 'GLI1', 'GLI2', 'GLI3', 'HHAT', 
                                     'SHH', 'MYCN', 'WNT1', 'WNT5A', 'NXN')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     31

with(data_setB, table(NXN))
# looking good


# WNTSHHgene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ WNTSHHgene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_WNTSHHgene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_incidence_7, "pathvar_WNTSHHgene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_incidence_7.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_incidence_7.csv')


# all variables + WNTSHHgene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+WNTSHHgene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_WNTSHHgene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_incidence_8, "pathvar_WNTSHHgene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_incidence_8.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_incidence_8.csv')


# GOIs alone
LR_9 <- glm(formula = cancer_ROOC_dummy ~ IHH+GPC6+WDR19+CHMP1A+ROR2+EVC+EVC2+ICK+IFT122+WDR35+GLI1+GLI2+GLI3+HHAT+SHH+MYCN+WNT1+WNT5A+NXN, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_WNTSHHgene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_incidence_9, "pathvar_WNTSHHgene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_incidence_9.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_incidence_9.csv')


# all variables + GOIs
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+IHH+GPC6+WDR19+CHMP1A+ROR2+EVC+EVC2+ICK+IFT122+WDR35+GLI1+GLI2+GLI3+HHAT+SHH+MYCN+WNT1+WNT5A+NXN, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_WNTSHHgene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_incidence_10, "pathvar_WNTSHHgene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_incidence_10.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_incidence_10.csv')










### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(WNTSHHgene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(WNTSHHgene, cancer_ROOC_dummy)))
# p-value = 0.04979
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9480613 1.0000000
# sample estimates:
#   odds ratio 
# 0.9737221

with(EntireCohort_pathvar, table(WNTSHHgene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(WNTSHHgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5194773 0.5634993
# sample estimates:
#   odds ratio 
# 0.5409295

with(EntireCohort_pathvar, table(WNTSHHgene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(WNTSHHgene, Ever_smoked_dummy)))
# p-value = 0.01143
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9488075 0.9933704
# sample estimates:
#   odds ratio 
# 0.9707981

with(EntireCohort_pathvar, table(WNTSHHgene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(WNTSHHgene, Sex_Male)))
# p-value = 0.6376

with(EntireCohort_pathvar, table(WNTSHHgene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(WNTSHHgene, University_education)))
# p-value = 0.6943

with(EntireCohort_pathvar, table(WNTSHHgene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(WNTSHHgene, Annual_income)))
# p-value = 0.03744
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9437004 0.9982904
# sample estimates:
#   odds ratio 
# 0.9706548

with(EntireCohort_pathvar, table(WNTSHHgene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(WNTSHHgene, Smoker_previous)))
# p-value = 0.1565

with(EntireCohort_pathvar, table(WNTSHHgene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(WNTSHHgene, Smoker_current)))
# p-value = 0.6406



### data_setA
with(data_setA, table(WNTSHHgene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(WNTSHHgene, cancer_ROOC_dummy)))
# p-value = 0.05195

with(data_setA, table(WNTSHHgene, Ethnicity_White))			
fisher.test(with(data_setA, table(WNTSHHgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5187480 0.5629275
# sample estimates:
#   odds ratio 
# 0.5402981

with(data_setA, table(WNTSHHgene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(WNTSHHgene, Ever_smoked_dummy)))
# p-value = 0.01493
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9497728 0.9945025
# sample estimates:
#   odds ratio 
# 0.9718634

with(data_setA, table(WNTSHHgene, Sex_Male))			
fisher.test(with(data_setA, table(WNTSHHgene, Sex_Male)))
# p-value = 0.5522



### data_setB
with(data_setB, table(WNTSHHgene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(WNTSHHgene, cancer_ROOC_dummy)))
# p-value = 0.1842

with(data_setB, table(WNTSHHgene, Ethnicity_White))			
fisher.test(with(data_setB, table(WNTSHHgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5170027 0.5702665
# sample estimates:
#   odds ratio 
# 0.5428731

with(data_setB, table(WNTSHHgene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(WNTSHHgene, Ever_smoked_dummy)))
# p-value = 0.1555

with(data_setB, table(WNTSHHgene, Sex_Male))			
fisher.test(with(data_setB, table(WNTSHHgene, Sex_Male)))
# p-value = 1

with(data_setB, table(WNTSHHgene, University_education))			
fisher.test(with(data_setB, table(WNTSHHgene, University_education)))
# p-value = 0.9101

with(data_setB, table(WNTSHHgene, Annual_income))			
fisher.test(with(data_setB, table(WNTSHHgene, Annual_income)))
# p-value = 0.04873
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9443418 0.9999011
# sample estimates:
#   odds ratio 
# 0.9717549 






### can plot the continuous variables: 

# make WNTSHHgene into factor
df_for_plot <- EntireCohort_pathvar
df_for_plot$WNTSHHgene <- as.factor(df_for_plot$WNTSHHgene)

# try a violin plot

# Age_at_recruitment  
df_for_plot %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ WNTSHHgene, data = df_for_plot)
# W = 7.273e+09, p-value = 6.93e-06

# BMI  
df_for_plot %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = BMI)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ WNTSHHgene, data = df_for_plot)
# W = 7035945218, p-value = 0.003661

# Standing_height 
df_for_plot %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ WNTSHHgene, data = df_for_plot)
# W = 7215211828, p-value = 6.174e-05

# make this a bar plot
# Weekly_exercise
df_for_plot %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(WNTSHHgene), scales = "free_y") +
  theme_light()
# 



### data_setA

# make WNTSHHgene into factor
data_setA$WNTSHHgene <- as.factor(data_setA$WNTSHHgene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ WNTSHHgene, data = data_setA)
# W = 7147105278, p-value = 9.553e-06



### data_setB

# make WNTSHHgene into factor
data_setB$WNTSHHgene <- as.factor(data_setB$WNTSHHgene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ WNTSHHgene, data = data_setB)
# W = 4.754e+09, p-value = 0.0005877

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ WNTSHHgene, data = data_setB)
# W = 4657286572, p-value = 0.03056

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ WNTSHHgene, data = data_setB)
# W = 4760102428, p-value = 0.0001492


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(WNTSHHgene), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# WNTSHHgene Type ------------------------------------------------------------

# make sure the execute the code in the 'WNTSHHgenes' section

## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_WNTSHHgene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_type_1, "pathvar_WNTSHHgene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_type_1.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_WNTSHHgene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_type_2, "pathvar_WNTSHHgene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_type_2.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_WNTSHHgene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_type_3, "pathvar_WNTSHHgene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_type_3.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_WNTSHHgene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_type_4, "pathvar_WNTSHHgene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_type_4.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_WNTSHHgene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_type_5, "pathvar_WNTSHHgene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_type_5.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_WNTSHHgene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_type_6, "pathvar_WNTSHHgene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_type_6.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_WNTSHHgene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_type_7, "pathvar_WNTSHHgene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_type_7.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_WNTSHHgene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_type_8, "pathvar_WNTSHHgene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_type_8.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_WNTSHHgene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_type_9, "pathvar_WNTSHHgene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_type_9.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_WNTSHHgene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_type_10, "pathvar_WNTSHHgene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_type_10.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_WNTSHHgene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_type_11, "pathvar_WNTSHHgene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_type_11.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_WNTSHHgene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_type_12, "pathvar_WNTSHHgene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_type_12.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_WNTSHHgene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_type_13, "pathvar_WNTSHHgene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_type_13.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_WNTSHHgene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_type_14, "pathvar_WNTSHHgene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_type_14.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_WNTSHHgene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_WNTSHHgene_type_15, "pathvar_WNTSHHgene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_type_15.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_type_15.csv')















# WNTSHHgene Age -------------------------------------------------------------

# make sure the execute the code in the 'WNTSHHgenes' section

### linear regression
# "WNTSHHgene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ WNTSHHgene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_WNTSHHgene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_WNTSHHgene_age_1, "pathvar_WNTSHHgene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_age_1.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ IHH+GPC6+WDR19+CHMP1A+ROR2+EVC+EVC2+ICK+IFT122+WDR35+GLI1+GLI2+GLI3+HHAT+SHH+MYCN+WNT1+WNT5A+NXN, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_WNTSHHgene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_WNTSHHgene_age_2, "pathvar_WNTSHHgene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_age_2.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "WNTSHHgene", 'IHH', 'GPC6', 'WDR19', 'CHMP1A', 'ROR2', 
                                     'EVC', 'EVC2', 'ICK', 'IFT122', 'WDR35', 'GLI1', 
                                     'GLI2', 'GLI3', 'HHAT', 'SHH', 'MYCN', 'WNT1', 
                                     'WNT5A', 'NXN')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     26




data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "WNTSHHgene", 'IHH', 'GPC6', 'WDR19', 
                                     'CHMP1A', 'ROR2', 'EVC', 'EVC2', 'ICK', 'IFT122', 
                                     'WDR35', 'GLI1', 'GLI2', 'GLI3', 'HHAT', 'SHH', 
                                     'MYCN', 'WNT1', 'WNT5A', 'NXN')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    31


### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ WNTSHHgene, data = data_setD)
# summary(LM_15)
pathvar_WNTSHHgene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_WNTSHHgene_age_3, "pathvar_WNTSHHgene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_age_3.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+WNTSHHgene, data = data_setD)
# summary(LM_16)
pathvar_WNTSHHgene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_WNTSHHgene_age_4, "pathvar_WNTSHHgene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_age_4.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ IHH+GPC6+WDR19+CHMP1A+ROR2+EVC+EVC2+ICK+IFT122+WDR35+GLI1+GLI2+GLI3+HHAT+SHH+MYCN+WNT1+WNT5A+NXN, data = data_setD)
# summary(LM_17)
pathvar_WNTSHHgene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_WNTSHHgene_age_5, "pathvar_WNTSHHgene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_age_5.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+IHH+GPC6+WDR19+CHMP1A+ROR2+EVC+EVC2+ICK+IFT122+WDR35+GLI1+GLI2+GLI3+HHAT+SHH+MYCN+WNT1+WNT5A+NXN, data = data_setD)
# summary(LM_18)
pathvar_WNTSHHgene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_WNTSHHgene_age_6, "pathvar_WNTSHHgene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_age_6.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ WNTSHHgene, data = data_setE)
# summary(LM_15)
pathvar_WNTSHHgene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_WNTSHHgene_age_7, "pathvar_WNTSHHgene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_age_7.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+WNTSHHgene, data = data_setE)
# summary(LM_16)
pathvar_WNTSHHgene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_WNTSHHgene_age_8, "pathvar_WNTSHHgene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_age_8.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ IHH+GPC6+WDR19+CHMP1A+ROR2+EVC+EVC2+ICK+IFT122+WDR35+GLI1+GLI2+GLI3+HHAT+SHH+MYCN+WNT1+WNT5A+NXN, data = data_setE)
# summary(LM_17)
pathvar_WNTSHHgene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_WNTSHHgene_age_9, "pathvar_WNTSHHgene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_age_9.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+IHH+GPC6+WDR19+CHMP1A+ROR2+EVC+EVC2+ICK+IFT122+WDR35+GLI1+GLI2+GLI3+HHAT+SHH+MYCN+WNT1+WNT5A+NXN, data = data_setE)
# summary(LM_18)
pathvar_WNTSHHgene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_WNTSHHgene_age_10, "pathvar_WNTSHHgene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_WNTSHHgene_age_10.csv --path Emily-folder/results_2/pathvar_WNTSHHgene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diganosis <- EntireCohort_pathvar
df_cancer_diganosis <- df_cancer_diganosis[complete.cases(df_cancer_diganosis$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(df_cancer_diganosis, table(WNTSHHgene, Ethnicity_White))			
fisher.test(with(df_cancer_diganosis, table(WNTSHHgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.4536508 0.5690805
# sample estimates:
#   odds ratio 
# 0.5075263

with(df_cancer_diganosis, table(WNTSHHgene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diganosis, table(WNTSHHgene, Ever_smoked_dummy)))
# p-value = 0.9409

with(df_cancer_diganosis, table(WNTSHHgene, Sex_Male))			
fisher.test(with(df_cancer_diganosis, table(WNTSHHgene, Sex_Male)))
# p-value = 0.8393

with(df_cancer_diganosis, table(WNTSHHgene, University_education))			
fisher.test(with(df_cancer_diganosis, table(WNTSHHgene, University_education)))
# p-value = 0.8655

with(df_cancer_diganosis, table(WNTSHHgene, Annual_income))			
fisher.test(with(df_cancer_diganosis, table(WNTSHHgene, Annual_income)))
# p-value = 0.3799


### data_setD
with(data_setD, table(WNTSHHgene, Ethnicity_White))			
fisher.test(with(data_setD, table(WNTSHHgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.452867 0.568675
# sample estimates:
#   odds ratio 
# 0.5069109

with(data_setD, table(WNTSHHgene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(WNTSHHgene, Ever_smoked_dummy)))
# p-value = 0.9309

with(data_setD, table(WNTSHHgene, Sex_Male))			
fisher.test(with(data_setD, table(WNTSHHgene, Sex_Male)))
# p-value = 0.8856


### data_setE
with(data_setE, table(WNTSHHgene, Ethnicity_White))			
fisher.test(with(data_setE, table(WNTSHHgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.4557069 0.6029668
# sample estimates:
#   odds ratio 
# 0.5232762

with(data_setE, table(WNTSHHgene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(WNTSHHgene, Ever_smoked_dummy)))
# p-value = 0.8248

with(data_setE, table(WNTSHHgene, Sex_Male))			
fisher.test(with(data_setE, table(WNTSHHgene, Sex_Male)))
# p-value = 0.3864

with(data_setE, table(WNTSHHgene, University_education))			
fisher.test(with(data_setE, table(WNTSHHgene, University_education)))
# p-value = 0.9775

with(data_setE, table(WNTSHHgene, Annual_income))			
fisher.test(with(data_setE, table(WNTSHHgene, Annual_income)))
# p-value = 0.3989



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 
df_for_plot <- EntireCohort_pathvar
df_for_plot <- df_for_plot[complete.cases(df_for_plot$Age_at_cancer_diagnosis_earliest), ]
# 110872

# make WNTSHHgene into factor
df_for_plot$WNTSHHgene <- as.factor(df_for_plot$WNTSHHgene)

# try a violin plot



# Age_at_cancer_diagnosis_earliest
df_for_plot %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_for_plot %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ WNTSHHgene, data = df_for_plot)
# 



# Age_at_recruitment  
df_for_plot %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_for_plot %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ WNTSHHgene, data = df_for_plot)
# W = 397221054, p-value = 0.06156



# BMI  
df_for_plot %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = BMI)) +
  geom_violin() +
  theme_light()

df_for_plot %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ WNTSHHgene, data = df_for_plot)
# W = 389479084, p-value = 0.9488



# Standing_height 
df_for_plot %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_for_plot %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ WNTSHHgene, data = df_for_plot)
# W = 395876644, p-value = 0.02914



# make this a bar plot
# Weekly_exercise
df_for_plot %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(WNTSHHgene), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$WNTSHHgene <- as.factor(data_setD$WNTSHHgene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ WNTSHHgene, data = data_setD)
# 



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ WNTSHHgene, data = data_setD)
# W = 390431536, p-value = 0.05882




### data_setE

data_setE$WNTSHHgene <- as.factor(data_setE$WNTSHHgene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ WNTSHHgene, data = data_setE)
# 



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ WNTSHHgene, data = data_setE)
# W = 253546570, p-value = 0.03698



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ WNTSHHgene, data = data_setE)
# W = 250034266, p-value = 0.8005



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = WNTSHHgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ WNTSHHgene, data = data_setE)
# W = 255473828, p-value = 0.002003


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(WNTSHHgene), scales = "free_y") +
  theme_light()
# 

















# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_WNTSHHgenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_WNTSHHgenes.R')
