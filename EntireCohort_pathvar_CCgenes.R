# 24/06/2024

# EntireCohort_pathvar_CCgenes.R

# this script is for the Cell cycle associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)


# CCgene
# PDCD6IP
# AMMECR1
# CDKN1C
# MCPH1
# NCAPD2
# NCAPD3
# NCAPH
# VRK1
# ESCO2
# ANKLE2
# DYRK1A
# ATR
# ATRIP
# DNA2
# FAM111A
# PRIM1
# PHC1
# KIF14







# CCgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')


# go thr all the genes to check that they're present
grep("KIF14", EntireCohort_pathvar$GENE)
# all seem to be present, except for:
# CDKN1C = 0
# FAM111A = 0




# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$CCgene <- NA
EntireCohort_pathvar$CCgene[EntireCohort_pathvar$GENE1 %in% c('PDCD6IP', 'AMMECR1', 'MCPH1', 'NCAPD2', 'NCAPD3', 'NCAPH', 'VRK1', 'ESCO2', 'ANKLE2', 'DYRK1A', 'ATR', 'ATRIP', 'DNA2', 'PRIM1', 'PHC1', 'KIF14')] <- 1
EntireCohort_pathvar$CCgene[EntireCohort_pathvar$GENE2 %in% c('PDCD6IP', 'AMMECR1', 'MCPH1', 'NCAPD2', 'NCAPD3', 'NCAPH', 'VRK1', 'ESCO2', 'ANKLE2', 'DYRK1A', 'ATR', 'ATRIP', 'DNA2', 'PRIM1', 'PHC1', 'KIF14')] <- 1
EntireCohort_pathvar$CCgene[EntireCohort_pathvar$GENE3 %in% c('PDCD6IP', 'AMMECR1', 'MCPH1', 'NCAPD2', 'NCAPD3', 'NCAPH', 'VRK1', 'ESCO2', 'ANKLE2', 'DYRK1A', 'ATR', 'ATRIP', 'DNA2', 'PRIM1', 'PHC1', 'KIF14')] <- 1
EntireCohort_pathvar$CCgene[EntireCohort_pathvar$GENE4 %in% c('PDCD6IP', 'AMMECR1', 'MCPH1', 'NCAPD2', 'NCAPD3', 'NCAPH', 'VRK1', 'ESCO2', 'ANKLE2', 'DYRK1A', 'ATR', 'ATRIP', 'DNA2', 'PRIM1', 'PHC1', 'KIF14')] <- 1
EntireCohort_pathvar$CCgene[EntireCohort_pathvar$GENE5 %in% c('PDCD6IP', 'AMMECR1', 'MCPH1', 'NCAPD2', 'NCAPD3', 'NCAPH', 'VRK1', 'ESCO2', 'ANKLE2', 'DYRK1A', 'ATR', 'ATRIP', 'DNA2', 'PRIM1', 'PHC1', 'KIF14')] <- 1
EntireCohort_pathvar$CCgene[EntireCohort_pathvar$GENE6 %in% c('PDCD6IP', 'AMMECR1', 'MCPH1', 'NCAPD2', 'NCAPD3', 'NCAPH', 'VRK1', 'ESCO2', 'ANKLE2', 'DYRK1A', 'ATR', 'ATRIP', 'DNA2', 'PRIM1', 'PHC1', 'KIF14')] <- 1
EntireCohort_pathvar$CCgene[EntireCohort_pathvar$GENE7 %in% c('PDCD6IP', 'AMMECR1', 'MCPH1', 'NCAPD2', 'NCAPD3', 'NCAPH', 'VRK1', 'ESCO2', 'ANKLE2', 'DYRK1A', 'ATR', 'ATRIP', 'DNA2', 'PRIM1', 'PHC1', 'KIF14')] <- 1
EntireCohort_pathvar$CCgene[EntireCohort_pathvar$GENE8 %in% c('PDCD6IP', 'AMMECR1', 'MCPH1', 'NCAPD2', 'NCAPD3', 'NCAPH', 'VRK1', 'ESCO2', 'ANKLE2', 'DYRK1A', 'ATR', 'ATRIP', 'DNA2', 'PRIM1', 'PHC1', 'KIF14')] <- 1
EntireCohort_pathvar$CCgene[EntireCohort_pathvar$GENE9 %in% c('PDCD6IP', 'AMMECR1', 'MCPH1', 'NCAPD2', 'NCAPD3', 'NCAPH', 'VRK1', 'ESCO2', 'ANKLE2', 'DYRK1A', 'ATR', 'ATRIP', 'DNA2', 'PRIM1', 'PHC1', 'KIF14')] <- 1


EntireCohort_pathvar$CCgene[is.na(EntireCohort_pathvar$CCgene)] <- 0
EntireCohort_pathvar$CCgene <- as.numeric(EntireCohort_pathvar$CCgene)
unique(EntireCohort_pathvar$CCgene)
# 0 1

table(EntireCohort_pathvar$CCgene)
# 0      1 
# 447493  22306 




# now to make a column for each of the CCgenes
### 'PDCD6IP'
EntireCohort_pathvar$PDCD6IP <- NA
EntireCohort_pathvar$PDCD6IP[EntireCohort_pathvar$GENE1 %in% c("PDCD6IP")] <- 1
EntireCohort_pathvar$PDCD6IP[EntireCohort_pathvar$GENE2 %in% c("PDCD6IP")] <- 1
EntireCohort_pathvar$PDCD6IP[EntireCohort_pathvar$GENE3 %in% c("PDCD6IP")] <- 1
EntireCohort_pathvar$PDCD6IP[EntireCohort_pathvar$GENE4 %in% c("PDCD6IP")] <- 1
EntireCohort_pathvar$PDCD6IP[EntireCohort_pathvar$GENE5 %in% c("PDCD6IP")] <- 1
EntireCohort_pathvar$PDCD6IP[EntireCohort_pathvar$GENE6 %in% c("PDCD6IP")] <- 1
EntireCohort_pathvar$PDCD6IP[EntireCohort_pathvar$GENE7 %in% c("PDCD6IP")] <- 1
EntireCohort_pathvar$PDCD6IP[EntireCohort_pathvar$GENE8 %in% c("PDCD6IP")] <- 1
EntireCohort_pathvar$PDCD6IP[EntireCohort_pathvar$GENE9 %in% c("PDCD6IP")] <- 1

EntireCohort_pathvar$PDCD6IP[is.na(EntireCohort_pathvar$PDCD6IP)] <- 0
EntireCohort_pathvar$PDCD6IP <- as.numeric(EntireCohort_pathvar$PDCD6IP)
unique(EntireCohort_pathvar$PDCD6IP)
# 0 1

table(EntireCohort_pathvar$PDCD6IP)
#  0      1 
# 468907    892 





### 'AMMECR1'
EntireCohort_pathvar$AMMECR1 <- NA
EntireCohort_pathvar$AMMECR1[EntireCohort_pathvar$GENE1 %in% c("AMMECR1")] <- 1
EntireCohort_pathvar$AMMECR1[EntireCohort_pathvar$GENE2 %in% c("AMMECR1")] <- 1
EntireCohort_pathvar$AMMECR1[EntireCohort_pathvar$GENE3 %in% c("AMMECR1")] <- 1
EntireCohort_pathvar$AMMECR1[EntireCohort_pathvar$GENE4 %in% c("AMMECR1")] <- 1
EntireCohort_pathvar$AMMECR1[EntireCohort_pathvar$GENE5 %in% c("AMMECR1")] <- 1
EntireCohort_pathvar$AMMECR1[EntireCohort_pathvar$GENE6 %in% c("AMMECR1")] <- 1
EntireCohort_pathvar$AMMECR1[EntireCohort_pathvar$GENE7 %in% c("AMMECR1")] <- 1
EntireCohort_pathvar$AMMECR1[EntireCohort_pathvar$GENE8 %in% c("AMMECR1")] <- 1
EntireCohort_pathvar$AMMECR1[EntireCohort_pathvar$GENE9 %in% c("AMMECR1")] <- 1

EntireCohort_pathvar$AMMECR1[is.na(EntireCohort_pathvar$AMMECR1)] <- 0
EntireCohort_pathvar$AMMECR1 <- as.numeric(EntireCohort_pathvar$AMMECR1)
unique(EntireCohort_pathvar$AMMECR1)
# 0 1

table(EntireCohort_pathvar$AMMECR1)
# 0      1 
# 469654    145




### 'CDKN1C'
# no variants present




### 'MCPH1'
EntireCohort_pathvar$MCPH1 <- NA
EntireCohort_pathvar$MCPH1[EntireCohort_pathvar$GENE1 %in% c("MCPH1")] <- 1
EntireCohort_pathvar$MCPH1[EntireCohort_pathvar$GENE2 %in% c("MCPH1")] <- 1
EntireCohort_pathvar$MCPH1[EntireCohort_pathvar$GENE3 %in% c("MCPH1")] <- 1
EntireCohort_pathvar$MCPH1[EntireCohort_pathvar$GENE4 %in% c("MCPH1")] <- 1
EntireCohort_pathvar$MCPH1[EntireCohort_pathvar$GENE5 %in% c("MCPH1")] <- 1
EntireCohort_pathvar$MCPH1[EntireCohort_pathvar$GENE6 %in% c("MCPH1")] <- 1
EntireCohort_pathvar$MCPH1[EntireCohort_pathvar$GENE7 %in% c("MCPH1")] <- 1
EntireCohort_pathvar$MCPH1[EntireCohort_pathvar$GENE8 %in% c("MCPH1")] <- 1
EntireCohort_pathvar$MCPH1[EntireCohort_pathvar$GENE9 %in% c("MCPH1")] <- 1

EntireCohort_pathvar$MCPH1[is.na(EntireCohort_pathvar$MCPH1)] <- 0
EntireCohort_pathvar$MCPH1 <- as.numeric(EntireCohort_pathvar$MCPH1)
unique(EntireCohort_pathvar$MCPH1)
# 0 1

table(EntireCohort_pathvar$MCPH1)
#   0      1 
# 467938   1861 




### 'NCAPD2'
EntireCohort_pathvar$NCAPD2 <- NA
EntireCohort_pathvar$NCAPD2[EntireCohort_pathvar$GENE1 %in% c("NCAPD2")] <- 1
EntireCohort_pathvar$NCAPD2[EntireCohort_pathvar$GENE2 %in% c("NCAPD2")] <- 1
EntireCohort_pathvar$NCAPD2[EntireCohort_pathvar$GENE3 %in% c("NCAPD2")] <- 1
EntireCohort_pathvar$NCAPD2[EntireCohort_pathvar$GENE4 %in% c("NCAPD2")] <- 1
EntireCohort_pathvar$NCAPD2[EntireCohort_pathvar$GENE5 %in% c("NCAPD2")] <- 1
EntireCohort_pathvar$NCAPD2[EntireCohort_pathvar$GENE6 %in% c("NCAPD2")] <- 1
EntireCohort_pathvar$NCAPD2[EntireCohort_pathvar$GENE7 %in% c("NCAPD2")] <- 1
EntireCohort_pathvar$NCAPD2[EntireCohort_pathvar$GENE8 %in% c("NCAPD2")] <- 1
EntireCohort_pathvar$NCAPD2[EntireCohort_pathvar$GENE9 %in% c("NCAPD2")] <- 1

EntireCohort_pathvar$NCAPD2[is.na(EntireCohort_pathvar$NCAPD2)] <- 0
EntireCohort_pathvar$NCAPD2 <- as.numeric(EntireCohort_pathvar$NCAPD2)
unique(EntireCohort_pathvar$NCAPD2)
# 0 1

table(EntireCohort_pathvar$NCAPD2)
#   0      1 
# 465301   4498




### 'NCAPD3' 
EntireCohort_pathvar$NCAPD3 <- NA
EntireCohort_pathvar$NCAPD3[EntireCohort_pathvar$GENE1 %in% c("NCAPD3")] <- 1
EntireCohort_pathvar$NCAPD3[EntireCohort_pathvar$GENE2 %in% c("NCAPD3")] <- 1
EntireCohort_pathvar$NCAPD3[EntireCohort_pathvar$GENE3 %in% c("NCAPD3")] <- 1
EntireCohort_pathvar$NCAPD3[EntireCohort_pathvar$GENE4 %in% c("NCAPD3")] <- 1
EntireCohort_pathvar$NCAPD3[EntireCohort_pathvar$GENE5 %in% c("NCAPD3")] <- 1
EntireCohort_pathvar$NCAPD3[EntireCohort_pathvar$GENE6 %in% c("NCAPD3")] <- 1
EntireCohort_pathvar$NCAPD3[EntireCohort_pathvar$GENE7 %in% c("NCAPD3")] <- 1
EntireCohort_pathvar$NCAPD3[EntireCohort_pathvar$GENE8 %in% c("NCAPD3")] <- 1
EntireCohort_pathvar$NCAPD3[EntireCohort_pathvar$GENE9 %in% c("NCAPD3")] <- 1

EntireCohort_pathvar$NCAPD3[is.na(EntireCohort_pathvar$NCAPD3)] <- 0
EntireCohort_pathvar$NCAPD3 <- as.numeric(EntireCohort_pathvar$NCAPD3)
unique(EntireCohort_pathvar$NCAPD3)
# 0 1

table(EntireCohort_pathvar$NCAPD3)
#   0      1 
# 468484   1315




### 'NCAPH'
EntireCohort_pathvar$NCAPH <- NA
EntireCohort_pathvar$NCAPH[EntireCohort_pathvar$GENE1 %in% c("NCAPH")] <- 1
EntireCohort_pathvar$NCAPH[EntireCohort_pathvar$GENE2 %in% c("NCAPH")] <- 1
EntireCohort_pathvar$NCAPH[EntireCohort_pathvar$GENE3 %in% c("NCAPH")] <- 1
EntireCohort_pathvar$NCAPH[EntireCohort_pathvar$GENE4 %in% c("NCAPH")] <- 1
EntireCohort_pathvar$NCAPH[EntireCohort_pathvar$GENE5 %in% c("NCAPH")] <- 1
EntireCohort_pathvar$NCAPH[EntireCohort_pathvar$GENE6 %in% c("NCAPH")] <- 1
EntireCohort_pathvar$NCAPH[EntireCohort_pathvar$GENE7 %in% c("NCAPH")] <- 1
EntireCohort_pathvar$NCAPH[EntireCohort_pathvar$GENE8 %in% c("NCAPH")] <- 1
EntireCohort_pathvar$NCAPH[EntireCohort_pathvar$GENE9 %in% c("NCAPH")] <- 1

EntireCohort_pathvar$NCAPH[is.na(EntireCohort_pathvar$NCAPH)] <- 0
EntireCohort_pathvar$NCAPH <- as.numeric(EntireCohort_pathvar$NCAPH)
unique(EntireCohort_pathvar$NCAPH)
# 0 1

table(EntireCohort_pathvar$NCAPH)
#  0      1 
# 468217   1582




### 'VRK1' 
EntireCohort_pathvar$VRK1 <- NA
EntireCohort_pathvar$VRK1[EntireCohort_pathvar$GENE1 %in% c("VRK1")] <- 1
EntireCohort_pathvar$VRK1[EntireCohort_pathvar$GENE2 %in% c("VRK1")] <- 1
EntireCohort_pathvar$VRK1[EntireCohort_pathvar$GENE3 %in% c("VRK1")] <- 1
EntireCohort_pathvar$VRK1[EntireCohort_pathvar$GENE4 %in% c("VRK1")] <- 1
EntireCohort_pathvar$VRK1[EntireCohort_pathvar$GENE5 %in% c("VRK1")] <- 1
EntireCohort_pathvar$VRK1[EntireCohort_pathvar$GENE6 %in% c("VRK1")] <- 1
EntireCohort_pathvar$VRK1[EntireCohort_pathvar$GENE7 %in% c("VRK1")] <- 1
EntireCohort_pathvar$VRK1[EntireCohort_pathvar$GENE8 %in% c("VRK1")] <- 1
EntireCohort_pathvar$VRK1[EntireCohort_pathvar$GENE9 %in% c("VRK1")] <- 1

EntireCohort_pathvar$VRK1[is.na(EntireCohort_pathvar$VRK1)] <- 0
EntireCohort_pathvar$VRK1 <- as.numeric(EntireCohort_pathvar$VRK1)
unique(EntireCohort_pathvar$VRK1)
# 0 1

table(EntireCohort_pathvar$VRK1)
# 0      1 
# 469368    431




### 'ESCO2'
EntireCohort_pathvar$ESCO2 <- NA
EntireCohort_pathvar$ESCO2[EntireCohort_pathvar$GENE1 %in% c("ESCO2")] <- 1
EntireCohort_pathvar$ESCO2[EntireCohort_pathvar$GENE2 %in% c("ESCO2")] <- 1
EntireCohort_pathvar$ESCO2[EntireCohort_pathvar$GENE3 %in% c("ESCO2")] <- 1
EntireCohort_pathvar$ESCO2[EntireCohort_pathvar$GENE4 %in% c("ESCO2")] <- 1
EntireCohort_pathvar$ESCO2[EntireCohort_pathvar$GENE5 %in% c("ESCO2")] <- 1
EntireCohort_pathvar$ESCO2[EntireCohort_pathvar$GENE6 %in% c("ESCO2")] <- 1
EntireCohort_pathvar$ESCO2[EntireCohort_pathvar$GENE7 %in% c("ESCO2")] <- 1
EntireCohort_pathvar$ESCO2[EntireCohort_pathvar$GENE8 %in% c("ESCO2")] <- 1
EntireCohort_pathvar$ESCO2[EntireCohort_pathvar$GENE9 %in% c("ESCO2")] <- 1

EntireCohort_pathvar$ESCO2[is.na(EntireCohort_pathvar$ESCO2)] <- 0
EntireCohort_pathvar$ESCO2 <- as.numeric(EntireCohort_pathvar$ESCO2)
unique(EntireCohort_pathvar$ESCO2)
# 0 1

table(EntireCohort_pathvar$ESCO2)
#  0      1 
# 469448    351





### 'ANKLE2'  
EntireCohort_pathvar$ANKLE2 <- NA
EntireCohort_pathvar$ANKLE2[EntireCohort_pathvar$GENE1 %in% c("ANKLE2")] <- 1
EntireCohort_pathvar$ANKLE2[EntireCohort_pathvar$GENE2 %in% c("ANKLE2")] <- 1
EntireCohort_pathvar$ANKLE2[EntireCohort_pathvar$GENE3 %in% c("ANKLE2")] <- 1
EntireCohort_pathvar$ANKLE2[EntireCohort_pathvar$GENE4 %in% c("ANKLE2")] <- 1
EntireCohort_pathvar$ANKLE2[EntireCohort_pathvar$GENE5 %in% c("ANKLE2")] <- 1
EntireCohort_pathvar$ANKLE2[EntireCohort_pathvar$GENE6 %in% c("ANKLE2")] <- 1
EntireCohort_pathvar$ANKLE2[EntireCohort_pathvar$GENE7 %in% c("ANKLE2")] <- 1
EntireCohort_pathvar$ANKLE2[EntireCohort_pathvar$GENE8 %in% c("ANKLE2")] <- 1
EntireCohort_pathvar$ANKLE2[EntireCohort_pathvar$GENE9 %in% c("ANKLE2")] <- 1

EntireCohort_pathvar$ANKLE2[is.na(EntireCohort_pathvar$ANKLE2)] <- 0
EntireCohort_pathvar$ANKLE2 <- as.numeric(EntireCohort_pathvar$ANKLE2)
unique(EntireCohort_pathvar$ANKLE2)
# 0 1

table(EntireCohort_pathvar$ANKLE2)
# 0      1 
# 468601   1198




### 'DYRK1A'
EntireCohort_pathvar$DYRK1A <- NA
EntireCohort_pathvar$DYRK1A[EntireCohort_pathvar$GENE1 %in% c("DYRK1A")] <- 1
EntireCohort_pathvar$DYRK1A[EntireCohort_pathvar$GENE2 %in% c("DYRK1A")] <- 1
EntireCohort_pathvar$DYRK1A[EntireCohort_pathvar$GENE3 %in% c("DYRK1A")] <- 1
EntireCohort_pathvar$DYRK1A[EntireCohort_pathvar$GENE4 %in% c("DYRK1A")] <- 1
EntireCohort_pathvar$DYRK1A[EntireCohort_pathvar$GENE5 %in% c("DYRK1A")] <- 1
EntireCohort_pathvar$DYRK1A[EntireCohort_pathvar$GENE6 %in% c("DYRK1A")] <- 1
EntireCohort_pathvar$DYRK1A[EntireCohort_pathvar$GENE7 %in% c("DYRK1A")] <- 1
EntireCohort_pathvar$DYRK1A[EntireCohort_pathvar$GENE8 %in% c("DYRK1A")] <- 1
EntireCohort_pathvar$DYRK1A[EntireCohort_pathvar$GENE9 %in% c("DYRK1A")] <- 1

EntireCohort_pathvar$DYRK1A[is.na(EntireCohort_pathvar$DYRK1A)] <- 0
EntireCohort_pathvar$DYRK1A <- as.numeric(EntireCohort_pathvar$DYRK1A)
unique(EntireCohort_pathvar$DYRK1A)
# 0 1

table(EntireCohort_pathvar$DYRK1A)
#  0      1 
# 468060   1739




### 'ATR'
EntireCohort_pathvar$ATR <- NA
EntireCohort_pathvar$ATR[EntireCohort_pathvar$GENE1 %in% c("ATR")] <- 1
EntireCohort_pathvar$ATR[EntireCohort_pathvar$GENE2 %in% c("ATR")] <- 1
EntireCohort_pathvar$ATR[EntireCohort_pathvar$GENE3 %in% c("ATR")] <- 1
EntireCohort_pathvar$ATR[EntireCohort_pathvar$GENE4 %in% c("ATR")] <- 1
EntireCohort_pathvar$ATR[EntireCohort_pathvar$GENE5 %in% c("ATR")] <- 1
EntireCohort_pathvar$ATR[EntireCohort_pathvar$GENE6 %in% c("ATR")] <- 1
EntireCohort_pathvar$ATR[EntireCohort_pathvar$GENE7 %in% c("ATR")] <- 1
EntireCohort_pathvar$ATR[EntireCohort_pathvar$GENE8 %in% c("ATR")] <- 1
EntireCohort_pathvar$ATR[EntireCohort_pathvar$GENE9 %in% c("ATR")] <- 1

EntireCohort_pathvar$ATR[is.na(EntireCohort_pathvar$ATR)] <- 0
EntireCohort_pathvar$ATR <- as.numeric(EntireCohort_pathvar$ATR)
unique(EntireCohort_pathvar$ATR)
# 0 1

table(EntireCohort_pathvar$ATR)
#  0      1 
# 466447   3352 




### 'ATRIP' 
EntireCohort_pathvar$ATRIP <- NA
EntireCohort_pathvar$ATRIP[EntireCohort_pathvar$GENE1 %in% c("ATRIP")] <- 1
EntireCohort_pathvar$ATRIP[EntireCohort_pathvar$GENE2 %in% c("ATRIP")] <- 1
EntireCohort_pathvar$ATRIP[EntireCohort_pathvar$GENE3 %in% c("ATRIP")] <- 1
EntireCohort_pathvar$ATRIP[EntireCohort_pathvar$GENE4 %in% c("ATRIP")] <- 1
EntireCohort_pathvar$ATRIP[EntireCohort_pathvar$GENE5 %in% c("ATRIP")] <- 1
EntireCohort_pathvar$ATRIP[EntireCohort_pathvar$GENE6 %in% c("ATRIP")] <- 1
EntireCohort_pathvar$ATRIP[EntireCohort_pathvar$GENE7 %in% c("ATRIP")] <- 1
EntireCohort_pathvar$ATRIP[EntireCohort_pathvar$GENE8 %in% c("ATRIP")] <- 1
EntireCohort_pathvar$ATRIP[EntireCohort_pathvar$GENE9 %in% c("ATRIP")] <- 1

EntireCohort_pathvar$ATRIP[is.na(EntireCohort_pathvar$ATRIP)] <- 0
EntireCohort_pathvar$ATRIP <- as.numeric(EntireCohort_pathvar$ATRIP)
unique(EntireCohort_pathvar$ATRIP)
# 0 1

table(EntireCohort_pathvar$ATRIP)
#   0      1 
# 469135    664




### 'DNA2' 
EntireCohort_pathvar$DNA2 <- NA
EntireCohort_pathvar$DNA2[EntireCohort_pathvar$GENE1 %in% c("DNA2")] <- 1
EntireCohort_pathvar$DNA2[EntireCohort_pathvar$GENE2 %in% c("DNA2")] <- 1
EntireCohort_pathvar$DNA2[EntireCohort_pathvar$GENE3 %in% c("DNA2")] <- 1
EntireCohort_pathvar$DNA2[EntireCohort_pathvar$GENE4 %in% c("DNA2")] <- 1
EntireCohort_pathvar$DNA2[EntireCohort_pathvar$GENE5 %in% c("DNA2")] <- 1
EntireCohort_pathvar$DNA2[EntireCohort_pathvar$GENE6 %in% c("DNA2")] <- 1
EntireCohort_pathvar$DNA2[EntireCohort_pathvar$GENE7 %in% c("DNA2")] <- 1
EntireCohort_pathvar$DNA2[EntireCohort_pathvar$GENE8 %in% c("DNA2")] <- 1
EntireCohort_pathvar$DNA2[EntireCohort_pathvar$GENE9 %in% c("DNA2")] <- 1

EntireCohort_pathvar$DNA2[is.na(EntireCohort_pathvar$DNA2)] <- 0
EntireCohort_pathvar$DNA2 <- as.numeric(EntireCohort_pathvar$DNA2)
unique(EntireCohort_pathvar$DNA2)
# 0 1

table(EntireCohort_pathvar$DNA2)
# 0      1 
# 467624   2175




### 'FAM111A'
# no variants present





### 'PRIM1' 
EntireCohort_pathvar$PRIM1 <- NA
EntireCohort_pathvar$PRIM1[EntireCohort_pathvar$GENE1 %in% c("PRIM1")] <- 1
EntireCohort_pathvar$PRIM1[EntireCohort_pathvar$GENE2 %in% c("PRIM1")] <- 1
EntireCohort_pathvar$PRIM1[EntireCohort_pathvar$GENE3 %in% c("PRIM1")] <- 1
EntireCohort_pathvar$PRIM1[EntireCohort_pathvar$GENE4 %in% c("PRIM1")] <- 1
EntireCohort_pathvar$PRIM1[EntireCohort_pathvar$GENE5 %in% c("PRIM1")] <- 1
EntireCohort_pathvar$PRIM1[EntireCohort_pathvar$GENE6 %in% c("PRIM1")] <- 1
EntireCohort_pathvar$PRIM1[EntireCohort_pathvar$GENE7 %in% c("PRIM1")] <- 1
EntireCohort_pathvar$PRIM1[EntireCohort_pathvar$GENE8 %in% c("PRIM1")] <- 1
EntireCohort_pathvar$PRIM1[EntireCohort_pathvar$GENE9 %in% c("PRIM1")] <- 1

EntireCohort_pathvar$PRIM1[is.na(EntireCohort_pathvar$PRIM1)] <- 0
EntireCohort_pathvar$PRIM1 <- as.numeric(EntireCohort_pathvar$PRIM1)
unique(EntireCohort_pathvar$PRIM1)
# 0 1

table(EntireCohort_pathvar$PRIM1)
#  0      1 
# 469308    491





### 'PHC1'
EntireCohort_pathvar$PHC1 <- NA
EntireCohort_pathvar$PHC1[EntireCohort_pathvar$GENE1 %in% c("PHC1")] <- 1
EntireCohort_pathvar$PHC1[EntireCohort_pathvar$GENE2 %in% c("PHC1")] <- 1
EntireCohort_pathvar$PHC1[EntireCohort_pathvar$GENE3 %in% c("PHC1")] <- 1
EntireCohort_pathvar$PHC1[EntireCohort_pathvar$GENE4 %in% c("PHC1")] <- 1
EntireCohort_pathvar$PHC1[EntireCohort_pathvar$GENE5 %in% c("PHC1")] <- 1
EntireCohort_pathvar$PHC1[EntireCohort_pathvar$GENE6 %in% c("PHC1")] <- 1
EntireCohort_pathvar$PHC1[EntireCohort_pathvar$GENE7 %in% c("PHC1")] <- 1
EntireCohort_pathvar$PHC1[EntireCohort_pathvar$GENE8 %in% c("PHC1")] <- 1
EntireCohort_pathvar$PHC1[EntireCohort_pathvar$GENE9 %in% c("PHC1")] <- 1

EntireCohort_pathvar$PHC1[is.na(EntireCohort_pathvar$PHC1)] <- 0
EntireCohort_pathvar$PHC1 <- as.numeric(EntireCohort_pathvar$PHC1)
unique(EntireCohort_pathvar$PHC1)
# 0 1

table(EntireCohort_pathvar$PHC1)
#   0      1 
# 469070    729




### 'KIF14'
EntireCohort_pathvar$KIF14 <- NA
EntireCohort_pathvar$KIF14[EntireCohort_pathvar$GENE1 %in% c("KIF14")] <- 1
EntireCohort_pathvar$KIF14[EntireCohort_pathvar$GENE2 %in% c("KIF14")] <- 1
EntireCohort_pathvar$KIF14[EntireCohort_pathvar$GENE3 %in% c("KIF14")] <- 1
EntireCohort_pathvar$KIF14[EntireCohort_pathvar$GENE4 %in% c("KIF14")] <- 1
EntireCohort_pathvar$KIF14[EntireCohort_pathvar$GENE5 %in% c("KIF14")] <- 1
EntireCohort_pathvar$KIF14[EntireCohort_pathvar$GENE6 %in% c("KIF14")] <- 1
EntireCohort_pathvar$KIF14[EntireCohort_pathvar$GENE7 %in% c("KIF14")] <- 1
EntireCohort_pathvar$KIF14[EntireCohort_pathvar$GENE8 %in% c("KIF14")] <- 1
EntireCohort_pathvar$KIF14[EntireCohort_pathvar$GENE9 %in% c("KIF14")] <- 1

EntireCohort_pathvar$KIF14[is.na(EntireCohort_pathvar$KIF14)] <- 0
EntireCohort_pathvar$KIF14 <- as.numeric(EntireCohort_pathvar$KIF14)
unique(EntireCohort_pathvar$KIF14)
# 0 1

table(EntireCohort_pathvar$KIF14)
#  0      1 
# 468442   1357 








# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)

















# CCgene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_CCgene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_CCgene_, "pathvar_CCgene_.csv", row.names=FALSE)
# system('dx upload pathvar_CCgene_.csv --path Emily-folder/results_2/pathvar_CCgene_.csv')


# CCgene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ CCgene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_CCgene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_incidence_1, "pathvar_CCgene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_incidence_1.csv --path Emily-folder/results_2/pathvar_CCgene_incidence_1.csv')


# GOIs
LR_2 <- glm(formula = cancer_ROOC_dummy ~ PDCD6IP+AMMECR1+MCPH1+NCAPD2+NCAPD3+NCAPH+VRK1+ESCO2+ANKLE2+DYRK1A+ATR+ATRIP+DNA2+PRIM1+PHC1+KIF14, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_CCgene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_incidence_2, "pathvar_CCgene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_incidence_2.csv --path Emily-folder/results_2/pathvar_CCgene_incidence_2.csv')







### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "CCgene", 'PDCD6IP', 'AMMECR1', 'MCPH1', 
                                     'NCAPD2', 'NCAPD3', 'NCAPH', 'VRK1', 'ESCO2', 'ANKLE2', 
                                     'DYRK1A', 'ATR', 'ATRIP', 'DNA2', 'PRIM1', 
                                     'PHC1', 'KIF14')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
#465881     23

with(data_setA, table(KIF14))




# CCgene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ CCgene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_CCgene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_incidence_3, "pathvar_CCgene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_incidence_3.csv --path Emily-folder/results_2/pathvar_CCgene_incidence_3.csv')


# age + ethnicity + smoking + sex + CCgene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_CCgene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_incidence_4, "pathvar_CCgene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_incidence_4.csv --path Emily-folder/results_2/pathvar_CCgene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ PDCD6IP+AMMECR1+MCPH1+NCAPD2+NCAPD3+NCAPH+VRK1+ESCO2+ANKLE2+DYRK1A+ATR+ATRIP+DNA2+PRIM1+PHC1+KIF14, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_CCgene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_incidence_5, "pathvar_CCgene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_incidence_5.csv --path Emily-folder/results_2/pathvar_CCgene_incidence_5.csv')


# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+PDCD6IP+AMMECR1+MCPH1+NCAPD2+NCAPD3+NCAPH+VRK1+ESCO2+ANKLE2+DYRK1A+ATR+ATRIP+DNA2+PRIM1+PHC1+KIF14, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_CCgene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_incidence_6, "pathvar_CCgene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_incidence_6.csv --path Emily-folder/results_2/pathvar_CCgene_incidence_6.csv')




### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "CCgene", 'PDCD6IP', 'AMMECR1', 'MCPH1', 
                                     'NCAPD2', 'NCAPD3', 'NCAPH', 'VRK1', 'ESCO2', 'ANKLE2', 
                                     'DYRK1A', 'ATR', 'ATRIP', 'DNA2', 'PRIM1', 
                                     'PHC1', 'KIF14')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     28

with(data_setB, table(KIF14))





# CCgene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ CCgene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_CCgene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_incidence_7, "pathvar_CCgene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_incidence_7.csv --path Emily-folder/results_2/pathvar_CCgene_incidence_7.csv')


# all variables + CCgene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+CCgene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_CCgene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_incidence_8, "pathvar_CCgene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_incidence_8.csv --path Emily-folder/results_2/pathvar_CCgene_incidence_8.csv')


# GOIs alone
LR_9 <- glm(formula = cancer_ROOC_dummy ~ PDCD6IP+AMMECR1+MCPH1+NCAPD2+NCAPD3+NCAPH+VRK1+ESCO2+ANKLE2+DYRK1A+ATR+ATRIP+DNA2+PRIM1+PHC1+KIF14, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_CCgene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_incidence_9, "pathvar_CCgene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_incidence_9.csv --path Emily-folder/results_2/pathvar_CCgene_incidence_9.csv')


# all variables + GOIs
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+PDCD6IP+AMMECR1+MCPH1+NCAPD2+NCAPD3+NCAPH+VRK1+ESCO2+ANKLE2+DYRK1A+ATR+ATRIP+DNA2+PRIM1+PHC1+KIF14, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_CCgene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_incidence_10, "pathvar_CCgene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_incidence_10.csv --path Emily-folder/results_2/pathvar_CCgene_incidence_10.csv')










### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(CCgene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(CCgene, cancer_ROOC_dummy)))
# p-value = 0.1906

with(EntireCohort_pathvar, table(CCgene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(CCgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5983833 0.6622905
# sample estimates:
#   odds ratio 
# 0.6293926 

with(EntireCohort_pathvar, table(CCgene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(CCgene, Ever_smoked_dummy)))
# p-value = 0.01554
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9404151 0.9936089
# sample estimates:
#   odds ratio 
# 0.9666076

with(EntireCohort_pathvar, table(CCgene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(CCgene, Sex_Male)))
# p-value = 0.3017

with(EntireCohort_pathvar, table(CCgene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(CCgene, University_education)))
# p-value = 0.04091
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.001150 1.060475
# sample estimates:
#   odds ratio 
# 1.030461

with(EntireCohort_pathvar, table(CCgene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(CCgene, Annual_income))	)
# p-value = 0.7088

with(EntireCohort_pathvar, table(CCgene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(CCgene, Smoker_previous)))
# p-value = 0.04239
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9435799 0.9990237
# sample estimates:
#   odds ratio 
# 0.9709472 

with(EntireCohort_pathvar, table(CCgene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(CCgene, Smoker_current))	)
# p-value = 0.8668



### data_setA
with(data_setA, table(CCgene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(CCgene, cancer_ROOC_dummy)))
# p-value = 0.1534

with(data_setA, table(CCgene, Ethnicity_White))			
fisher.test(with(data_setA, table(CCgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5968255 0.6607669
# sample estimates:
#   odds ratio 
# 0.6278522

with(data_setA, table(CCgene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(CCgene, Ever_smoked_dummy)))
# p-value = 0.01758
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9409775 0.9942930
# sample estimates:
#   odds ratio 
# 0.9672287

with(data_setA, table(CCgene, Sex_Male))			
fisher.test(with(data_setA, table(CCgene, Sex_Male)))
# p-value = 0.2569



### data_setB
with(data_setB, table(CCgene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(CCgene, cancer_ROOC_dummy)))
# p-value = 0.23

with(data_setB, table(CCgene, Ethnicity_White))			
fisher.test(with(data_setB, table(CCgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5769521 0.6511081
# sample estimates:
#   odds ratio 
# 0.6126878

with(data_setB, table(CCgene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(CCgene, Ever_smoked_dummy)))
# p-value = 0.07955

with(data_setB, table(CCgene, Sex_Male))			
fisher.test(with(data_setB, table(CCgene, Sex_Male)))
# p-value = 0.2495

with(data_setB, table(CCgene, University_education))			
fisher.test(with(data_setB, table(CCgene, University_education)))
# p-value = 0.06878

with(data_setB, table(CCgene, Annual_income))			
fisher.test(with(data_setB, table(CCgene, Annual_income)))
# p-value = 0.7823






### can plot the continuous variables: 

# make CCgene into factor
df_for_plots <- EntireCohort_pathvar
df_for_plots$CCgene <- as.factor(df_for_plots$CCgene)

# try a violin plot

# Age_at_recruitment  
df_for_plots %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ CCgene, data = df_for_plots)
# W = 5005409698, p-value = 0.4623

# BMI  
df_for_plots %>% 
  ggplot(mapping = aes(x = CCgene, y = BMI)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = CCgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ CCgene, data = df_for_plots)
# W = 4939468240, p-value = 0.6967

# Standing_height 
df_for_plots %>% 
  ggplot(mapping = aes(x = CCgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = CCgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ CCgene, data = df_for_plots)
# W = 4.997e+09, p-value = 0.05183

# make this a bar plot
# Weekly_exercise
df_for_plots %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(CCgene), scales = "free_y") +
  theme_light()
# 



### data_setA

# make CCgene into factor
data_setA$CCgene <- as.factor(data_setA$CCgene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ CCgene, data = data_setA)
# W = 4923803254, p-value = 0.4626



### data_setB

# make CCgene into factor
data_setB$CCgene <- as.factor(data_setB$CCgene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ CCgene, data = data_setB)
# W = 3287225058, p-value = 0.7115

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = CCgene, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = CCgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ CCgene, data = data_setB)
# W = 3275328196, p-value = 0.6495

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = CCgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = CCgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ CCgene, data = data_setB)
# W = 3309369496, p-value = 0.05691


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(CCgene), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# CCgene Type ------------------------------------------------------------

# make sure the execute the code in the 'CCgenes' section

## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_CCgene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_type_1, "pathvar_CCgene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_type_1.csv --path Emily-folder/results_2/pathvar_CCgene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_CCgene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_type_2, "pathvar_CCgene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_type_2.csv --path Emily-folder/results_2/pathvar_CCgene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_CCgene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_type_3, "pathvar_CCgene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_type_3.csv --path Emily-folder/results_2/pathvar_CCgene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_CCgene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_type_4, "pathvar_CCgene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_type_4.csv --path Emily-folder/results_2/pathvar_CCgene_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_CCgene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_type_5, "pathvar_CCgene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_type_5.csv --path Emily-folder/results_2/pathvar_CCgene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_CCgene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_type_6, "pathvar_CCgene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_type_6.csv --path Emily-folder/results_2/pathvar_CCgene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_CCgene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_type_7, "pathvar_CCgene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_type_7.csv --path Emily-folder/results_2/pathvar_CCgene_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_CCgene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_type_8, "pathvar_CCgene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_type_8.csv --path Emily-folder/results_2/pathvar_CCgene_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_CCgene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_type_9, "pathvar_CCgene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_type_9.csv --path Emily-folder/results_2/pathvar_CCgene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_CCgene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_type_10, "pathvar_CCgene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_type_10.csv --path Emily-folder/results_2/pathvar_CCgene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_CCgene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_type_11, "pathvar_CCgene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_type_11.csv --path Emily-folder/results_2/pathvar_CCgene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_CCgene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_type_12, "pathvar_CCgene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_type_12.csv --path Emily-folder/results_2/pathvar_CCgene_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_CCgene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_type_13, "pathvar_CCgene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_type_13.csv --path Emily-folder/results_2/pathvar_CCgene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_CCgene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_type_14, "pathvar_CCgene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_type_14.csv --path Emily-folder/results_2/pathvar_CCgene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_CCgene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CCgene_type_15, "pathvar_CCgene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_type_15.csv --path Emily-folder/results_2/pathvar_CCgene_type_15.csv')










# CCgene Age -------------------------------------------------------------

# make sure the execute the code in the 'CCgenes' section

### linear regression
# "CCgene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ CCgene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_CCgene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_CCgene_age_1, "pathvar_CCgene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_age_1.csv --path Emily-folder/results_2/pathvar_CCgene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ PDCD6IP+AMMECR1+MCPH1+NCAPD2+NCAPD3+NCAPH+VRK1+ESCO2+ANKLE2+DYRK1A+ATR+ATRIP+DNA2+PRIM1+PHC1+KIF14, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_CCgene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_CCgene_age_2, "pathvar_CCgene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_age_2.csv --path Emily-folder/results_2/pathvar_CCgene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "CCgene", 'PDCD6IP', 'AMMECR1', 'MCPH1', 'NCAPD2', 
                                     'NCAPD3', 'NCAPH', 'VRK1', 'ESCO2', 'ANKLE2', 'DYRK1A', 
                                     'ATR', 'ATRIP', 'DNA2', 'PRIM1', 'PHC1', 
                                     'KIF14')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     23

with(data_setD, table(KIF14))
# AMMECR1 = 30
# VRK1 = 90
# ESCO2 = 77



data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "CCgene", 'PDCD6IP', 'AMMECR1', 
                                     'MCPH1', 'NCAPD2', 'NCAPD3', 'NCAPH', 'VRK1', 'ESCO2', 'ANKLE2', 
                                     'DYRK1A', 'ATR', 'ATRIP', 'DNA2', 'PRIM1', 'PHC1', 
                                     'KIF14')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    28

with(data_setE, table(KIF14))
# AMMECR1 = 27
# VRK1 = 69
# ESCO2 = 61
# PRIM1 = 87



### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ CCgene, data = data_setD)
# summary(LM_15)
pathvar_CCgene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_CCgene_age_3, "pathvar_CCgene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_age_3.csv --path Emily-folder/results_2/pathvar_CCgene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CCgene, data = data_setD)
# summary(LM_16)
pathvar_CCgene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_CCgene_age_4, "pathvar_CCgene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_age_4.csv --path Emily-folder/results_2/pathvar_CCgene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ PDCD6IP+AMMECR1+MCPH1+NCAPD2+NCAPD3+NCAPH+VRK1+ESCO2+ANKLE2+DYRK1A+ATR+ATRIP+DNA2+PRIM1+PHC1+KIF14, data = data_setD)
# summary(LM_17)
pathvar_CCgene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_CCgene_age_5, "pathvar_CCgene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_age_5.csv --path Emily-folder/results_2/pathvar_CCgene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+PDCD6IP+AMMECR1+MCPH1+NCAPD2+NCAPD3+NCAPH+VRK1+ESCO2+ANKLE2+DYRK1A+ATR+ATRIP+DNA2+PRIM1+PHC1+KIF14, data = data_setD)
# summary(LM_18)
pathvar_CCgene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_CCgene_age_6, "pathvar_CCgene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_age_6.csv --path Emily-folder/results_2/pathvar_CCgene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ CCgene, data = data_setE)
# summary(LM_15)
pathvar_CCgene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_CCgene_age_7, "pathvar_CCgene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_age_7.csv --path Emily-folder/results_2/pathvar_CCgene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+CCgene, data = data_setE)
# summary(LM_16)
pathvar_CCgene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_CCgene_age_8, "pathvar_CCgene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_age_8.csv --path Emily-folder/results_2/pathvar_CCgene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ PDCD6IP+AMMECR1+MCPH1+NCAPD2+NCAPD3+NCAPH+VRK1+ESCO2+ANKLE2+DYRK1A+ATR+ATRIP+DNA2+PRIM1+PHC1+KIF14, data = data_setE)
# summary(LM_17)
pathvar_CCgene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_CCgene_age_9, "pathvar_CCgene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_age_9.csv --path Emily-folder/results_2/pathvar_CCgene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+PDCD6IP+AMMECR1+MCPH1+NCAPD2+NCAPD3+NCAPH+VRK1+ESCO2+ANKLE2+DYRK1A+ATR+ATRIP+DNA2+PRIM1+PHC1+KIF14, data = data_setE)
# summary(LM_18)
pathvar_CCgene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_CCgene_age_10, "pathvar_CCgene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_CCgene_age_10.csv --path Emily-folder/results_2/pathvar_CCgene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(df_cancer_diagnosis, table(CCgene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(CCgene, Ethnicity_White)))
# p-value = 3.916e-09
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5494978 0.7355796
# sample estimates:
#   odds ratio 
# 0.6344078

with(df_cancer_diagnosis, table(CCgene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(CCgene, Ever_smoked_dummy)))
# p-value = 0.5639

with(df_cancer_diagnosis, table(CCgene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(CCgene, Sex_Male)))
# p-value = 0.9886

with(df_cancer_diagnosis, table(CCgene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(CCgene, University_education)))
# p-value = 0.4822

with(df_cancer_diagnosis, table(CCgene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(CCgene, Annual_income)))
# p-value = 0.5033


### data_setD
with(data_setD, table(CCgene, Ethnicity_White))			
fisher.test(with(data_setD, table(CCgene, Ethnicity_White)))
# p-value = 1.907e-09
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5438609 0.7281435
# sample estimates:
#   odds ratio 
# 0.6279445

with(data_setD, table(CCgene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(CCgene, Ever_smoked_dummy)))
# p-value = 0.5433

with(data_setD, table(CCgene, Sex_Male))			
fisher.test(with(data_setD, table(CCgene, Sex_Male)))
# p-value = 0.8973


### data_setE
with(data_setE, table(CCgene, Ethnicity_White))			
fisher.test(with(data_setE, table(CCgene, Ethnicity_White)))
# p-value = 2.087e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5184411 0.7378525
# sample estimates:
#   odds ratio 
# 0.6165207

with(data_setE, table(CCgene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(CCgene, Ever_smoked_dummy)))
# p-value = 0.8038

with(data_setE, table(CCgene, Sex_Male))			
fisher.test(with(data_setE, table(CCgene, Sex_Male)))
# p-value = 0.9618

with(data_setE, table(CCgene, University_education))			
fisher.test(with(data_setE, table(CCgene, University_education)))
# p-value = 0.3435

with(data_setE, table(CCgene, Annual_income))			
fisher.test(with(data_setE, table(CCgene, Annual_income)))
# p-value = 0.4966



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 

# make CCgene into factor
df_cancer_diagnosis$CCgene <- as.factor(df_cancer_diagnosis$CCgene)

# try a violin plot



# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ CCgene, data = df_cancer_diagnosis)
# 



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ CCgene, data = df_cancer_diagnosis)
# W = 274201204, p-value = 0.8909



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CCgene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CCgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ CCgene, data = df_cancer_diagnosis)
# W = 273322386, p-value = 0.5492



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CCgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CCgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ CCgene, data = df_cancer_diagnosis)
# W = 273957423, p-value = 0.5338



# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(CCgene), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$CCgene <- as.factor(data_setD$CCgene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ CCgene, data = data_setD)
# 



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ CCgene, data = data_setD)
# W = 269249414, p-value = 0.9022




### data_setE

data_setE$CCgene <- as.factor(data_setE$CCgene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ CCgene, data = data_setE)
# 



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = CCgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ CCgene, data = data_setE)
# W = 174460818, p-value = 0.8323



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = CCgene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = CCgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ CCgene, data = data_setE)
# W = 174571506, p-value = 0.7791



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = CCgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = CCgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ CCgene, data = data_setE)
# W = 175498938, p-value = 0.3899


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(CCgene), scales = "free_y") +
  theme_light()
# 










# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_CCgenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_CCgenes.R')
