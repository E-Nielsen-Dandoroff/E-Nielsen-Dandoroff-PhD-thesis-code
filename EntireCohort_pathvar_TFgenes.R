# 09/07/2024

# EntireCohort_pathvar_TFgenes.R

# this script is for the Transcription factor associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)


# TFgene
# ERCC2
# CREB3L1
# NKX3_2
# ARX
# HMGB3
# RAD21
# SMC1A
# SMC3
# NIPBL
# SMARCA2
# SMARCE1
# SRCAP
# DNMT3A
# HDAC8
# PQBP1
# FOXG1
# HMGA2
# LARP7
# SHOX
# SOX11
# TAF13
# TBX15
# TBX3
# SP7
# TRPS1
# ZNF335
# GSC









# TFgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')


# go thr all the genes to check that they're present
grep("GSC", EntireCohort_pathvar$GENE)
# all seem to be present, except for:
# ARX = 0
# SMARCA2 = 0
# SMARCE1 = 0
# SRCAP = 0
# DNMT3A = 0




# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$TFgene <- NA
EntireCohort_pathvar$TFgene[EntireCohort_pathvar$GENE1 %in% c('ERCC2', 'CREB3L1', 'NKX3_2', 'HMGB3', 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 'HDAC8', 'PQBP1', 'FOXG1', 'HMGA2', 'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 'SP7', 'TRPS1', 'ZNF335', 'GSC')] <- 1
EntireCohort_pathvar$TFgene[EntireCohort_pathvar$GENE2 %in% c('ERCC2', 'CREB3L1', 'NKX3_2', 'HMGB3', 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 'HDAC8', 'PQBP1', 'FOXG1', 'HMGA2', 'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 'SP7', 'TRPS1', 'ZNF335', 'GSC')] <- 1
EntireCohort_pathvar$TFgene[EntireCohort_pathvar$GENE3 %in% c('ERCC2', 'CREB3L1', 'NKX3_2', 'HMGB3', 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 'HDAC8', 'PQBP1', 'FOXG1', 'HMGA2', 'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 'SP7', 'TRPS1', 'ZNF335', 'GSC')] <- 1
EntireCohort_pathvar$TFgene[EntireCohort_pathvar$GENE4 %in% c('ERCC2', 'CREB3L1', 'NKX3_2', 'HMGB3', 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 'HDAC8', 'PQBP1', 'FOXG1', 'HMGA2', 'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 'SP7', 'TRPS1', 'ZNF335', 'GSC')] <- 1
EntireCohort_pathvar$TFgene[EntireCohort_pathvar$GENE5 %in% c('ERCC2', 'CREB3L1', 'NKX3_2', 'HMGB3', 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 'HDAC8', 'PQBP1', 'FOXG1', 'HMGA2', 'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 'SP7', 'TRPS1', 'ZNF335', 'GSC')] <- 1
EntireCohort_pathvar$TFgene[EntireCohort_pathvar$GENE6 %in% c('ERCC2', 'CREB3L1', 'NKX3_2', 'HMGB3', 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 'HDAC8', 'PQBP1', 'FOXG1', 'HMGA2', 'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 'SP7', 'TRPS1', 'ZNF335', 'GSC')] <- 1
EntireCohort_pathvar$TFgene[EntireCohort_pathvar$GENE7 %in% c('ERCC2', 'CREB3L1', 'NKX3_2', 'HMGB3', 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 'HDAC8', 'PQBP1', 'FOXG1', 'HMGA2', 'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 'SP7', 'TRPS1', 'ZNF335', 'GSC')] <- 1
EntireCohort_pathvar$TFgene[EntireCohort_pathvar$GENE8 %in% c('ERCC2', 'CREB3L1', 'NKX3_2', 'HMGB3', 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 'HDAC8', 'PQBP1', 'FOXG1', 'HMGA2', 'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 'SP7', 'TRPS1', 'ZNF335', 'GSC')] <- 1
EntireCohort_pathvar$TFgene[EntireCohort_pathvar$GENE9 %in% c('ERCC2', 'CREB3L1', 'NKX3_2', 'HMGB3', 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 'HDAC8', 'PQBP1', 'FOXG1', 'HMGA2', 'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 'SP7', 'TRPS1', 'ZNF335', 'GSC')] <- 1


EntireCohort_pathvar$TFgene[is.na(EntireCohort_pathvar$TFgene)] <- 0
EntireCohort_pathvar$TFgene <- as.numeric(EntireCohort_pathvar$TFgene)
unique(EntireCohort_pathvar$TFgene)
# 0 1

table(EntireCohort_pathvar$TFgene)
#    0      1 
# 451428  18371 




# now to make a column for each of the TFgenes
### 'ERCC2'
EntireCohort_pathvar$ERCC2 <- NA
EntireCohort_pathvar$ERCC2[EntireCohort_pathvar$GENE1 %in% c("ERCC2")] <- 1
EntireCohort_pathvar$ERCC2[EntireCohort_pathvar$GENE2 %in% c("ERCC2")] <- 1
EntireCohort_pathvar$ERCC2[EntireCohort_pathvar$GENE3 %in% c("ERCC2")] <- 1
EntireCohort_pathvar$ERCC2[EntireCohort_pathvar$GENE4 %in% c("ERCC2")] <- 1
EntireCohort_pathvar$ERCC2[EntireCohort_pathvar$GENE5 %in% c("ERCC2")] <- 1
EntireCohort_pathvar$ERCC2[EntireCohort_pathvar$GENE6 %in% c("ERCC2")] <- 1
EntireCohort_pathvar$ERCC2[EntireCohort_pathvar$GENE7 %in% c("ERCC2")] <- 1
EntireCohort_pathvar$ERCC2[EntireCohort_pathvar$GENE8 %in% c("ERCC2")] <- 1
EntireCohort_pathvar$ERCC2[EntireCohort_pathvar$GENE9 %in% c("ERCC2")] <- 1

EntireCohort_pathvar$ERCC2[is.na(EntireCohort_pathvar$ERCC2)] <- 0
EntireCohort_pathvar$ERCC2 <- as.numeric(EntireCohort_pathvar$ERCC2)
unique(EntireCohort_pathvar$ERCC2)
# 0 1

table(EntireCohort_pathvar$ERCC2)
#  0      1 
# 464657   5142





### 'CREB3L1'
EntireCohort_pathvar$CREB3L1 <- NA
EntireCohort_pathvar$CREB3L1[EntireCohort_pathvar$GENE1 %in% c("CREB3L1")] <- 1
EntireCohort_pathvar$CREB3L1[EntireCohort_pathvar$GENE2 %in% c("CREB3L1")] <- 1
EntireCohort_pathvar$CREB3L1[EntireCohort_pathvar$GENE3 %in% c("CREB3L1")] <- 1
EntireCohort_pathvar$CREB3L1[EntireCohort_pathvar$GENE4 %in% c("CREB3L1")] <- 1
EntireCohort_pathvar$CREB3L1[EntireCohort_pathvar$GENE5 %in% c("CREB3L1")] <- 1
EntireCohort_pathvar$CREB3L1[EntireCohort_pathvar$GENE6 %in% c("CREB3L1")] <- 1
EntireCohort_pathvar$CREB3L1[EntireCohort_pathvar$GENE7 %in% c("CREB3L1")] <- 1
EntireCohort_pathvar$CREB3L1[EntireCohort_pathvar$GENE8 %in% c("CREB3L1")] <- 1
EntireCohort_pathvar$CREB3L1[EntireCohort_pathvar$GENE9 %in% c("CREB3L1")] <- 1

EntireCohort_pathvar$CREB3L1[is.na(EntireCohort_pathvar$CREB3L1)] <- 0
EntireCohort_pathvar$CREB3L1 <- as.numeric(EntireCohort_pathvar$CREB3L1)
unique(EntireCohort_pathvar$CREB3L1)
#  0 1

table(EntireCohort_pathvar$CREB3L1)
#  0      1 
# 468855    944




### 'NKX3_2'
EntireCohort_pathvar$NKX3_2 <- NA
EntireCohort_pathvar$NKX3_2[EntireCohort_pathvar$GENE1 %in% c("NKX3_2")] <- 1
EntireCohort_pathvar$NKX3_2[EntireCohort_pathvar$GENE2 %in% c("NKX3_2")] <- 1
EntireCohort_pathvar$NKX3_2[EntireCohort_pathvar$GENE3 %in% c("NKX3_2")] <- 1
EntireCohort_pathvar$NKX3_2[EntireCohort_pathvar$GENE4 %in% c("NKX3_2")] <- 1
EntireCohort_pathvar$NKX3_2[EntireCohort_pathvar$GENE5 %in% c("NKX3_2")] <- 1
EntireCohort_pathvar$NKX3_2[EntireCohort_pathvar$GENE6 %in% c("NKX3_2")] <- 1
EntireCohort_pathvar$NKX3_2[EntireCohort_pathvar$GENE7 %in% c("NKX3_2")] <- 1
EntireCohort_pathvar$NKX3_2[EntireCohort_pathvar$GENE8 %in% c("NKX3_2")] <- 1
EntireCohort_pathvar$NKX3_2[EntireCohort_pathvar$GENE9 %in% c("NKX3_2")] <- 1

EntireCohort_pathvar$NKX3_2[is.na(EntireCohort_pathvar$NKX3_2)] <- 0
EntireCohort_pathvar$NKX3_2 <- as.numeric(EntireCohort_pathvar$NKX3_2)
unique(EntireCohort_pathvar$NKX3_2)
#  0 1

table(EntireCohort_pathvar$NKX3_2)
#   0      1 
# 469294    505 




### 'ARX'
# no variants present




### 'HMGB3'
EntireCohort_pathvar$HMGB3 <- NA
EntireCohort_pathvar$HMGB3[EntireCohort_pathvar$GENE1 %in% c("HMGB3")] <- 1
EntireCohort_pathvar$HMGB3[EntireCohort_pathvar$GENE2 %in% c("HMGB3")] <- 1
EntireCohort_pathvar$HMGB3[EntireCohort_pathvar$GENE3 %in% c("HMGB3")] <- 1
EntireCohort_pathvar$HMGB3[EntireCohort_pathvar$GENE4 %in% c("HMGB3")] <- 1
EntireCohort_pathvar$HMGB3[EntireCohort_pathvar$GENE5 %in% c("HMGB3")] <- 1
EntireCohort_pathvar$HMGB3[EntireCohort_pathvar$GENE6 %in% c("HMGB3")] <- 1
EntireCohort_pathvar$HMGB3[EntireCohort_pathvar$GENE7 %in% c("HMGB3")] <- 1
EntireCohort_pathvar$HMGB3[EntireCohort_pathvar$GENE8 %in% c("HMGB3")] <- 1
EntireCohort_pathvar$HMGB3[EntireCohort_pathvar$GENE9 %in% c("HMGB3")] <- 1

EntireCohort_pathvar$HMGB3[is.na(EntireCohort_pathvar$HMGB3)] <- 0
EntireCohort_pathvar$HMGB3 <- as.numeric(EntireCohort_pathvar$HMGB3)
unique(EntireCohort_pathvar$HMGB3)
#  0 1

table(EntireCohort_pathvar$HMGB3)
#  0      1 
# 469742     57




### 'RAD21'
EntireCohort_pathvar$RAD21 <- NA
EntireCohort_pathvar$RAD21[EntireCohort_pathvar$GENE1 %in% c("RAD21")] <- 1
EntireCohort_pathvar$RAD21[EntireCohort_pathvar$GENE2 %in% c("RAD21")] <- 1
EntireCohort_pathvar$RAD21[EntireCohort_pathvar$GENE3 %in% c("RAD21")] <- 1
EntireCohort_pathvar$RAD21[EntireCohort_pathvar$GENE4 %in% c("RAD21")] <- 1
EntireCohort_pathvar$RAD21[EntireCohort_pathvar$GENE5 %in% c("RAD21")] <- 1
EntireCohort_pathvar$RAD21[EntireCohort_pathvar$GENE6 %in% c("RAD21")] <- 1
EntireCohort_pathvar$RAD21[EntireCohort_pathvar$GENE7 %in% c("RAD21")] <- 1
EntireCohort_pathvar$RAD21[EntireCohort_pathvar$GENE8 %in% c("RAD21")] <- 1
EntireCohort_pathvar$RAD21[EntireCohort_pathvar$GENE9 %in% c("RAD21")] <- 1

EntireCohort_pathvar$RAD21[is.na(EntireCohort_pathvar$RAD21)] <- 0
EntireCohort_pathvar$RAD21 <- as.numeric(EntireCohort_pathvar$RAD21)
unique(EntireCohort_pathvar$RAD21)
# 0 1

table(EntireCohort_pathvar$RAD21)
#  0      1 
# 468809    990 




### 'SMC1A'
EntireCohort_pathvar$SMC1A <- NA
EntireCohort_pathvar$SMC1A[EntireCohort_pathvar$GENE1 %in% c("SMC1A")] <- 1
EntireCohort_pathvar$SMC1A[EntireCohort_pathvar$GENE2 %in% c("SMC1A")] <- 1
EntireCohort_pathvar$SMC1A[EntireCohort_pathvar$GENE3 %in% c("SMC1A")] <- 1
EntireCohort_pathvar$SMC1A[EntireCohort_pathvar$GENE4 %in% c("SMC1A")] <- 1
EntireCohort_pathvar$SMC1A[EntireCohort_pathvar$GENE5 %in% c("SMC1A")] <- 1
EntireCohort_pathvar$SMC1A[EntireCohort_pathvar$GENE6 %in% c("SMC1A")] <- 1
EntireCohort_pathvar$SMC1A[EntireCohort_pathvar$GENE7 %in% c("SMC1A")] <- 1
EntireCohort_pathvar$SMC1A[EntireCohort_pathvar$GENE8 %in% c("SMC1A")] <- 1
EntireCohort_pathvar$SMC1A[EntireCohort_pathvar$GENE9 %in% c("SMC1A")] <- 1

EntireCohort_pathvar$SMC1A[is.na(EntireCohort_pathvar$SMC1A)] <- 0
EntireCohort_pathvar$SMC1A <- as.numeric(EntireCohort_pathvar$SMC1A)
unique(EntireCohort_pathvar$SMC1A)
# 0 1

table(EntireCohort_pathvar$SMC1A)
#  0      1 
# 469779     20




### 'SMC3'
EntireCohort_pathvar$SMC3 <- NA
EntireCohort_pathvar$SMC3[EntireCohort_pathvar$GENE1 %in% c("SMC3")] <- 1
EntireCohort_pathvar$SMC3[EntireCohort_pathvar$GENE2 %in% c("SMC3")] <- 1
EntireCohort_pathvar$SMC3[EntireCohort_pathvar$GENE3 %in% c("SMC3")] <- 1
EntireCohort_pathvar$SMC3[EntireCohort_pathvar$GENE4 %in% c("SMC3")] <- 1
EntireCohort_pathvar$SMC3[EntireCohort_pathvar$GENE5 %in% c("SMC3")] <- 1
EntireCohort_pathvar$SMC3[EntireCohort_pathvar$GENE6 %in% c("SMC3")] <- 1
EntireCohort_pathvar$SMC3[EntireCohort_pathvar$GENE7 %in% c("SMC3")] <- 1
EntireCohort_pathvar$SMC3[EntireCohort_pathvar$GENE8 %in% c("SMC3")] <- 1
EntireCohort_pathvar$SMC3[EntireCohort_pathvar$GENE9 %in% c("SMC3")] <- 1

EntireCohort_pathvar$SMC3[is.na(EntireCohort_pathvar$SMC3)] <- 0
EntireCohort_pathvar$SMC3 <- as.numeric(EntireCohort_pathvar$SMC3)
unique(EntireCohort_pathvar$SMC3)
# 0 1

table(EntireCohort_pathvar$SMC3)
#   0      1 
# 469469    330 




### 'NIPBL'
EntireCohort_pathvar$NIPBL <- NA
EntireCohort_pathvar$NIPBL[EntireCohort_pathvar$GENE1 %in% c("NIPBL")] <- 1
EntireCohort_pathvar$NIPBL[EntireCohort_pathvar$GENE2 %in% c("NIPBL")] <- 1
EntireCohort_pathvar$NIPBL[EntireCohort_pathvar$GENE3 %in% c("NIPBL")] <- 1
EntireCohort_pathvar$NIPBL[EntireCohort_pathvar$GENE4 %in% c("NIPBL")] <- 1
EntireCohort_pathvar$NIPBL[EntireCohort_pathvar$GENE5 %in% c("NIPBL")] <- 1
EntireCohort_pathvar$NIPBL[EntireCohort_pathvar$GENE6 %in% c("NIPBL")] <- 1
EntireCohort_pathvar$NIPBL[EntireCohort_pathvar$GENE7 %in% c("NIPBL")] <- 1
EntireCohort_pathvar$NIPBL[EntireCohort_pathvar$GENE8 %in% c("NIPBL")] <- 1
EntireCohort_pathvar$NIPBL[EntireCohort_pathvar$GENE9 %in% c("NIPBL")] <- 1

EntireCohort_pathvar$NIPBL[is.na(EntireCohort_pathvar$NIPBL)] <- 0
EntireCohort_pathvar$NIPBL <- as.numeric(EntireCohort_pathvar$NIPBL)
unique(EntireCohort_pathvar$NIPBL)
#  0 1

table(EntireCohort_pathvar$NIPBL)
#  0      1 
# 468536   1263 





### 'SMARCA2'
# no variants present




### 'SMARCE1'
# no varinats present





### 'SRCAP'
# no variants present





### 'DNMT3A'
# no variants present




### 'HDAC8'
EntireCohort_pathvar$HDAC8 <- NA
EntireCohort_pathvar$HDAC8[EntireCohort_pathvar$GENE1 %in% c("HDAC8")] <- 1
EntireCohort_pathvar$HDAC8[EntireCohort_pathvar$GENE2 %in% c("HDAC8")] <- 1
EntireCohort_pathvar$HDAC8[EntireCohort_pathvar$GENE3 %in% c("HDAC8")] <- 1
EntireCohort_pathvar$HDAC8[EntireCohort_pathvar$GENE4 %in% c("HDAC8")] <- 1
EntireCohort_pathvar$HDAC8[EntireCohort_pathvar$GENE5 %in% c("HDAC8")] <- 1
EntireCohort_pathvar$HDAC8[EntireCohort_pathvar$GENE6 %in% c("HDAC8")] <- 1
EntireCohort_pathvar$HDAC8[EntireCohort_pathvar$GENE7 %in% c("HDAC8")] <- 1
EntireCohort_pathvar$HDAC8[EntireCohort_pathvar$GENE8 %in% c("HDAC8")] <- 1
EntireCohort_pathvar$HDAC8[EntireCohort_pathvar$GENE9 %in% c("HDAC8")] <- 1

EntireCohort_pathvar$HDAC8[is.na(EntireCohort_pathvar$HDAC8)] <- 0
EntireCohort_pathvar$HDAC8 <- as.numeric(EntireCohort_pathvar$HDAC8)
unique(EntireCohort_pathvar$HDAC8)
#  0 1

table(EntireCohort_pathvar$HDAC8)
#  0      1 
# 469756     43




### 'PQBP1'
EntireCohort_pathvar$PQBP1 <- NA
EntireCohort_pathvar$PQBP1[EntireCohort_pathvar$GENE1 %in% c("PQBP1")] <- 1
EntireCohort_pathvar$PQBP1[EntireCohort_pathvar$GENE2 %in% c("PQBP1")] <- 1
EntireCohort_pathvar$PQBP1[EntireCohort_pathvar$GENE3 %in% c("PQBP1")] <- 1
EntireCohort_pathvar$PQBP1[EntireCohort_pathvar$GENE4 %in% c("PQBP1")] <- 1
EntireCohort_pathvar$PQBP1[EntireCohort_pathvar$GENE5 %in% c("PQBP1")] <- 1
EntireCohort_pathvar$PQBP1[EntireCohort_pathvar$GENE6 %in% c("PQBP1")] <- 1
EntireCohort_pathvar$PQBP1[EntireCohort_pathvar$GENE7 %in% c("PQBP1")] <- 1
EntireCohort_pathvar$PQBP1[EntireCohort_pathvar$GENE8 %in% c("PQBP1")] <- 1
EntireCohort_pathvar$PQBP1[EntireCohort_pathvar$GENE9 %in% c("PQBP1")] <- 1

EntireCohort_pathvar$PQBP1[is.na(EntireCohort_pathvar$PQBP1)] <- 0
EntireCohort_pathvar$PQBP1 <- as.numeric(EntireCohort_pathvar$PQBP1)
unique(EntireCohort_pathvar$PQBP1)
#  0 1

table(EntireCohort_pathvar$PQBP1)
#  0      1 
# 469674    125 




### 'FOXG1'
EntireCohort_pathvar$FOXG1 <- NA
EntireCohort_pathvar$FOXG1[EntireCohort_pathvar$GENE1 %in% c("FOXG1")] <- 1
EntireCohort_pathvar$FOXG1[EntireCohort_pathvar$GENE2 %in% c("FOXG1")] <- 1
EntireCohort_pathvar$FOXG1[EntireCohort_pathvar$GENE3 %in% c("FOXG1")] <- 1
EntireCohort_pathvar$FOXG1[EntireCohort_pathvar$GENE4 %in% c("FOXG1")] <- 1
EntireCohort_pathvar$FOXG1[EntireCohort_pathvar$GENE5 %in% c("FOXG1")] <- 1
EntireCohort_pathvar$FOXG1[EntireCohort_pathvar$GENE6 %in% c("FOXG1")] <- 1
EntireCohort_pathvar$FOXG1[EntireCohort_pathvar$GENE7 %in% c("FOXG1")] <- 1
EntireCohort_pathvar$FOXG1[EntireCohort_pathvar$GENE8 %in% c("FOXG1")] <- 1
EntireCohort_pathvar$FOXG1[EntireCohort_pathvar$GENE9 %in% c("FOXG1")] <- 1

EntireCohort_pathvar$FOXG1[is.na(EntireCohort_pathvar$FOXG1)] <- 0
EntireCohort_pathvar$FOXG1 <- as.numeric(EntireCohort_pathvar$FOXG1)
unique(EntireCohort_pathvar$FOXG1)
#  0 1

table(EntireCohort_pathvar$FOXG1)
#  0      1 
# 469694    105 




### 'HMGA2' 
EntireCohort_pathvar$HMGA2 <- NA
EntireCohort_pathvar$HMGA2[EntireCohort_pathvar$GENE1 %in% c("HMGA2")] <- 1
EntireCohort_pathvar$HMGA2[EntireCohort_pathvar$GENE2 %in% c("HMGA2")] <- 1
EntireCohort_pathvar$HMGA2[EntireCohort_pathvar$GENE3 %in% c("HMGA2")] <- 1
EntireCohort_pathvar$HMGA2[EntireCohort_pathvar$GENE4 %in% c("HMGA2")] <- 1
EntireCohort_pathvar$HMGA2[EntireCohort_pathvar$GENE5 %in% c("HMGA2")] <- 1
EntireCohort_pathvar$HMGA2[EntireCohort_pathvar$GENE6 %in% c("HMGA2")] <- 1
EntireCohort_pathvar$HMGA2[EntireCohort_pathvar$GENE7 %in% c("HMGA2")] <- 1
EntireCohort_pathvar$HMGA2[EntireCohort_pathvar$GENE8 %in% c("HMGA2")] <- 1
EntireCohort_pathvar$HMGA2[EntireCohort_pathvar$GENE9 %in% c("HMGA2")] <- 1

EntireCohort_pathvar$HMGA2[is.na(EntireCohort_pathvar$HMGA2)] <- 0
EntireCohort_pathvar$HMGA2 <- as.numeric(EntireCohort_pathvar$HMGA2)
unique(EntireCohort_pathvar$HMGA2)
# 0 1

table(EntireCohort_pathvar$HMGA2)
#  0      1 
# 469728     71 




### 'LARP7'
EntireCohort_pathvar$LARP7 <- NA
EntireCohort_pathvar$LARP7[EntireCohort_pathvar$GENE1 %in% c("LARP7")] <- 1
EntireCohort_pathvar$LARP7[EntireCohort_pathvar$GENE2 %in% c("LARP7")] <- 1
EntireCohort_pathvar$LARP7[EntireCohort_pathvar$GENE3 %in% c("LARP7")] <- 1
EntireCohort_pathvar$LARP7[EntireCohort_pathvar$GENE4 %in% c("LARP7")] <- 1
EntireCohort_pathvar$LARP7[EntireCohort_pathvar$GENE5 %in% c("LARP7")] <- 1
EntireCohort_pathvar$LARP7[EntireCohort_pathvar$GENE6 %in% c("LARP7")] <- 1
EntireCohort_pathvar$LARP7[EntireCohort_pathvar$GENE7 %in% c("LARP7")] <- 1
EntireCohort_pathvar$LARP7[EntireCohort_pathvar$GENE8 %in% c("LARP7")] <- 1
EntireCohort_pathvar$LARP7[EntireCohort_pathvar$GENE9 %in% c("LARP7")] <- 1

EntireCohort_pathvar$LARP7[is.na(EntireCohort_pathvar$LARP7)] <- 0
EntireCohort_pathvar$LARP7 <- as.numeric(EntireCohort_pathvar$LARP7)
unique(EntireCohort_pathvar$LARP7)
# 0 1

table(EntireCohort_pathvar$LARP7)
#  0      1 
# 469202    597 




### 'SHOX' 
EntireCohort_pathvar$SHOX <- NA
EntireCohort_pathvar$SHOX[EntireCohort_pathvar$GENE1 %in% c("SHOX")] <- 1
EntireCohort_pathvar$SHOX[EntireCohort_pathvar$GENE2 %in% c("SHOX")] <- 1
EntireCohort_pathvar$SHOX[EntireCohort_pathvar$GENE3 %in% c("SHOX")] <- 1
EntireCohort_pathvar$SHOX[EntireCohort_pathvar$GENE4 %in% c("SHOX")] <- 1
EntireCohort_pathvar$SHOX[EntireCohort_pathvar$GENE5 %in% c("SHOX")] <- 1
EntireCohort_pathvar$SHOX[EntireCohort_pathvar$GENE6 %in% c("SHOX")] <- 1
EntireCohort_pathvar$SHOX[EntireCohort_pathvar$GENE7 %in% c("SHOX")] <- 1
EntireCohort_pathvar$SHOX[EntireCohort_pathvar$GENE8 %in% c("SHOX")] <- 1
EntireCohort_pathvar$SHOX[EntireCohort_pathvar$GENE9 %in% c("SHOX")] <- 1

EntireCohort_pathvar$SHOX[is.na(EntireCohort_pathvar$SHOX)] <- 0
EntireCohort_pathvar$SHOX <- as.numeric(EntireCohort_pathvar$SHOX)
unique(EntireCohort_pathvar$SHOX)
# 0 1

table(EntireCohort_pathvar$SHOX)
#  0      1 
# 469672    127




### 'SOX11'
EntireCohort_pathvar$SOX11 <- NA
EntireCohort_pathvar$SOX11[EntireCohort_pathvar$GENE1 %in% c("SOX11")] <- 1
EntireCohort_pathvar$SOX11[EntireCohort_pathvar$GENE2 %in% c("SOX11")] <- 1
EntireCohort_pathvar$SOX11[EntireCohort_pathvar$GENE3 %in% c("SOX11")] <- 1
EntireCohort_pathvar$SOX11[EntireCohort_pathvar$GENE4 %in% c("SOX11")] <- 1
EntireCohort_pathvar$SOX11[EntireCohort_pathvar$GENE5 %in% c("SOX11")] <- 1
EntireCohort_pathvar$SOX11[EntireCohort_pathvar$GENE6 %in% c("SOX11")] <- 1
EntireCohort_pathvar$SOX11[EntireCohort_pathvar$GENE7 %in% c("SOX11")] <- 1
EntireCohort_pathvar$SOX11[EntireCohort_pathvar$GENE8 %in% c("SOX11")] <- 1
EntireCohort_pathvar$SOX11[EntireCohort_pathvar$GENE9 %in% c("SOX11")] <- 1

EntireCohort_pathvar$SOX11[is.na(EntireCohort_pathvar$SOX11)] <- 0
EntireCohort_pathvar$SOX11 <- as.numeric(EntireCohort_pathvar$SOX11)
unique(EntireCohort_pathvar$SOX11)
#  0 1

table(EntireCohort_pathvar$SOX11)
#  0      1 
# 469473    326 





### 'TAF13'    
EntireCohort_pathvar$TAF13 <- NA
EntireCohort_pathvar$TAF13[EntireCohort_pathvar$GENE1 %in% c("TAF13")] <- 1
EntireCohort_pathvar$TAF13[EntireCohort_pathvar$GENE2 %in% c("TAF13")] <- 1
EntireCohort_pathvar$TAF13[EntireCohort_pathvar$GENE3 %in% c("TAF13")] <- 1
EntireCohort_pathvar$TAF13[EntireCohort_pathvar$GENE4 %in% c("TAF13")] <- 1
EntireCohort_pathvar$TAF13[EntireCohort_pathvar$GENE5 %in% c("TAF13")] <- 1
EntireCohort_pathvar$TAF13[EntireCohort_pathvar$GENE6 %in% c("TAF13")] <- 1
EntireCohort_pathvar$TAF13[EntireCohort_pathvar$GENE7 %in% c("TAF13")] <- 1
EntireCohort_pathvar$TAF13[EntireCohort_pathvar$GENE8 %in% c("TAF13")] <- 1
EntireCohort_pathvar$TAF13[EntireCohort_pathvar$GENE9 %in% c("TAF13")] <- 1

EntireCohort_pathvar$TAF13[is.na(EntireCohort_pathvar$TAF13)] <- 0
EntireCohort_pathvar$TAF13 <- as.numeric(EntireCohort_pathvar$TAF13)
unique(EntireCohort_pathvar$TAF13)
#  0 1

table(EntireCohort_pathvar$TAF13)
#  0      1 
# 469387    412 




### 'TBX15'
EntireCohort_pathvar$TBX15 <- NA
EntireCohort_pathvar$TBX15[EntireCohort_pathvar$GENE1 %in% c("TBX15")] <- 1
EntireCohort_pathvar$TBX15[EntireCohort_pathvar$GENE2 %in% c("TBX15")] <- 1
EntireCohort_pathvar$TBX15[EntireCohort_pathvar$GENE3 %in% c("TBX15")] <- 1
EntireCohort_pathvar$TBX15[EntireCohort_pathvar$GENE4 %in% c("TBX15")] <- 1
EntireCohort_pathvar$TBX15[EntireCohort_pathvar$GENE5 %in% c("TBX15")] <- 1
EntireCohort_pathvar$TBX15[EntireCohort_pathvar$GENE6 %in% c("TBX15")] <- 1
EntireCohort_pathvar$TBX15[EntireCohort_pathvar$GENE7 %in% c("TBX15")] <- 1
EntireCohort_pathvar$TBX15[EntireCohort_pathvar$GENE8 %in% c("TBX15")] <- 1
EntireCohort_pathvar$TBX15[EntireCohort_pathvar$GENE9 %in% c("TBX15")] <- 1

EntireCohort_pathvar$TBX15[is.na(EntireCohort_pathvar$TBX15)] <- 0
EntireCohort_pathvar$TBX15 <- as.numeric(EntireCohort_pathvar$TBX15)
unique(EntireCohort_pathvar$TBX15)
#  0 1

table(EntireCohort_pathvar$TBX15)
#  0      1 
# 468209   1590 




### 'TBX3' 
EntireCohort_pathvar$TBX3 <- NA
EntireCohort_pathvar$TBX3[EntireCohort_pathvar$GENE1 %in% c("TBX3")] <- 1
EntireCohort_pathvar$TBX3[EntireCohort_pathvar$GENE2 %in% c("TBX3")] <- 1
EntireCohort_pathvar$TBX3[EntireCohort_pathvar$GENE3 %in% c("TBX3")] <- 1
EntireCohort_pathvar$TBX3[EntireCohort_pathvar$GENE4 %in% c("TBX3")] <- 1
EntireCohort_pathvar$TBX3[EntireCohort_pathvar$GENE5 %in% c("TBX3")] <- 1
EntireCohort_pathvar$TBX3[EntireCohort_pathvar$GENE6 %in% c("TBX3")] <- 1
EntireCohort_pathvar$TBX3[EntireCohort_pathvar$GENE7 %in% c("TBX3")] <- 1
EntireCohort_pathvar$TBX3[EntireCohort_pathvar$GENE8 %in% c("TBX3")] <- 1
EntireCohort_pathvar$TBX3[EntireCohort_pathvar$GENE9 %in% c("TBX3")] <- 1

EntireCohort_pathvar$TBX3[is.na(EntireCohort_pathvar$TBX3)] <- 0
EntireCohort_pathvar$TBX3 <- as.numeric(EntireCohort_pathvar$TBX3)
unique(EntireCohort_pathvar$TBX3)
#  0 1

table(EntireCohort_pathvar$TBX3)
#  0      1 
# 468697   1102




### 'SP7' 
EntireCohort_pathvar$SP7 <- NA
EntireCohort_pathvar$SP7[EntireCohort_pathvar$GENE1 %in% c("SP7")] <- 1
EntireCohort_pathvar$SP7[EntireCohort_pathvar$GENE2 %in% c("SP7")] <- 1
EntireCohort_pathvar$SP7[EntireCohort_pathvar$GENE3 %in% c("SP7")] <- 1
EntireCohort_pathvar$SP7[EntireCohort_pathvar$GENE4 %in% c("SP7")] <- 1
EntireCohort_pathvar$SP7[EntireCohort_pathvar$GENE5 %in% c("SP7")] <- 1
EntireCohort_pathvar$SP7[EntireCohort_pathvar$GENE6 %in% c("SP7")] <- 1
EntireCohort_pathvar$SP7[EntireCohort_pathvar$GENE7 %in% c("SP7")] <- 1
EntireCohort_pathvar$SP7[EntireCohort_pathvar$GENE8 %in% c("SP7")] <- 1
EntireCohort_pathvar$SP7[EntireCohort_pathvar$GENE9 %in% c("SP7")] <- 1

EntireCohort_pathvar$SP7[is.na(EntireCohort_pathvar$SP7)] <- 0
EntireCohort_pathvar$SP7 <- as.numeric(EntireCohort_pathvar$SP7)
unique(EntireCohort_pathvar$SP7)
#  0 1

table(EntireCohort_pathvar$SP7)
#  0      1 
# 469529    270




### 'TRPS1' 
EntireCohort_pathvar$TRPS1 <- NA
EntireCohort_pathvar$TRPS1[EntireCohort_pathvar$GENE1 %in% c("TRPS1")] <- 1
EntireCohort_pathvar$TRPS1[EntireCohort_pathvar$GENE2 %in% c("TRPS1")] <- 1
EntireCohort_pathvar$TRPS1[EntireCohort_pathvar$GENE3 %in% c("TRPS1")] <- 1
EntireCohort_pathvar$TRPS1[EntireCohort_pathvar$GENE4 %in% c("TRPS1")] <- 1
EntireCohort_pathvar$TRPS1[EntireCohort_pathvar$GENE5 %in% c("TRPS1")] <- 1
EntireCohort_pathvar$TRPS1[EntireCohort_pathvar$GENE6 %in% c("TRPS1")] <- 1
EntireCohort_pathvar$TRPS1[EntireCohort_pathvar$GENE7 %in% c("TRPS1")] <- 1
EntireCohort_pathvar$TRPS1[EntireCohort_pathvar$GENE8 %in% c("TRPS1")] <- 1
EntireCohort_pathvar$TRPS1[EntireCohort_pathvar$GENE9 %in% c("TRPS1")] <- 1

EntireCohort_pathvar$TRPS1[is.na(EntireCohort_pathvar$TRPS1)] <- 0
EntireCohort_pathvar$TRPS1 <- as.numeric(EntireCohort_pathvar$TRPS1)
unique(EntireCohort_pathvar$TRPS1)
# 0 1

table(EntireCohort_pathvar$TRPS1)
#  0      1 
# 467669   2130 




### 'ZNF335' 
EntireCohort_pathvar$ZNF335 <- NA
EntireCohort_pathvar$ZNF335[EntireCohort_pathvar$GENE1 %in% c("ZNF335")] <- 1
EntireCohort_pathvar$ZNF335[EntireCohort_pathvar$GENE2 %in% c("ZNF335")] <- 1
EntireCohort_pathvar$ZNF335[EntireCohort_pathvar$GENE3 %in% c("ZNF335")] <- 1
EntireCohort_pathvar$ZNF335[EntireCohort_pathvar$GENE4 %in% c("ZNF335")] <- 1
EntireCohort_pathvar$ZNF335[EntireCohort_pathvar$GENE5 %in% c("ZNF335")] <- 1
EntireCohort_pathvar$ZNF335[EntireCohort_pathvar$GENE6 %in% c("ZNF335")] <- 1
EntireCohort_pathvar$ZNF335[EntireCohort_pathvar$GENE7 %in% c("ZNF335")] <- 1
EntireCohort_pathvar$ZNF335[EntireCohort_pathvar$GENE8 %in% c("ZNF335")] <- 1
EntireCohort_pathvar$ZNF335[EntireCohort_pathvar$GENE9 %in% c("ZNF335")] <- 1

EntireCohort_pathvar$ZNF335[is.na(EntireCohort_pathvar$ZNF335)] <- 0
EntireCohort_pathvar$ZNF335 <- as.numeric(EntireCohort_pathvar$ZNF335)
unique(EntireCohort_pathvar$ZNF335)
# 0 1

table(EntireCohort_pathvar$ZNF335)
#  0      1 
# 467799   2000





### 'GSC' 
EntireCohort_pathvar$GSC <- NA
EntireCohort_pathvar$GSC[EntireCohort_pathvar$GENE1 %in% c("GSC")] <- 1
EntireCohort_pathvar$GSC[EntireCohort_pathvar$GENE2 %in% c("GSC")] <- 1
EntireCohort_pathvar$GSC[EntireCohort_pathvar$GENE3 %in% c("GSC")] <- 1
EntireCohort_pathvar$GSC[EntireCohort_pathvar$GENE4 %in% c("GSC")] <- 1
EntireCohort_pathvar$GSC[EntireCohort_pathvar$GENE5 %in% c("GSC")] <- 1
EntireCohort_pathvar$GSC[EntireCohort_pathvar$GENE6 %in% c("GSC")] <- 1
EntireCohort_pathvar$GSC[EntireCohort_pathvar$GENE7 %in% c("GSC")] <- 1
EntireCohort_pathvar$GSC[EntireCohort_pathvar$GENE8 %in% c("GSC")] <- 1
EntireCohort_pathvar$GSC[EntireCohort_pathvar$GENE9 %in% c("GSC")] <- 1

EntireCohort_pathvar$GSC[is.na(EntireCohort_pathvar$GSC)] <- 0
EntireCohort_pathvar$GSC <- as.numeric(EntireCohort_pathvar$GSC)
unique(EntireCohort_pathvar$GSC)
#  0 1

table(EntireCohort_pathvar$GSC)
#  0      1 
# 469267    532 





# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)

















# TFgene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_TFgene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_TFgene_, "pathvar_TFgene_.csv", row.names=FALSE)
# system('dx upload pathvar_TFgene_.csv --path Emily-folder/results_2/pathvar_TFgene_.csv')


# TFgene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ TFgene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_TFgene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_incidence_1, "pathvar_TFgene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_incidence_1.csv --path Emily-folder/results_2/pathvar_TFgene_incidence_1.csv')


# GOIs
LR_2 <- glm(formula = cancer_ROOC_dummy ~ ERCC2+CREB3L1+NKX3_2+HMGB3+RAD21+SMC1A+SMC3+NIPBL+HDAC8+PQBP1+FOXG1+HMGA2+LARP7+SHOX+SOX11+TAF13+TBX15+TBX3+SP7+TRPS1+ZNF335+GSC, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_TFgene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_incidence_2, "pathvar_TFgene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_incidence_2.csv --path Emily-folder/results_2/pathvar_TFgene_incidence_2.csv')







### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "TFgene", 'ERCC2', 'CREB3L1', 'NKX3_2', 
                                     'HMGB3', 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 'HDAC8', 'PQBP1', 
                                     'FOXG1', 'HMGA2', 'LARP7', 'SHOX', 'SOX11', 'TAF13', 
                                     'TBX15', 'TBX3', 'SP7', 'TRPS1', 'ZNF335', 'GSC')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     29

with(data_setA, table(GSC))
# HMGB3 = 56
# SMC1A = 20
# HDAC8 = 42
# HMGA2 = 70




# TFgene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ TFgene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_TFgene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_incidence_3, "pathvar_TFgene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_incidence_3.csv --path Emily-folder/results_2/pathvar_TFgene_incidence_3.csv')


# age + ethnicity + smoking + sex + TFgene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_TFgene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_incidence_4, "pathvar_TFgene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_incidence_4.csv --path Emily-folder/results_2/pathvar_TFgene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ ERCC2+CREB3L1+NKX3_2+HMGB3+RAD21+SMC1A+SMC3+NIPBL+HDAC8+PQBP1+FOXG1+HMGA2+LARP7+SHOX+SOX11+TAF13+TBX15+TBX3+SP7+TRPS1+ZNF335+GSC, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_TFgene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_incidence_5, "pathvar_TFgene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_incidence_5.csv --path Emily-folder/results_2/pathvar_TFgene_incidence_5.csv')


# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC2+CREB3L1+NKX3_2+HMGB3+RAD21+SMC1A+SMC3+NIPBL+HDAC8+PQBP1+FOXG1+HMGA2+LARP7+SHOX+SOX11+TAF13+TBX15+TBX3+SP7+TRPS1+ZNF335+GSC, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_TFgene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_incidence_6, "pathvar_TFgene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_incidence_6.csv --path Emily-folder/results_2/pathvar_TFgene_incidence_6.csv')




### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "TFgene", 'ERCC2', 'CREB3L1', 'NKX3_2', 'HMGB3', 
                                     'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 
                                     'HDAC8', 'PQBP1', 'FOXG1', 'HMGA2', 
                                     'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 
                                     'SP7', 'TRPS1', 'ZNF335', 'GSC')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     34

with(data_setB, table(GSC))
# HMGB3 = 40
# SMC1A = 15
# HDAC8 = 34
# PQBP1 = 95
# FOXG1 = 82
# HMGA2 = 52




# TFgene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ TFgene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_TFgene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_incidence_7, "pathvar_TFgene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_incidence_7.csv --path Emily-folder/results_2/pathvar_TFgene_incidence_7.csv')


# all variables + TFgene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+TFgene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_TFgene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_incidence_8, "pathvar_TFgene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_incidence_8.csv --path Emily-folder/results_2/pathvar_TFgene_incidence_8.csv')


# GOIs alone
LR_9 <- glm(formula = cancer_ROOC_dummy ~ ERCC2+CREB3L1+NKX3_2+HMGB3+RAD21+SMC1A+SMC3+NIPBL+HDAC8+PQBP1+FOXG1+HMGA2+LARP7+SHOX+SOX11+TAF13+TBX15+TBX3+SP7+TRPS1+ZNF335+GSC, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_TFgene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_incidence_9, "pathvar_TFgene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_incidence_9.csv --path Emily-folder/results_2/pathvar_TFgene_incidence_9.csv')


# all variables + GOIs
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+ERCC2+CREB3L1+NKX3_2+HMGB3+RAD21+SMC1A+SMC3+NIPBL+HDAC8+PQBP1+FOXG1+HMGA2+LARP7+SHOX+SOX11+TAF13+TBX15+TBX3+SP7+TRPS1+ZNF335+GSC, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_TFgene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_incidence_10, "pathvar_TFgene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_incidence_10.csv --path Emily-folder/results_2/pathvar_TFgene_incidence_10.csv')










### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(TFgene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(TFgene, cancer_ROOC_dummy)))
# p-value = 0.9929

with(EntireCohort_pathvar, table(TFgene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(TFgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5702495 0.6359565
# sample estimates:
#   odds ratio 
# 0.6020417 

with(EntireCohort_pathvar, table(TFgene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(TFgene, Ever_smoked_dummy)))
# p-value = 0.02411
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9372086 0.9955373
# sample estimates:
#   odds ratio 
# 0.9658962 

with(EntireCohort_pathvar, table(TFgene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(TFgene, Sex_Male)))
# p-value = 0.4871

with(EntireCohort_pathvar, table(TFgene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(TFgene, University_education)))
# p-value = 0.6803

with(EntireCohort_pathvar, table(TFgene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(TFgene, Annual_income))	)
# p-value = 0.2746

with(EntireCohort_pathvar, table(TFgene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(TFgene, Smoker_previous)))
# p-value = 0.006234
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9276060 0.9877475
# sample estimates:
#   odds ratio 
# 0.9572453 

with(EntireCohort_pathvar, table(TFgene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(TFgene, Smoker_current))	)
# p-value = 0.9315



### data_setA
with(data_setA, table(TFgene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(TFgene, cancer_ROOC_dummy)))
# p-value = 0.9291

with(data_setA, table(TFgene, Ethnicity_White))			
fisher.test(with(data_setA, table(TFgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5673680 0.6329686
# sample estimates:
#   odds ratio 
# 0.5991246 

with(data_setA, table(TFgene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(TFgene, Ever_smoked_dummy)))
# p-value = 0.02296
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9367777 0.9951697
# sample estimates:
#   odds ratio 
# 0.9655025

with(data_setA, table(TFgene, Sex_Male))			
fisher.test(with(data_setA, table(TFgene, Sex_Male)))
# p-value = 0.5093



### data_setB
with(data_setB, table(TFgene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(TFgene, cancer_ROOC_dummy)))
# p-value = 0.9366

with(data_setB, table(TFgene, Ethnicity_White))			
fisher.test(with(data_setB, table(TFgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5784216 0.6607484
# sample estimates:
#   odds ratio 
# 0.6179541

with(data_setB, table(TFgene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(TFgene, Ever_smoked_dummy)))
# p-value = 0.1696

with(data_setB, table(TFgene, Sex_Male))			
fisher.test(with(data_setB, table(TFgene, Sex_Male)))
# p-value = 0.6868

with(data_setB, table(TFgene, University_education))			
fisher.test(with(data_setB, table(TFgene, University_education)))
# p-value = 0.9163

with(data_setB, table(TFgene, Annual_income))			
fisher.test(with(data_setB, table(TFgene, Annual_income)))
# p-value = 0.4188






### can plot the continuous variables: 


# make TFgene into factor
df_for_plots <- EntireCohort_pathvar
df_for_plots$TFgene <- as.factor(df_for_plots$TFgene)

# try a violin plot

# Age_at_recruitment  
df_for_plots %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ TFgene, data = df_for_plots)
# W = 4215078102, p-value = 0.0001427

# BMI  
df_for_plots %>% 
  ggplot(mapping = aes(x = TFgene, y = BMI)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = TFgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ TFgene, data = df_for_plots)
# W = 4098916592, p-value = 0.6692

# Standing_height 
df_for_plots %>% 
  ggplot(mapping = aes(x = TFgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = TFgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ TFgene, data = df_for_plots)
# W = 4178092457, p-value = 0.0005548

# make this a bar plot
# Weekly_exercise
df_for_plots %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(TFgene), scales = "free_y") +
  theme_light()
# 



### data_setA

# make TFgene into factor
data_setA$TFgene <- as.factor(data_setA$TFgene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ TFgene, data = data_setA)
# W = 4.145e+09, p-value = 0.0001346



### data_setB

# make TFgene into factor
data_setB$TFgene <- as.factor(data_setB$TFgene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ TFgene, data = data_setB)
# W = 2745707127, p-value = 0.001603

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = TFgene, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = TFgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ TFgene, data = data_setB)
# W = 2.707e+09, p-value = 0.8434

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = TFgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = TFgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ TFgene, data = data_setB)
# W = 2739806210, p-value = 0.006846


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(TFgene), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# TFgene Type ------------------------------------------------------------

# make sure the execute the code in the 'TFgenes' section

## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_TFgene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_type_1, "pathvar_TFgene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_type_1.csv --path Emily-folder/results_2/pathvar_TFgene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_TFgene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_type_2, "pathvar_TFgene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_type_2.csv --path Emily-folder/results_2/pathvar_TFgene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_TFgene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_type_3, "pathvar_TFgene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_type_3.csv --path Emily-folder/results_2/pathvar_TFgene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_TFgene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_type_4, "pathvar_TFgene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_type_4.csv --path Emily-folder/results_2/pathvar_TFgene_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_TFgene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_type_5, "pathvar_TFgene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_type_5.csv --path Emily-folder/results_2/pathvar_TFgene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_TFgene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_type_6, "pathvar_TFgene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_type_6.csv --path Emily-folder/results_2/pathvar_TFgene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_TFgene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_type_7, "pathvar_TFgene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_type_7.csv --path Emily-folder/results_2/pathvar_TFgene_type_7.csv')





LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_TFgene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_type_8, "pathvar_TFgene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_type_8.csv --path Emily-folder/results_2/pathvar_TFgene_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_TFgene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_type_9, "pathvar_TFgene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_type_9.csv --path Emily-folder/results_2/pathvar_TFgene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_TFgene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_type_10, "pathvar_TFgene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_type_10.csv --path Emily-folder/results_2/pathvar_TFgene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_TFgene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_type_11, "pathvar_TFgene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_type_11.csv --path Emily-folder/results_2/pathvar_TFgene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_TFgene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_type_12, "pathvar_TFgene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_type_12.csv --path Emily-folder/results_2/pathvar_TFgene_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_TFgene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_type_13, "pathvar_TFgene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_type_13.csv --path Emily-folder/results_2/pathvar_TFgene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_TFgene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_type_14, "pathvar_TFgene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_type_14.csv --path Emily-folder/results_2/pathvar_TFgene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_TFgene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TFgene_type_15, "pathvar_TFgene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_type_15.csv --path Emily-folder/results_2/pathvar_TFgene_type_15.csv')













# SHOX Type -------------------------------------------------------------
# repeat the type analysis for just SHOX



LR_oral_SHOX <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SHOX, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_SHOX)
# exp(cbind(OR = coef(LR_oral_SHOX), confint(LR_oral_SHOX)))
pathvar_SHOX_type_1 <- tidy(LR_oral_SHOX, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SHOX_type_1, "pathvar_SHOX_type_1.csv", row.names=FALSE)
system('dx upload pathvar_SHOX_type_1.csv --path Emily-folder/results_2/pathvar_SHOX_type_1.csv')



LR_digestive_SHOX <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SHOX, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_SHOX)
# exp(cbind(OR = coef(LR_digestive_SHOX), confint(LR_digestive_SHOX)))
pathvar_SHOX_type_2 <- tidy(LR_digestive_SHOX, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SHOX_type_2, "pathvar_SHOX_type_2.csv", row.names=FALSE)
system('dx upload pathvar_SHOX_type_2.csv --path Emily-folder/results_2/pathvar_SHOX_type_2.csv')



LR_respiratory_SHOX <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SHOX, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_SHOX)
# exp(cbind(OR = coef(LR_respiratory_SHOX), confint(LR_respiratory_SHOX)))
pathvar_SHOX_type_3 <- tidy(LR_respiratory_SHOX, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SHOX_type_3, "pathvar_SHOX_type_3.csv", row.names=FALSE)
system('dx upload pathvar_SHOX_type_3.csv --path Emily-folder/results_2/pathvar_SHOX_type_3.csv')



LR_bone_SHOX <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SHOX, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_SHOX)
# exp(cbind(OR = coef(LR_bone_SHOX), confint(LR_bone_SHOX)))
pathvar_SHOX_type_4 <- tidy(LR_bone_SHOX, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SHOX_type_4, "pathvar_SHOX_type_4.csv", row.names=FALSE)
system('dx upload pathvar_SHOX_type_4.csv --path Emily-folder/results_2/pathvar_SHOX_type_4.csv')



LR_skin_SHOX <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SHOX, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_SHOX)
# exp(cbind(OR = coef(LR_skin_SHOX), confint(LR_skin_SHOX)))
pathvar_SHOX_type_5 <- tidy(LR_skin_SHOX, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SHOX_type_5, "pathvar_SHOX_type_5.csv", row.names=FALSE)
system('dx upload pathvar_SHOX_type_5.csv --path Emily-folder/results_2/pathvar_SHOX_type_5.csv')



LR_mesothelium_SHOX <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SHOX, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_SHOX)
# exp(cbind(OR = coef(LR_mesothelium_SHOX), confint(LR_mesothelium_SHOX)))
pathvar_SHOX_type_6 <- tidy(LR_mesothelium_SHOX, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SHOX_type_6, "pathvar_SHOX_type_6.csv", row.names=FALSE)
system('dx upload pathvar_SHOX_type_6.csv --path Emily-folder/results_2/pathvar_SHOX_type_6.csv')



LR_breast_SHOX <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SHOX, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_SHOX)
# exp(cbind(OR = coef(LR_breast_SHOX), confint(LR_breast_SHOX)))
pathvar_SHOX_type_7 <- tidy(LR_breast_SHOX, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SHOX_type_7, "pathvar_SHOX_type_7.csv", row.names=FALSE)
system('dx upload pathvar_SHOX_type_7.csv --path Emily-folder/results_2/pathvar_SHOX_type_7.csv')



LR_femalegenital_SHOX <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SHOX, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_SHOX)
# exp(cbind(OR = coef(LR_femalegenital_SHOX), confint(LR_femalegenital_SHOX)))
pathvar_SHOX_type_8 <- tidy(LR_femalegenital_SHOX, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SHOX_type_8, "pathvar_SHOX_type_8.csv", row.names=FALSE)
system('dx upload pathvar_SHOX_type_8.csv --path Emily-folder/results_2/pathvar_SHOX_type_8.csv')



LR_malegenital_SHOX <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SHOX, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_SHOX)
# exp(cbind(OR = coef(LR_malegenital_SHOX), confint(LR_malegenital_SHOX)))
pathvar_SHOX_type_9 <- tidy(LR_malegenital_SHOX, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SHOX_type_9, "pathvar_SHOX_type_9.csv", row.names=FALSE)
system('dx upload pathvar_SHOX_type_9.csv --path Emily-folder/results_2/pathvar_SHOX_type_9.csv')



LR_urinary_SHOX <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SHOX, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_SHOX)
# exp(cbind(OR = coef(LR_urinary_SHOX), confint(LR_urinary_SHOX)))
pathvar_SHOX_type_10 <- tidy(LR_urinary_SHOX, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SHOX_type_10, "pathvar_SHOX_type_10.csv", row.names=FALSE)
system('dx upload pathvar_SHOX_type_10.csv --path Emily-folder/results_2/pathvar_SHOX_type_10.csv')



LR_cns_SHOX <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SHOX, 
                   data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_SHOX)
# exp(cbind(OR = coef(LR_cns_SHOX), confint(LR_cns_SHOX)))
pathvar_SHOX_type_11 <- tidy(LR_cns_SHOX, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SHOX_type_11, "pathvar_SHOX_type_11.csv", row.names=FALSE)
system('dx upload pathvar_SHOX_type_11.csv --path Emily-folder/results_2/pathvar_SHOX_type_11.csv')



LR_endocrine_SHOX <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SHOX, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_SHOX)
# exp(cbind(OR = coef(LR_endocrine_SHOX), confint(LR_endocrine_SHOX)))
pathvar_SHOX_type_12 <- tidy(LR_endocrine_SHOX, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SHOX_type_12, "pathvar_SHOX_type_12.csv", row.names=FALSE)
system('dx upload pathvar_SHOX_type_12.csv --path Emily-folder/results_2/pathvar_SHOX_type_12.csv')



LR_lymphatic_SHOX <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SHOX, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_SHOX)
# exp(cbind(OR = coef(LR_lymphatic_SHOX), confint(LR_lymphatic_SHOX)))
pathvar_SHOX_type_13 <- tidy(LR_lymphatic_SHOX, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SHOX_type_13, "pathvar_SHOX_type_13.csv", row.names=FALSE)
system('dx upload pathvar_SHOX_type_13.csv --path Emily-folder/results_2/pathvar_SHOX_type_13.csv')



LR_secondary_SHOX <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SHOX, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_SHOX)
# exp(cbind(OR = coef(LR_secondary_SHOX), confint(LR_secondary_SHOX)))
pathvar_SHOX_type_14 <- tidy(LR_secondary_SHOX, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SHOX_type_14, "pathvar_SHOX_type_14.csv", row.names=FALSE)
system('dx upload pathvar_SHOX_type_14.csv --path Emily-folder/results_2/pathvar_SHOX_type_14.csv')



LR_other_SHOX <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SHOX, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_SHOX)
# exp(cbind(OR = coef(LR_other_SHOX), confint(LR_other_SHOX)))
pathvar_SHOX_type_15 <- tidy(LR_other_SHOX, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SHOX_type_15, "pathvar_SHOX_type_15.csv", row.names=FALSE)
system('dx upload pathvar_SHOX_type_15.csv --path Emily-folder/results_2/pathvar_SHOX_type_15.csv')












# TFgene Age -------------------------------------------------------------

# make sure the execute the code in the 'TFgenes' section

### linear regression
# "TFgene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ TFgene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_TFgene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_TFgene_age_1, "pathvar_TFgene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_age_1.csv --path Emily-folder/results_2/pathvar_TFgene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ ERCC2+CREB3L1+NKX3_2+HMGB3+RAD21+SMC1A+SMC3+NIPBL+HDAC8+PQBP1+FOXG1+HMGA2+LARP7+SHOX+SOX11+TAF13+TBX15+TBX3+SP7+TRPS1+ZNF335+GSC, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_TFgene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_TFgene_age_2, "pathvar_TFgene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_age_2.csv --path Emily-folder/results_2/pathvar_TFgene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "TFgene", 'ERCC2', 'CREB3L1', 'NKX3_2', 'HMGB3', 
                                     'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 
                                     'HDAC8', 'PQBP1', 'FOXG1', 'HMGA2', 
                                     'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 'SP7', 
                                     'TRPS1', 'ZNF335', 'GSC')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     29

with(data_setD, table(GSC))
# HMGB3 = 13
# SMC1A = 6
# SMC3 = 88
# HDAC8 = 13
# PQBP1 = 32
# FOXG1 = 30
# HMGA2 = 15
# SHOX = 39
# SOX11 = 83
# TAF13 = 86
# SP7 = 61



data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "TFgene", 'ERCC2', 'CREB3L1', 'NKX3_2', 
                                     'HMGB3', 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 'HDAC8', 'PQBP1', 'FOXG1', 
                                     'HMGA2', 'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 
                                     'SP7', 'TRPS1', 'ZNF335', 'GSC')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    34

with(data_setE, table(GSC))
# HMGB3 = 11
# SMC1A = 4
# SMC3 = 66
# HDAC8 = 9
# PQBP1 = 22
# FOXG1 = 22
# HMGA2 = 11
# SHOX = 33
# SOX11 = 71
# TAF13 = 63
# SP7 = 52
# GSC = 90




### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ TFgene, data = data_setD)
# summary(LM_15)
pathvar_TFgene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_TFgene_age_3, "pathvar_TFgene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_age_3.csv --path Emily-folder/results_2/pathvar_TFgene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TFgene, data = data_setD)
# summary(LM_16)
pathvar_TFgene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_TFgene_age_4, "pathvar_TFgene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_age_4.csv --path Emily-folder/results_2/pathvar_TFgene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ ERCC2+CREB3L1+NKX3_2+HMGB3+RAD21+SMC1A+SMC3+NIPBL+HDAC8+PQBP1+FOXG1+HMGA2+LARP7+SHOX+SOX11+TAF13+TBX15+TBX3+SP7+TRPS1+ZNF335+GSC, data = data_setD)
# summary(LM_17)
pathvar_TFgene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_TFgene_age_5, "pathvar_TFgene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_age_5.csv --path Emily-folder/results_2/pathvar_TFgene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC2+CREB3L1+NKX3_2+HMGB3+RAD21+SMC1A+SMC3+NIPBL+HDAC8+PQBP1+FOXG1+HMGA2+LARP7+SHOX+SOX11+TAF13+TBX15+TBX3+SP7+TRPS1+ZNF335+GSC, data = data_setD)
# summary(LM_18)
pathvar_TFgene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_TFgene_age_6, "pathvar_TFgene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_age_6.csv --path Emily-folder/results_2/pathvar_TFgene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ TFgene, data = data_setE)
# summary(LM_15)
pathvar_TFgene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_TFgene_age_7, "pathvar_TFgene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_age_7.csv --path Emily-folder/results_2/pathvar_TFgene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+TFgene, data = data_setE)
# summary(LM_16)
pathvar_TFgene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_TFgene_age_8, "pathvar_TFgene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_age_8.csv --path Emily-folder/results_2/pathvar_TFgene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ ERCC2+CREB3L1+NKX3_2+HMGB3+RAD21+SMC1A+SMC3+NIPBL+HDAC8+PQBP1+FOXG1+HMGA2+LARP7+SHOX+SOX11+TAF13+TBX15+TBX3+SP7+TRPS1+ZNF335+GSC, data = data_setE)
# summary(LM_17)
pathvar_TFgene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_TFgene_age_9, "pathvar_TFgene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_age_9.csv --path Emily-folder/results_2/pathvar_TFgene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+ERCC2+CREB3L1+NKX3_2+HMGB3+RAD21+SMC1A+SMC3+NIPBL+HDAC8+PQBP1+FOXG1+HMGA2+LARP7+SHOX+SOX11+TAF13+TBX15+TBX3+SP7+TRPS1+ZNF335+GSC, data = data_setE)
# summary(LM_18)
pathvar_TFgene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_TFgene_age_10, "pathvar_TFgene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_TFgene_age_10.csv --path Emily-folder/results_2/pathvar_TFgene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(df_cancer_diagnosis, table(TFgene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(TFgene, Ethnicity_White)))
# p-value = 3.657e-11
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5001555 0.6789628
# sample estimates:
#   odds ratio 
# 0.5813482

with(df_cancer_diagnosis, table(TFgene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(TFgene, Ever_smoked_dummy)))
# p-value = 0.4018

with(df_cancer_diagnosis, table(TFgene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(TFgene, Sex_Male)))
# p-value = 0.1276

with(df_cancer_diagnosis, table(TFgene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(TFgene, University_education)))
# p-value = 0.7335

with(df_cancer_diagnosis, table(TFgene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(TFgene, Annual_income)))
# p-value = 0.9335


### data_setD
with(data_setD, table(TFgene, Ethnicity_White))			
fisher.test(with(data_setD, table(TFgene, Ethnicity_White)))
# p-value = 2.863e-11
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.4970897 0.6753782
# sample estimates:
#   odds ratio 
# 0.5780249 

with(data_setD, table(TFgene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(TFgene, Ever_smoked_dummy)))
# p-value = 0.4576

with(data_setD, table(TFgene, Sex_Male))			
fisher.test(with(data_setD, table(TFgene, Sex_Male)))
# p-value = 0.1075


### data_setE
with(data_setE, table(TFgene, Ethnicity_White))			
fisher.test(with(data_setE, table(TFgene, Ethnicity_White)))
# p-value = 6.563e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5015529 0.7349218
# sample estimates:
#   odds ratio 
# 0.6048001

with(data_setE, table(TFgene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(TFgene, Ever_smoked_dummy)))
# p-value = 0.8138

with(data_setE, table(TFgene, Sex_Male))			
fisher.test(with(data_setE, table(TFgene, Sex_Male)))
# p-value = 0.138

with(data_setE, table(TFgene, University_education))			
fisher.test(with(data_setE, table(TFgene, University_education)))
# p-value = 0.8099

with(data_setE, table(TFgene, Annual_income))			
fisher.test(with(data_setE, table(TFgene, Annual_income)))
# p-value = 0.9492



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 

# make TFgene into factor
df_cancer_diagnosis$TFgene <- as.factor(df_cancer_diagnosis$TFgene)

# try a violin plot



# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ TFgene, data = df_cancer_diagnosis)
# 



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ TFgene, data = df_cancer_diagnosis)
# W = 233698434, p-value = 0.1861



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TFgene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TFgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ TFgene, data = df_cancer_diagnosis)
# W = 227920156, p-value = 0.6288



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TFgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TFgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ TFgene, data = df_cancer_diagnosis)
# W = 234323352, p-value = 0.02298



# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(TFgene), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$TFgene <- as.factor(data_setD$TFgene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ TFgene, data = data_setD)
# 



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ TFgene, data = data_setD)
# W = 229185550, p-value = 0.234




### data_setE

data_setE$TFgene <- as.factor(data_setE$TFgene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ TFgene, data = data_setE)
# W = 145502012, p-value = 0.9918



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = TFgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ TFgene, data = data_setE)
# W = 146533722, p-value = 0.4872



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = TFgene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = TFgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ TFgene, data = data_setE)
# W = 145185965, p-value = 0.8211



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = TFgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = TFgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ TFgene, data = data_setE)
# W = 148193174, p-value = 0.06759


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(TFgene), scales = "free_y") +
  theme_light()
# 












# SHOX extras -----------------------------------------------------------



# look at the demographics of the SHOX participants

### Fisher test on these datasets

### EntireCohort_pathvar
with(EntireCohort_pathvar, table(SHOX, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(SHOX, cancer_ROOC_dummy)))
# p-value = 0.07477

with(EntireCohort_pathvar, table(SHOX, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(SHOX, Ethnicity_White)))
# p-value = 6.776e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.1636250 0.4349493
# sample estimates:
#   odds ratio 
# 0.2612552 

with(EntireCohort_pathvar, table(SHOX, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(SHOX, Ever_smoked_dummy)))
# p-value = 0.5263

with(EntireCohort_pathvar, table(SHOX, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(SHOX, Sex_Male)))
# p-value = 0.5941

with(EntireCohort_pathvar, table(SHOX, University_education))			
fisher.test(with(EntireCohort_pathvar, table(SHOX, University_education)))
# p-value = 0.1788

with(EntireCohort_pathvar, table(SHOX, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(SHOX, Annual_income))	)
# p-value = 0.4425

with(EntireCohort_pathvar, table(SHOX, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(SHOX, Smoker_previous)))
# p-value = 0.5746

with(EntireCohort_pathvar, table(SHOX, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(SHOX, Smoker_current))	)
# p-value = 0.6625



### data_setA
with(data_setA, table(SHOX, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(SHOX, cancer_ROOC_dummy)))
# p-value = 0.0571

with(data_setA, table(SHOX, Ethnicity_White))			
fisher.test(with(data_setA, table(SHOX, Ethnicity_White)))
# p-value = 6.287e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.1628621 0.4328779
# sample estimates:
#   odds ratio 
# 0.2600142

with(data_setA, table(SHOX, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(SHOX, Ever_smoked_dummy)))
# p-value = 0.4665

with(data_setA, table(SHOX, Sex_Male))			
fisher.test(with(data_setA, table(SHOX, Sex_Male)))
# p-value = 0.7202



### data_setB
with(data_setB, table(SHOX, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(SHOX, cancer_ROOC_dummy)))
# p-value = 0.04963
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9726016 2.3203887
# sample estimates:
#   odds ratio 
# 1.516513 

with(data_setB, table(SHOX, Ethnicity_White))			
fisher.test(with(data_setB, table(SHOX, Ethnicity_White)))
# p-value = 1.535e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.1497136 0.4686896
# sample estimates:
#   odds ratio 
# 0.2567725

with(data_setB, table(SHOX, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(SHOX, Ever_smoked_dummy)))
# p-value = 0.5495

with(data_setB, table(SHOX, Sex_Male))			
fisher.test(with(data_setB, table(SHOX, Sex_Male)))
# p-value = 0.5581

with(data_setB, table(SHOX, University_education))			
fisher.test(with(data_setB, table(SHOX, University_education)))
# p-value = 0.08396

with(data_setB, table(SHOX, Annual_income))			
fisher.test(with(data_setB, table(SHOX, Annual_income)))
# p-value = 0.4401




### can plot the continuous variables: 

# make SHOX into factor
df_to_plot <- EntireCohort_pathvar
df_to_plot$SHOX <- as.factor(df_to_plot$SHOX)

# try a violin plot

# Age_at_recruitment  
df_to_plot %>% 
  ggplot(mapping = aes(x = SHOX, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = SHOX, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ SHOX, data = df_to_plot)
# W = 34055083, p-value = 0.005596

# BMI  
df_to_plot %>% 
  ggplot(mapping = aes(x = SHOX, y = BMI)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = SHOX, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ SHOX, data = df_to_plot)
# W = 28380274, p-value = 0.3845

# Standing_height 
df_to_plot %>% 
  ggplot(mapping = aes(x = SHOX, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = SHOX, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ SHOX, data = df_to_plot)
# W = 37843464, p-value = 1.021e-07

# make this a bar plot
# Weekly_exercise
df_to_plot %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(SHOX), scales = "free_y") +
  theme_light()
# 



### data_setA

# make SHOX into factor
data_setA$SHOX <- as.factor(data_setA$SHOX)


# jpeg(file="SHOX_Age_v_A.jpeg", width = 300, height = 300)
# 
# dev.off()
# system('dx upload SHOX_Age_v_A.jpeg --path Emily-folder/pathvar/plots/SHOX_Age_v_A.jpeg')



# Age_at_recruitment 

jpeg(file="SHOX_Age_v_A.jpeg", width = 300, height = 300)
data_setA %>% 
  ggplot(mapping = aes(x = SHOX, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
dev.off()
system('dx upload SHOX_Age_v_A.jpeg --path Emily-folder/pathvar/plots/SHOX_Age_v_A.jpeg')

jpeg(file="SHOX_Age_b_A.jpeg", width = 300, height = 300)
data_setA %>% 
  ggplot(mapping = aes(x = SHOX, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
dev.off()
system('dx upload SHOX_Age_b_A.jpeg --path Emily-folder/pathvar/plots/SHOX_Age_b_A.jpeg')

wilcox.test(Age_at_recruitment ~ SHOX, data = data_setA)
# W = 33220380, p-value = 0.006218



### data_setB

# make SHOX into factor
data_setB$SHOX <- as.factor(data_setB$SHOX)

# Age_at_recruitment

jpeg(file="SHOX_Age_v_B.jpeg", width = 300, height = 300)
data_setB %>% 
  ggplot(mapping = aes(x = SHOX, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
dev.off()
system('dx upload SHOX_Age_v_B.jpeg --path Emily-folder/pathvar/plots/SHOX_Age_v_B.jpeg')

jpeg(file="SHOX_Age_b_B.jpeg", width = 300, height = 300)
data_setB %>% 
  ggplot(mapping = aes(x = SHOX, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
dev.off()
system('dx upload SHOX_Age_b_B.jpeg --path Emily-folder/pathvar/plots/SHOX_Age_b_B.jpeg')

wilcox.test(Age_at_recruitment ~ SHOX, data = data_setB)
# W = 22957410, p-value = 0.008784



# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = SHOX, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = SHOX, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ SHOX, data = data_setB)
# W = 18410530, p-value = 0.1572



# Standing_height 

jpeg(file="SHOX_Height_v_B.jpeg", width = 300, height = 300)
data_setB %>% 
  ggplot(mapping = aes(x = SHOX, y = Standing_height)) +
  geom_violin() +
  theme_light()
dev.off()
system('dx upload SHOX_Height_v_B.jpeg --path Emily-folder/pathvar/plots/SHOX_Height_v_B.jpeg')

jpeg(file="SHOX_Height_b_B.jpeg", width = 300, height = 300)
data_setB %>% 
  ggplot(mapping = aes(x = SHOX, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
dev.off()
system('dx upload SHOX_Height_b_B.jpeg --path Emily-folder/pathvar/plots/SHOX_Height_b_B.jpeg')

wilcox.test(Standing_height ~ SHOX, data = data_setB)
# W = 25758372, p-value = 3.3e-07


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(SHOX), scales = "free_y") +
  theme_light()
# similar











## looking at proportions of variants in the SHOX group

# all vars
SHOX_df <- EntireCohort_pathvar[EntireCohort_pathvar$SHOX == 1,]
ALL <- SHOX_df$VAR


# all setA vars
data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "TFgene", 'ERCC2', 'CREB3L1', 'NKX3_2', 
                                     'HMGB3', 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 'HDAC8', 'PQBP1', 
                                     'FOXG1', 'HMGA2', 'LARP7', 'SHOX', 'SOX11', 'TAF13', 
                                     'TBX15', 'TBX3', 'SP7', 'TRPS1', 'ZNF335', 'GSC', 'VAR')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     29

SHOX_df_A <- data_setA[data_setA$SHOX == 1,]
setA <- SHOX_df_A$VAR


# all setB vars
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "TFgene", 'ERCC2', 'CREB3L1', 'NKX3_2', 'HMGB3', 
                                     'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 
                                     'HDAC8', 'PQBP1', 'FOXG1', 'HMGA2', 
                                     'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 
                                     'SP7', 'TRPS1', 'ZNF335', 'GSC', 'VAR')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     34

SHOX_df_B <- data_setB[data_setB$SHOX == 1,]
setB <- SHOX_df_B$VAR


# all setD vars
data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "TFgene", 'ERCC2', 'CREB3L1', 'NKX3_2', 'HMGB3', 
                                     'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 
                                     'HDAC8', 'PQBP1', 'FOXG1', 'HMGA2', 
                                     'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 'SP7', 
                                     'TRPS1', 'ZNF335', 'GSC', 'VAR')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     29

SHOX_df_D <- data_setD[data_setD$SHOX == 1,]
setD <- SHOX_df_D$VAR


# all setE vars
data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "TFgene", 'ERCC2', 'CREB3L1', 'NKX3_2', 
                                     'HMGB3', 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 'HDAC8', 'PQBP1', 'FOXG1', 
                                     'HMGA2', 'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 
                                     'SP7', 'TRPS1', 'ZNF335', 'GSC', 'VAR')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    34

SHOX_df_E <- data_setE[data_setE$SHOX == 1,]
setE <- SHOX_df_E$VAR


# save these as csv files to work with in R on my own computer
write.csv(ALL, "SHOX_VAR_all.csv", row.names=FALSE)
system('dx upload SHOX_VAR_all.csv --path Emily-folder/results_2/SHOX_VAR_all.csv')

write.csv(setA, "SHOX_VAR_setA.csv", row.names=FALSE)
system('dx upload SHOX_VAR_setA.csv --path Emily-folder/results_2/SHOX_VAR_setA.csv')

write.csv(setB, "SHOX_VAR_setB.csv", row.names=FALSE)
system('dx upload SHOX_VAR_setB.csv --path Emily-folder/results_2/SHOX_VAR_setB.csv')

write.csv(setD, "SHOX_VAR_setD.csv", row.names=FALSE)
system('dx upload SHOX_VAR_setD.csv --path Emily-folder/results_2/SHOX_VAR_setD.csv')

write.csv(setE, "SHOX_VAR_setE.csv", row.names=FALSE)
system('dx upload SHOX_VAR_setE.csv --path Emily-folder/results_2/SHOX_VAR_setE.csv')











# SHOX 23:644586:G:A_A ----------------------------------------------------------


# -	Look at the 23:644586:G:A_A variant (the most frequent SHOX var)

unique(EntireCohort_pathvar$SHOX)
# 0 1

SHOX_df <- EntireCohort_pathvar[EntireCohort_pathvar$SHOX == 1,]
# 127 obs. 

# separate the var names into 9 columns and then use these to make a new column
SHOX_df <- separate(data = SHOX_df, col = VAR, into = c("VARIANT1","VARIANT2","VARIANT3","VARIANT4","VARIANT5","VARIANT6","VARIANT7","VARIANT8","VARIANT9"), sep = ";")

# find participants with the 23:644586:G:A_A variant
VARIANT1_df <- SHOX_df[SHOX_df$VARIANT1 == "23:644586:G:A_A",]
VARIANT1_df$eid
# 5845855 1872011 2180873 1859775 3838654 5463619 4430562 1965179 2396447 4362882 3896349 5701127 5174604

VARIANT2_df <- SHOX_df[SHOX_df$VARIANT2 == "23:644586:G:A_A",]
unique(VARIANT2_df$eid)
# 3579348 2299615 3914238 4917673 1458025 4808136 1557470 2644184 1968661 1822258 2907843 4523948 2585201 4843947 4428031 2016775 2495310 2278427 3509272 3368190 5188109 4931690

VARIANT3_df <- SHOX_df[SHOX_df$VARIANT3 == "23:644586:G:A_A",]
unique(VARIANT3_df$eid)
# 4405402 1484444 4422499 5924248 5480175 3697854 2241624 2275418 4520301

VARIANT4_df <- SHOX_df[SHOX_df$VARIANT4 == "23:644586:G:A_A",]
unique(VARIANT4_df$eid)
# 4453640 2558972 1558812


# none in the rest of the VAR columns



pop_var_eid <- c(5845855, 1872011, 2180873, 1859775, 3838654, 5463619, 4430562, 1965179, 2396447, 4362882, 
                 3896349, 5701127, 5174604, 3579348, 2299615, 3914238, 4917673, 1458025, 4808136, 1557470, 
                 2644184, 1968661, 1822258, 2907843, 4523948, 2585201, 4843947, 4428031, 2016775, 2495310, 
                 2278427, 3509272, 3368190, 5188109, 4931690, 4405402, 1484444, 4422499, 5924248, 5480175, 
                 3697854, 2241624, 2275418, 4520301, 4453640, 2558972, 1558812)

# 47


# annotate the EntireCohort_pathvar df with this info
pop_var_df <- data.frame(eid = pop_var_eid,
                         pop_var = 1)

# merge this with the EntireCohort_pathvar df
EntireCohort_pathvar_pop_var <- full_join(EntireCohort_pathvar, pop_var_df, by = 'eid')

EntireCohort_pathvar_pop_var$pop_var[is.na(EntireCohort_pathvar_pop_var$pop_var)] <- 0

unique(EntireCohort_pathvar_pop_var$pop_var)
# 0 1







# check incidence

# EntireCohort
LR_1 <- glm(formula = cancer_ROOC_dummy ~ pop_var, data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_1)
# Coefficients:
#   Estimate Std. Error  z value Pr(>|z|)    
# (Intercept) -1.174742   0.003436 -341.884   <2e-16 ***
#   pop_var     -0.010882   0.344526   -0.032    0.975
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
# OR     2.5 %    97.5 %
#   (Intercept) 0.3088988 0.3068243 0.3109849
# pop_var     0.9891769 0.4799358 1.8786417


data_setA <- EntireCohort_pathvar_pop_var[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                             "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                             "TFgene", 'ERCC2', 'CREB3L1', 'NKX3_2', 
                                             'HMGB3', 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 'HDAC8', 'PQBP1', 
                                             'FOXG1', 'HMGA2', 'LARP7', 'SHOX', 'SOX11', 'TAF13', 
                                             'TBX15', 'TBX3', 'SP7', 'TRPS1', 'ZNF335', 'GSC', 'pop_var')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     30


# data set A
LR_1 <- glm(formula = cancer_ROOC_dummy ~ pop_var, data = data_setA, family = binomial)
# summary(LR_1)
# Coefficients:
#   Estimate Std. Error  z value Pr(>|z|)    
# (Intercept) -1.174877   0.003451 -340.482   <2e-16 ***
#   pop_var     -0.010747   0.344527   -0.031    0.975
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
# OR     2.5 %    97.5 %
#   (Intercept) 0.3088570 0.3067740 0.3109517
# pop_var     0.9893106 0.4800006 1.8788963



data_setB <- EntireCohort_pathvar_pop_var[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                             "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                             "BMI", "Standing_height", "University_education",
                                             "Annual_income", "Weekly_exercise",
                                             "TFgene", 'ERCC2', 'CREB3L1', 'NKX3_2', 'HMGB3', 
                                             'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 
                                             'HDAC8', 'PQBP1', 'FOXG1', 'HMGA2', 
                                             'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 
                                             'SP7', 'TRPS1', 'ZNF335', 'GSC', 'pop_var')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     35

# data set A
LR_1 <- glm(formula = cancer_ROOC_dummy ~ pop_var, data = data_setB, family = binomial)
# summary(LR_1)
# Coefficients:
#   Estimate Std. Error  z value Pr(>|z|)    
# (Intercept) -1.196447   0.003837 -311.830   <2e-16 ***
#   pop_var      0.002525   0.360995    0.007    0.994
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
# OR     2.5 %    97.5 %
#   (Intercept) 0.3022661 0.3000001 0.3045462
# pop_var     1.0025281 0.4684376 1.9601108





# look at height

# make SHOX into factor
df_to_plot <- EntireCohort_pathvar_pop_var
df_to_plot$pop_var <- as.factor(df_to_plot$pop_var)



# Standing_height 
df_to_plot %>% 
  ggplot(mapping = aes(x = pop_var, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = pop_var, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ pop_var, data = df_to_plot)
# W = 13990905, p-value = 0.001279











# SHOX 23:644455:C:T_T ----------------------------------------------------------


# -	Look at the 23:644455:C:T_T variant (the most frequent SHOX var)

unique(EntireCohort_pathvar$SHOX)
# 0 1

SHOX_df <- EntireCohort_pathvar[EntireCohort_pathvar$SHOX == 1,]
# 127 obs. 

# separate the var names into 9 columns and then use these to make a new column
SHOX_df <- separate(data = SHOX_df, col = VAR, into = c("VARIANT1","VARIANT2","VARIANT3","VARIANT4","VARIANT5","VARIANT6","VARIANT7","VARIANT8","VARIANT9"), sep = ";")

# find participants with the 23:644455:C:T_T variant
VARIANT1_df <- SHOX_df[SHOX_df$VARIANT1 == "23:644455:C:T_T",]
VARIANT1_df$eid
# 4351995 4402667 1771248

VARIANT2_df <- SHOX_df[SHOX_df$VARIANT2 == "23:644455:C:T_T",]
unique(VARIANT2_df$eid)
# 5619060 1468072 4581079 2668466 4443684 3842223 1950528 1572697

VARIANT3_df <- SHOX_df[SHOX_df$VARIANT3 == "23:644455:C:T_T",]
unique(VARIANT3_df$eid)
# 1379548 5760459 4144558 3010888 1415402 5735149 2897500

VARIANT4_df <- SHOX_df[SHOX_df$VARIANT4 == "23:644455:C:T_T",]
unique(VARIANT4_df$eid)
# 3797350 1115113 5068239 3012988

VARIANT5_df <- SHOX_df[SHOX_df$VARIANT5 == "23:644455:C:T_T",]
unique(VARIANT5_df$eid)
# 2715902

VARIANT6_df <- SHOX_df[SHOX_df$VARIANT6 == "23:644455:C:T_T",]
unique(VARIANT6_df$eid)
# 3035712


# none in the rest of the VAR columns



pop_var2_eid <- c(4351995, 4402667, 1771248, 5619060, 1468072, 4581079, 2668466, 4443684, 3842223, 
                  1950528, 1572697,1379548, 5760459, 4144558, 3010888, 1415402, 5735149, 2897500,
                  3797350, 1115113, 5068239, 3012988,2715902,3035712)

# 24


# annotate the EntireCohort_pathvar df with this info
pop_var2_df <- data.frame(eid = pop_var2_eid,
                          pop_var2 = 1)

# merge this with the EntireCohort_pathvar df
EntireCohort_pathvar_pop_var2 <- full_join(EntireCohort_pathvar, pop_var2_df, by = 'eid')

EntireCohort_pathvar_pop_var2$pop_var2[is.na(EntireCohort_pathvar_pop_var2$pop_var2)] <- 0

unique(EntireCohort_pathvar_pop_var2$pop_var2)
# 0 1







# check incidence

# EntireCohort
LR_1 <- glm(formula = cancer_ROOC_dummy ~ pop_var2, data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_1)
# Coefficients:
# Estimate Std. Error  z value Pr(>|z|)    
# (Intercept) -1.174782   0.003436 -341.901   <2e-16 ***
#   pop_var2     0.663956   0.421651    1.575    0.115
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
#                   OR     2.5 %    97.5 %
# (Intercept) 0.3088863 0.3068119 0.3109723
# pop_var2    1.9424623 0.8159087 4.3645101


data_setA <- EntireCohort_pathvar_pop_var2[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                              "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                              "TFgene", 'ERCC2', 'CREB3L1', 'NKX3_2', 
                                              'HMGB3', 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 'HDAC8', 'PQBP1', 
                                              'FOXG1', 'HMGA2', 'LARP7', 'SHOX', 'SOX11', 'TAF13', 
                                              'TBX15', 'TBX3', 'SP7', 'TRPS1', 'ZNF335', 'GSC', 'pop_var2')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     30


# data set A
LR_1 <- glm(formula = cancer_ROOC_dummy ~ pop_var2, data = data_setA, family = binomial)
# summary(LR_1)
# Coefficients:
# Estimate Std. Error  z value Pr(>|z|)    
# (Intercept) -1.174923   0.003451 -340.501   <2e-16 ***
#   pop_var2     0.807198   0.433643    1.861   0.0627 . 
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
#                   OR     2.5 %    97.5 %
# (Intercept) 0.3088427 0.3067598 0.3109372
# pop_var2    2.2416191 0.9249484 5.1958306



data_setB <- EntireCohort_pathvar_pop_var2[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                              "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                              "BMI", "Standing_height", "University_education",
                                              "Annual_income", "Weekly_exercise",
                                              "TFgene", 'ERCC2', 'CREB3L1', 'NKX3_2', 'HMGB3', 
                                              'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 
                                              'HDAC8', 'PQBP1', 'FOXG1', 'HMGA2', 
                                              'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 
                                              'SP7', 'TRPS1', 'ZNF335', 'GSC', 'pop_var2')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     35

# data set A
LR_1 <- glm(formula = cancer_ROOC_dummy ~ pop_var2, data = data_setB, family = binomial)
# summary(LR_1)
# Coefficients:
# Estimate Std. Error  z value Pr(>|z|)    
# (Intercept) -1.196477   0.003837 -311.846   <2e-16 ***
#   pop_var2     0.590342   0.507534    1.163    0.245
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
#                    OR     2.5 %    97.5 %
# (Intercept) 0.3022571 0.2999911 0.3045371
# pop_var2    1.8046046 0.6216251 4.7429070





# look at height

# make SHOX into factor
df_to_plot <- EntireCohort_pathvar_pop_var2
df_to_plot$pop_var2 <- as.factor(df_to_plot$pop_var2)



# Standing_height 
df_to_plot %>% 
  ggplot(mapping = aes(x = pop_var2, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = pop_var2, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ pop_var2, data = df_to_plot)
# W = 7120764, p-value = 0.02349








# all non-pop_var SHOX variants -------------------------------------------

# identify all SHOX carriers that aren't in pop_var_eid or pop_var2_eid

SHOX_df <- EntireCohort_pathvar[EntireCohort_pathvar$SHOX == 1,]
# 127 obs. 

non_pop_var_eid <- SHOX_df[ ! SHOX_df$eid %in% pop_var_eid, ]
# 80
non_pop_var2_eid <- non_pop_var_eid[ ! non_pop_var_eid$eid %in% pop_var2_eid, ]
# 56
non_pop_var2_eid <- non_pop_var2_eid$eid




# annotate the EntireCohort_pathvar df with this info
non_pop_var_df <- data.frame(eid = non_pop_var2_eid,
                             non_pop_var = 1)

# merge this with the EntireCohort_pathvar df
EntireCohort_pathvar_non_pop_var <- full_join(EntireCohort_pathvar, non_pop_var_df, by = 'eid')

EntireCohort_pathvar_non_pop_var$non_pop_var[is.na(EntireCohort_pathvar_non_pop_var$non_pop_var)] <- 0

unique(EntireCohort_pathvar_non_pop_var$non_pop_var)
# 0 1







# check incidence

# EntireCohort
LR_1 <- glm(formula = cancer_ROOC_dummy ~ non_pop_var, data = EntireCohort_pathvar_non_pop_var, family = binomial)
# summary(LR_1)
# Coefficients:
# Estimate Std. Error  z value Pr(>|z|)    
# (Intercept) -1.174811   0.003436 -341.895   <2e-16 ***
#   non_pop_var  0.508332   0.282259    1.801   0.0717 . 
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
#                    OR     2.5 %    97.5 %
# (Intercept) 0.3088774 0.3068029 0.3109634
# non_pop_var 1.6625158 0.9369972 2.8538873


data_setA <- EntireCohort_pathvar_non_pop_var[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                                 "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                                 "TFgene", 'ERCC2', 'CREB3L1', 'NKX3_2', 
                                                 'HMGB3', 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 'HDAC8', 'PQBP1', 
                                                 'FOXG1', 'HMGA2', 'LARP7', 'SHOX', 'SOX11', 'TAF13', 
                                                 'TBX15', 'TBX3', 'SP7', 'TRPS1', 'ZNF335', 'GSC', 'non_pop_var')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     30


# data set A
LR_1 <- glm(formula = cancer_ROOC_dummy ~ non_pop_var, data = data_setA, family = binomial)
# summary(LR_1)
# Coefficients:
# Estimate Std. Error  z value Pr(>|z|)    
# (Intercept) -1.174947   0.003451 -340.493   <2e-16 ***
#   non_pop_var  0.508468   0.282260    1.801   0.0716 .  
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
#                    OR     2.5 %   97.5 %
# (Intercept) 0.3088354 0.3067525 0.310930
# non_pop_var 1.6627416 0.9371241 2.854276



data_setB <- EntireCohort_pathvar_non_pop_var[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                                 "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                                 "BMI", "Standing_height", "University_education",
                                                 "Annual_income", "Weekly_exercise",
                                                 "TFgene", 'ERCC2', 'CREB3L1', 'NKX3_2', 'HMGB3', 
                                                 'RAD21', 'SMC1A', 'SMC3', 'NIPBL', 
                                                 'HDAC8', 'PQBP1', 'FOXG1', 'HMGA2', 
                                                 'LARP7', 'SHOX', 'SOX11', 'TAF13', 'TBX15', 'TBX3', 
                                                 'SP7', 'TRPS1', 'ZNF335', 'GSC', 'non_pop_var')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     35

# data set A
LR_1 <- glm(formula = cancer_ROOC_dummy ~ non_pop_var, data = data_setB, family = binomial)
# summary(LR_1)
# Coefficients:
# Estimate Std. Error  z value Pr(>|z|)    
# (Intercept) -1.196544   0.003837 -311.846   <2e-16 ***
#   non_pop_var  0.697552   0.307494    2.269   0.0233 *  
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
#                    OR     2.5 %   97.5 %
# (Intercept) 0.3022371 0.2999711 0.304517
# non_pop_var 2.0088300 1.0781675 3.632645





# look at height

# make SHOX into factor
df_to_plot <- EntireCohort_pathvar_non_pop_var
df_to_plot$non_pop_var <- as.factor(df_to_plot$non_pop_var)



# Standing_height 
df_to_plot %>% 
  ggplot(mapping = aes(x = non_pop_var, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = non_pop_var, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ non_pop_var, data = df_to_plot)
# W = 16736900, p-value = 0.0003412
















# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_TFgenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_TFgenes.R')
