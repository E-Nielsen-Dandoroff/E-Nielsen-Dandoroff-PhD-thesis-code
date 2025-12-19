# 09/07/2024

# EntireCohort_pathvar_CAgenes.R

# this script is for the Cancer associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)


# CAgene
# BRAF
# BUB1B
# CBL
# CREBBP
# EP300
# ERCC3
# FANCA
# KDM6A
# KRAS
# LIG4
# LZTR1
# NRAS
# POLE
# PTPN11
# RAF1
# RECQL3
# RECQL4
# RIT1
# RRAS2
# SBDS
# SMAD4
# SMARCAL1
# SOS1
# TALDO1
# TRIM37
# XRCC4


#'BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1','TRIM37','XRCC4'





# CAgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')

# go thr all the genes to check that they're present
grep("XRCC4", EntireCohort_pathvar$GENE)
# all seem to be present exccept for RAF1

# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$CAgene <- NA
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE1 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1',  'TRIM37','XRCC4')] <- 1
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE2 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1',  'TRIM37','XRCC4')] <- 1
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE3 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1',  'TRIM37','XRCC4')] <- 1
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE4 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1',  'TRIM37','XRCC4')] <- 1
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE5 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1',  'TRIM37','XRCC4')] <- 1
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE6 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1',  'TRIM37','XRCC4')] <- 1
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE7 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1',  'TRIM37','XRCC4')] <- 1
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE8 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1',  'TRIM37','XRCC4')] <- 1
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE9 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1',  'TRIM37','XRCC4')] <- 1


EntireCohort_pathvar$CAgene[is.na(EntireCohort_pathvar$CAgene)] <- 0
EntireCohort_pathvar$CAgene <- as.numeric(EntireCohort_pathvar$CAgene)
unique(EntireCohort_pathvar$CAgene)
# 0 1

table(EntireCohort_pathvar$CAgene)
#      0      1 
# 436600  33199 
# as of 18/9/25



# now to make a column for each of the CAgenes
### 'BRAF'
EntireCohort_pathvar$BRAF <- NA
EntireCohort_pathvar$BRAF[EntireCohort_pathvar$GENE1 %in% c("BRAF")] <- 1
EntireCohort_pathvar$BRAF[EntireCohort_pathvar$GENE2 %in% c("BRAF")] <- 1
EntireCohort_pathvar$BRAF[EntireCohort_pathvar$GENE3 %in% c("BRAF")] <- 1
EntireCohort_pathvar$BRAF[EntireCohort_pathvar$GENE4 %in% c("BRAF")] <- 1
EntireCohort_pathvar$BRAF[EntireCohort_pathvar$GENE5 %in% c("BRAF")] <- 1
EntireCohort_pathvar$BRAF[EntireCohort_pathvar$GENE6 %in% c("BRAF")] <- 1
EntireCohort_pathvar$BRAF[EntireCohort_pathvar$GENE7 %in% c("BRAF")] <- 1
EntireCohort_pathvar$BRAF[EntireCohort_pathvar$GENE8 %in% c("BRAF")] <- 1
EntireCohort_pathvar$BRAF[EntireCohort_pathvar$GENE9 %in% c("BRAF")] <- 1

EntireCohort_pathvar$BRAF[is.na(EntireCohort_pathvar$BRAF)] <- 0
EntireCohort_pathvar$BRAF <- as.numeric(EntireCohort_pathvar$BRAF)
unique(EntireCohort_pathvar$BRAF)
# 0 1

table(EntireCohort_pathvar$BRAF)
#      0      1 
# 469798      1





### 'BUB1B'
EntireCohort_pathvar$BUB1B <- NA
EntireCohort_pathvar$BUB1B[EntireCohort_pathvar$GENE1 %in% c("BUB1B")] <- 1
EntireCohort_pathvar$BUB1B[EntireCohort_pathvar$GENE2 %in% c("BUB1B")] <- 1
EntireCohort_pathvar$BUB1B[EntireCohort_pathvar$GENE3 %in% c("BUB1B")] <- 1
EntireCohort_pathvar$BUB1B[EntireCohort_pathvar$GENE4 %in% c("BUB1B")] <- 1
EntireCohort_pathvar$BUB1B[EntireCohort_pathvar$GENE5 %in% c("BUB1B")] <- 1
EntireCohort_pathvar$BUB1B[EntireCohort_pathvar$GENE6 %in% c("BUB1B")] <- 1
EntireCohort_pathvar$BUB1B[EntireCohort_pathvar$GENE7 %in% c("BUB1B")] <- 1
EntireCohort_pathvar$BUB1B[EntireCohort_pathvar$GENE8 %in% c("BUB1B")] <- 1
EntireCohort_pathvar$BUB1B[EntireCohort_pathvar$GENE9 %in% c("BUB1B")] <- 1

EntireCohort_pathvar$BUB1B[is.na(EntireCohort_pathvar$BUB1B)] <- 0
EntireCohort_pathvar$BUB1B <- as.numeric(EntireCohort_pathvar$BUB1B)
unique(EntireCohort_pathvar$BUB1B)
# 0 1

table(EntireCohort_pathvar$BUB1B)
#      0      1 
# 468589   1210




### 'CBL'
EntireCohort_pathvar$CBL <- NA
EntireCohort_pathvar$CBL[EntireCohort_pathvar$GENE1 %in% c("CBL")] <- 1
EntireCohort_pathvar$CBL[EntireCohort_pathvar$GENE2 %in% c("CBL")] <- 1
EntireCohort_pathvar$CBL[EntireCohort_pathvar$GENE3 %in% c("CBL")] <- 1
EntireCohort_pathvar$CBL[EntireCohort_pathvar$GENE4 %in% c("CBL")] <- 1
EntireCohort_pathvar$CBL[EntireCohort_pathvar$GENE5 %in% c("CBL")] <- 1
EntireCohort_pathvar$CBL[EntireCohort_pathvar$GENE6 %in% c("CBL")] <- 1
EntireCohort_pathvar$CBL[EntireCohort_pathvar$GENE7 %in% c("CBL")] <- 1
EntireCohort_pathvar$CBL[EntireCohort_pathvar$GENE8 %in% c("CBL")] <- 1
EntireCohort_pathvar$CBL[EntireCohort_pathvar$GENE9 %in% c("CBL")] <- 1

EntireCohort_pathvar$CBL[is.na(EntireCohort_pathvar$CBL)] <- 0
EntireCohort_pathvar$CBL <- as.numeric(EntireCohort_pathvar$CBL)
unique(EntireCohort_pathvar$CBL)
# 0 1

table(EntireCohort_pathvar$CBL)
#      0      1 
# 469780     19




### 'CREBBP'
EntireCohort_pathvar$CREBBP <- NA
EntireCohort_pathvar$CREBBP[EntireCohort_pathvar$GENE1 %in% c("CREBBP")] <- 1
EntireCohort_pathvar$CREBBP[EntireCohort_pathvar$GENE2 %in% c("CREBBP")] <- 1
EntireCohort_pathvar$CREBBP[EntireCohort_pathvar$GENE3 %in% c("CREBBP")] <- 1
EntireCohort_pathvar$CREBBP[EntireCohort_pathvar$GENE4 %in% c("CREBBP")] <- 1
EntireCohort_pathvar$CREBBP[EntireCohort_pathvar$GENE5 %in% c("CREBBP")] <- 1
EntireCohort_pathvar$CREBBP[EntireCohort_pathvar$GENE6 %in% c("CREBBP")] <- 1
EntireCohort_pathvar$CREBBP[EntireCohort_pathvar$GENE7 %in% c("CREBBP")] <- 1
EntireCohort_pathvar$CREBBP[EntireCohort_pathvar$GENE8 %in% c("CREBBP")] <- 1
EntireCohort_pathvar$CREBBP[EntireCohort_pathvar$GENE9 %in% c("CREBBP")] <- 1

EntireCohort_pathvar$CREBBP[is.na(EntireCohort_pathvar$CREBBP)] <- 0
EntireCohort_pathvar$CREBBP <- as.numeric(EntireCohort_pathvar$CREBBP)
unique(EntireCohort_pathvar$CREBBP)
# 0 1

table(EntireCohort_pathvar$CREBBP)
#      0      1 
# 467486   2313




### 'EP300' 
EntireCohort_pathvar$EP300 <- NA
EntireCohort_pathvar$EP300[EntireCohort_pathvar$GENE1 %in% c("EP300")] <- 1
EntireCohort_pathvar$EP300[EntireCohort_pathvar$GENE2 %in% c("EP300")] <- 1
EntireCohort_pathvar$EP300[EntireCohort_pathvar$GENE3 %in% c("EP300")] <- 1
EntireCohort_pathvar$EP300[EntireCohort_pathvar$GENE4 %in% c("EP300")] <- 1
EntireCohort_pathvar$EP300[EntireCohort_pathvar$GENE5 %in% c("EP300")] <- 1
EntireCohort_pathvar$EP300[EntireCohort_pathvar$GENE6 %in% c("EP300")] <- 1
EntireCohort_pathvar$EP300[EntireCohort_pathvar$GENE7 %in% c("EP300")] <- 1
EntireCohort_pathvar$EP300[EntireCohort_pathvar$GENE8 %in% c("EP300")] <- 1
EntireCohort_pathvar$EP300[EntireCohort_pathvar$GENE9 %in% c("EP300")] <- 1

EntireCohort_pathvar$EP300[is.na(EntireCohort_pathvar$EP300)] <- 0
EntireCohort_pathvar$EP300 <- as.numeric(EntireCohort_pathvar$EP300)
unique(EntireCohort_pathvar$EP300)
# 0 1

table(EntireCohort_pathvar$EP300)
#      0      1 
# 467607   2192




### 'ERCC3' 
EntireCohort_pathvar$ERCC3 <- NA
EntireCohort_pathvar$ERCC3[EntireCohort_pathvar$GENE1 %in% c("ERCC3")] <- 1
EntireCohort_pathvar$ERCC3[EntireCohort_pathvar$GENE2 %in% c("ERCC3")] <- 1
EntireCohort_pathvar$ERCC3[EntireCohort_pathvar$GENE3 %in% c("ERCC3")] <- 1
EntireCohort_pathvar$ERCC3[EntireCohort_pathvar$GENE4 %in% c("ERCC3")] <- 1
EntireCohort_pathvar$ERCC3[EntireCohort_pathvar$GENE5 %in% c("ERCC3")] <- 1
EntireCohort_pathvar$ERCC3[EntireCohort_pathvar$GENE6 %in% c("ERCC3")] <- 1
EntireCohort_pathvar$ERCC3[EntireCohort_pathvar$GENE7 %in% c("ERCC3")] <- 1
EntireCohort_pathvar$ERCC3[EntireCohort_pathvar$GENE8 %in% c("ERCC3")] <- 1
EntireCohort_pathvar$ERCC3[EntireCohort_pathvar$GENE9 %in% c("ERCC3")] <- 1

EntireCohort_pathvar$ERCC3[is.na(EntireCohort_pathvar$ERCC3)] <- 0
EntireCohort_pathvar$ERCC3 <- as.numeric(EntireCohort_pathvar$ERCC3)
unique(EntireCohort_pathvar$ERCC3)
# 0 1

table(EntireCohort_pathvar$ERCC3)
#      0      1 
# 465544   4255




### 'FANCA' 
EntireCohort_pathvar$FANCA <- NA
EntireCohort_pathvar$FANCA[EntireCohort_pathvar$GENE1 %in% c("FANCA")] <- 1
EntireCohort_pathvar$FANCA[EntireCohort_pathvar$GENE2 %in% c("FANCA")] <- 1
EntireCohort_pathvar$FANCA[EntireCohort_pathvar$GENE3 %in% c("FANCA")] <- 1
EntireCohort_pathvar$FANCA[EntireCohort_pathvar$GENE4 %in% c("FANCA")] <- 1
EntireCohort_pathvar$FANCA[EntireCohort_pathvar$GENE5 %in% c("FANCA")] <- 1
EntireCohort_pathvar$FANCA[EntireCohort_pathvar$GENE6 %in% c("FANCA")] <- 1
EntireCohort_pathvar$FANCA[EntireCohort_pathvar$GENE7 %in% c("FANCA")] <- 1
EntireCohort_pathvar$FANCA[EntireCohort_pathvar$GENE8 %in% c("FANCA")] <- 1
EntireCohort_pathvar$FANCA[EntireCohort_pathvar$GENE9 %in% c("FANCA")] <- 1

EntireCohort_pathvar$FANCA[is.na(EntireCohort_pathvar$FANCA)] <- 0
EntireCohort_pathvar$FANCA <- as.numeric(EntireCohort_pathvar$FANCA)
unique(EntireCohort_pathvar$FANCA)
# 0 1

table(EntireCohort_pathvar$FANCA)
#      0      1 
# 467829   1970




### 'KDM6A' 
EntireCohort_pathvar$KDM6A <- NA
EntireCohort_pathvar$KDM6A[EntireCohort_pathvar$GENE1 %in% c("KDM6A")] <- 1
EntireCohort_pathvar$KDM6A[EntireCohort_pathvar$GENE2 %in% c("KDM6A")] <- 1
EntireCohort_pathvar$KDM6A[EntireCohort_pathvar$GENE3 %in% c("KDM6A")] <- 1
EntireCohort_pathvar$KDM6A[EntireCohort_pathvar$GENE4 %in% c("KDM6A")] <- 1
EntireCohort_pathvar$KDM6A[EntireCohort_pathvar$GENE5 %in% c("KDM6A")] <- 1
EntireCohort_pathvar$KDM6A[EntireCohort_pathvar$GENE6 %in% c("KDM6A")] <- 1
EntireCohort_pathvar$KDM6A[EntireCohort_pathvar$GENE7 %in% c("KDM6A")] <- 1
EntireCohort_pathvar$KDM6A[EntireCohort_pathvar$GENE8 %in% c("KDM6A")] <- 1
EntireCohort_pathvar$KDM6A[EntireCohort_pathvar$GENE9 %in% c("KDM6A")] <- 1

EntireCohort_pathvar$KDM6A[is.na(EntireCohort_pathvar$KDM6A)] <- 0
EntireCohort_pathvar$KDM6A <- as.numeric(EntireCohort_pathvar$KDM6A)
unique(EntireCohort_pathvar$KDM6A)
# 0 1

table(EntireCohort_pathvar$KDM6A)
#      0      1 
# 469141    658




### 'KRAS' 
EntireCohort_pathvar$KRAS <- NA
EntireCohort_pathvar$KRAS[EntireCohort_pathvar$GENE1 %in% c("KRAS")] <- 1
EntireCohort_pathvar$KRAS[EntireCohort_pathvar$GENE2 %in% c("KRAS")] <- 1
EntireCohort_pathvar$KRAS[EntireCohort_pathvar$GENE3 %in% c("KRAS")] <- 1
EntireCohort_pathvar$KRAS[EntireCohort_pathvar$GENE4 %in% c("KRAS")] <- 1
EntireCohort_pathvar$KRAS[EntireCohort_pathvar$GENE5 %in% c("KRAS")] <- 1
EntireCohort_pathvar$KRAS[EntireCohort_pathvar$GENE6 %in% c("KRAS")] <- 1
EntireCohort_pathvar$KRAS[EntireCohort_pathvar$GENE7 %in% c("KRAS")] <- 1
EntireCohort_pathvar$KRAS[EntireCohort_pathvar$GENE8 %in% c("KRAS")] <- 1
EntireCohort_pathvar$KRAS[EntireCohort_pathvar$GENE9 %in% c("KRAS")] <- 1

EntireCohort_pathvar$KRAS[is.na(EntireCohort_pathvar$KRAS)] <- 0
EntireCohort_pathvar$KRAS <- as.numeric(EntireCohort_pathvar$KRAS)
unique(EntireCohort_pathvar$KRAS)
# 0 1

table(EntireCohort_pathvar$KRAS)
#      0      1 
# 469794      5





### 'LIG4' 
EntireCohort_pathvar$LIG4 <- NA
EntireCohort_pathvar$LIG4[EntireCohort_pathvar$GENE1 %in% c("LIG4")] <- 1
EntireCohort_pathvar$LIG4[EntireCohort_pathvar$GENE2 %in% c("LIG4")] <- 1
EntireCohort_pathvar$LIG4[EntireCohort_pathvar$GENE3 %in% c("LIG4")] <- 1
EntireCohort_pathvar$LIG4[EntireCohort_pathvar$GENE4 %in% c("LIG4")] <- 1
EntireCohort_pathvar$LIG4[EntireCohort_pathvar$GENE5 %in% c("LIG4")] <- 1
EntireCohort_pathvar$LIG4[EntireCohort_pathvar$GENE6 %in% c("LIG4")] <- 1
EntireCohort_pathvar$LIG4[EntireCohort_pathvar$GENE7 %in% c("LIG4")] <- 1
EntireCohort_pathvar$LIG4[EntireCohort_pathvar$GENE8 %in% c("LIG4")] <- 1
EntireCohort_pathvar$LIG4[EntireCohort_pathvar$GENE9 %in% c("LIG4")] <- 1

EntireCohort_pathvar$LIG4[is.na(EntireCohort_pathvar$LIG4)] <- 0
EntireCohort_pathvar$LIG4 <- as.numeric(EntireCohort_pathvar$LIG4)
unique(EntireCohort_pathvar$LIG4)
# 0 1

table(EntireCohort_pathvar$LIG4)
#      0      1 
# 468426   1373




### 'LZTR1' 
EntireCohort_pathvar$LZTR1 <- NA
EntireCohort_pathvar$LZTR1[EntireCohort_pathvar$GENE1 %in% c("LZTR1")] <- 1
EntireCohort_pathvar$LZTR1[EntireCohort_pathvar$GENE2 %in% c("LZTR1")] <- 1
EntireCohort_pathvar$LZTR1[EntireCohort_pathvar$GENE3 %in% c("LZTR1")] <- 1
EntireCohort_pathvar$LZTR1[EntireCohort_pathvar$GENE4 %in% c("LZTR1")] <- 1
EntireCohort_pathvar$LZTR1[EntireCohort_pathvar$GENE5 %in% c("LZTR1")] <- 1
EntireCohort_pathvar$LZTR1[EntireCohort_pathvar$GENE6 %in% c("LZTR1")] <- 1
EntireCohort_pathvar$LZTR1[EntireCohort_pathvar$GENE7 %in% c("LZTR1")] <- 1
EntireCohort_pathvar$LZTR1[EntireCohort_pathvar$GENE8 %in% c("LZTR1")] <- 1
EntireCohort_pathvar$LZTR1[EntireCohort_pathvar$GENE9 %in% c("LZTR1")] <- 1

EntireCohort_pathvar$LZTR1[is.na(EntireCohort_pathvar$LZTR1)] <- 0
EntireCohort_pathvar$LZTR1 <- as.numeric(EntireCohort_pathvar$LZTR1)
unique(EntireCohort_pathvar$LZTR1)
# 0 1

table(EntireCohort_pathvar$LZTR1)
#      0      1 
# 467150   2649





### 'NRAS' 
EntireCohort_pathvar$NRAS <- NA
EntireCohort_pathvar$NRAS[EntireCohort_pathvar$GENE1 %in% c("NRAS")] <- 1
EntireCohort_pathvar$NRAS[EntireCohort_pathvar$GENE2 %in% c("NRAS")] <- 1
EntireCohort_pathvar$NRAS[EntireCohort_pathvar$GENE3 %in% c("NRAS")] <- 1
EntireCohort_pathvar$NRAS[EntireCohort_pathvar$GENE4 %in% c("NRAS")] <- 1
EntireCohort_pathvar$NRAS[EntireCohort_pathvar$GENE5 %in% c("NRAS")] <- 1
EntireCohort_pathvar$NRAS[EntireCohort_pathvar$GENE6 %in% c("NRAS")] <- 1
EntireCohort_pathvar$NRAS[EntireCohort_pathvar$GENE7 %in% c("NRAS")] <- 1
EntireCohort_pathvar$NRAS[EntireCohort_pathvar$GENE8 %in% c("NRAS")] <- 1
EntireCohort_pathvar$NRAS[EntireCohort_pathvar$GENE9 %in% c("NRAS")] <- 1

EntireCohort_pathvar$NRAS[is.na(EntireCohort_pathvar$NRAS)] <- 0
EntireCohort_pathvar$NRAS <- as.numeric(EntireCohort_pathvar$NRAS)
unique(EntireCohort_pathvar$NRAS)
# 0 1

table(EntireCohort_pathvar$NRAS)
#      0      1 
# 469786     13





### 'POLE'
EntireCohort_pathvar$POLE <- NA
EntireCohort_pathvar$POLE[EntireCohort_pathvar$GENE1 %in% c("POLE")] <- 1
EntireCohort_pathvar$POLE[EntireCohort_pathvar$GENE2 %in% c("POLE")] <- 1
EntireCohort_pathvar$POLE[EntireCohort_pathvar$GENE3 %in% c("POLE")] <- 1
EntireCohort_pathvar$POLE[EntireCohort_pathvar$GENE4 %in% c("POLE")] <- 1
EntireCohort_pathvar$POLE[EntireCohort_pathvar$GENE5 %in% c("POLE")] <- 1
EntireCohort_pathvar$POLE[EntireCohort_pathvar$GENE6 %in% c("POLE")] <- 1
EntireCohort_pathvar$POLE[EntireCohort_pathvar$GENE7 %in% c("POLE")] <- 1
EntireCohort_pathvar$POLE[EntireCohort_pathvar$GENE8 %in% c("POLE")] <- 1
EntireCohort_pathvar$POLE[EntireCohort_pathvar$GENE9 %in% c("POLE")] <- 1

EntireCohort_pathvar$POLE[is.na(EntireCohort_pathvar$POLE)] <- 0
EntireCohort_pathvar$POLE <- as.numeric(EntireCohort_pathvar$POLE)
unique(EntireCohort_pathvar$POLE)
# 0 1

table(EntireCohort_pathvar$POLE)
#      0      1 
# 465246   4553




### 'PTPN11'
EntireCohort_pathvar$PTPN11 <- NA
EntireCohort_pathvar$PTPN11[EntireCohort_pathvar$GENE1 %in% c("PTPN11")] <- 1
EntireCohort_pathvar$PTPN11[EntireCohort_pathvar$GENE2 %in% c("PTPN11")] <- 1
EntireCohort_pathvar$PTPN11[EntireCohort_pathvar$GENE3 %in% c("PTPN11")] <- 1
EntireCohort_pathvar$PTPN11[EntireCohort_pathvar$GENE4 %in% c("PTPN11")] <- 1
EntireCohort_pathvar$PTPN11[EntireCohort_pathvar$GENE5 %in% c("PTPN11")] <- 1
EntireCohort_pathvar$PTPN11[EntireCohort_pathvar$GENE6 %in% c("PTPN11")] <- 1
EntireCohort_pathvar$PTPN11[EntireCohort_pathvar$GENE7 %in% c("PTPN11")] <- 1
EntireCohort_pathvar$PTPN11[EntireCohort_pathvar$GENE8 %in% c("PTPN11")] <- 1
EntireCohort_pathvar$PTPN11[EntireCohort_pathvar$GENE9 %in% c("PTPN11")] <- 1

EntireCohort_pathvar$PTPN11[is.na(EntireCohort_pathvar$PTPN11)] <- 0
EntireCohort_pathvar$PTPN11 <- as.numeric(EntireCohort_pathvar$PTPN11)
unique(EntireCohort_pathvar$PTPN11)
# 0 1

table(EntireCohort_pathvar$PTPN11)
#      0      1 
# 469714     85




### 'RAF1'
EntireCohort_pathvar$RAF1 <- NA
EntireCohort_pathvar$RAF1[EntireCohort_pathvar$GENE1 %in% c("RAF1")] <- 1
EntireCohort_pathvar$RAF1[EntireCohort_pathvar$GENE2 %in% c("RAF1")] <- 1
EntireCohort_pathvar$RAF1[EntireCohort_pathvar$GENE3 %in% c("RAF1")] <- 1
EntireCohort_pathvar$RAF1[EntireCohort_pathvar$GENE4 %in% c("RAF1")] <- 1
EntireCohort_pathvar$RAF1[EntireCohort_pathvar$GENE5 %in% c("RAF1")] <- 1
EntireCohort_pathvar$RAF1[EntireCohort_pathvar$GENE6 %in% c("RAF1")] <- 1
EntireCohort_pathvar$RAF1[EntireCohort_pathvar$GENE7 %in% c("RAF1")] <- 1
EntireCohort_pathvar$RAF1[EntireCohort_pathvar$GENE8 %in% c("RAF1")] <- 1
EntireCohort_pathvar$RAF1[EntireCohort_pathvar$GENE9 %in% c("RAF1")] <- 1

EntireCohort_pathvar$RAF1[is.na(EntireCohort_pathvar$RAF1)] <- 0
EntireCohort_pathvar$RAF1 <- as.numeric(EntireCohort_pathvar$RAF1)
unique(EntireCohort_pathvar$RAF1)
# 0

table(EntireCohort_pathvar$RAF1)
# 0 
# 469799




### 'RECQL3'
EntireCohort_pathvar$RECQL3 <- NA
EntireCohort_pathvar$RECQL3[EntireCohort_pathvar$GENE1 %in% c("RECQL3")] <- 1
EntireCohort_pathvar$RECQL3[EntireCohort_pathvar$GENE2 %in% c("RECQL3")] <- 1
EntireCohort_pathvar$RECQL3[EntireCohort_pathvar$GENE3 %in% c("RECQL3")] <- 1
EntireCohort_pathvar$RECQL3[EntireCohort_pathvar$GENE4 %in% c("RECQL3")] <- 1
EntireCohort_pathvar$RECQL3[EntireCohort_pathvar$GENE5 %in% c("RECQL3")] <- 1
EntireCohort_pathvar$RECQL3[EntireCohort_pathvar$GENE6 %in% c("RECQL3")] <- 1
EntireCohort_pathvar$RECQL3[EntireCohort_pathvar$GENE7 %in% c("RECQL3")] <- 1
EntireCohort_pathvar$RECQL3[EntireCohort_pathvar$GENE8 %in% c("RECQL3")] <- 1
EntireCohort_pathvar$RECQL3[EntireCohort_pathvar$GENE9 %in% c("RECQL3")] <- 1

EntireCohort_pathvar$RECQL3[is.na(EntireCohort_pathvar$RECQL3)] <- 0
EntireCohort_pathvar$RECQL3 <- as.numeric(EntireCohort_pathvar$RECQL3)
unique(EntireCohort_pathvar$RECQL3)
# 0 1

table(EntireCohort_pathvar$RECQL3)
#      0      1 
# 468224   1575




### 'RECQL4'
EntireCohort_pathvar$RECQL4 <- NA
EntireCohort_pathvar$RECQL4[EntireCohort_pathvar$GENE1 %in% c("RECQL4")] <- 1
EntireCohort_pathvar$RECQL4[EntireCohort_pathvar$GENE2 %in% c("RECQL4")] <- 1
EntireCohort_pathvar$RECQL4[EntireCohort_pathvar$GENE3 %in% c("RECQL4")] <- 1
EntireCohort_pathvar$RECQL4[EntireCohort_pathvar$GENE4 %in% c("RECQL4")] <- 1
EntireCohort_pathvar$RECQL4[EntireCohort_pathvar$GENE5 %in% c("RECQL4")] <- 1
EntireCohort_pathvar$RECQL4[EntireCohort_pathvar$GENE6 %in% c("RECQL4")] <- 1
EntireCohort_pathvar$RECQL4[EntireCohort_pathvar$GENE7 %in% c("RECQL4")] <- 1
EntireCohort_pathvar$RECQL4[EntireCohort_pathvar$GENE8 %in% c("RECQL4")] <- 1
EntireCohort_pathvar$RECQL4[EntireCohort_pathvar$GENE9 %in% c("RECQL4")] <- 1

EntireCohort_pathvar$RECQL4[is.na(EntireCohort_pathvar$RECQL4)] <- 0
EntireCohort_pathvar$RECQL4 <- as.numeric(EntireCohort_pathvar$RECQL4)
unique(EntireCohort_pathvar$RECQL4)
# 0 1

table(EntireCohort_pathvar$RECQL4)
#      0      1 
# 467240   2559





### 'RIT1'
EntireCohort_pathvar$RIT1 <- NA
EntireCohort_pathvar$RIT1[EntireCohort_pathvar$GENE1 %in% c("RIT1")] <- 1
EntireCohort_pathvar$RIT1[EntireCohort_pathvar$GENE2 %in% c("RIT1")] <- 1
EntireCohort_pathvar$RIT1[EntireCohort_pathvar$GENE3 %in% c("RIT1")] <- 1
EntireCohort_pathvar$RIT1[EntireCohort_pathvar$GENE4 %in% c("RIT1")] <- 1
EntireCohort_pathvar$RIT1[EntireCohort_pathvar$GENE5 %in% c("RIT1")] <- 1
EntireCohort_pathvar$RIT1[EntireCohort_pathvar$GENE6 %in% c("RIT1")] <- 1
EntireCohort_pathvar$RIT1[EntireCohort_pathvar$GENE7 %in% c("RIT1")] <- 1
EntireCohort_pathvar$RIT1[EntireCohort_pathvar$GENE8 %in% c("RIT1")] <- 1
EntireCohort_pathvar$RIT1[EntireCohort_pathvar$GENE9 %in% c("RIT1")] <- 1

EntireCohort_pathvar$RIT1[is.na(EntireCohort_pathvar$RIT1)] <- 0
EntireCohort_pathvar$RIT1 <- as.numeric(EntireCohort_pathvar$RIT1)
unique(EntireCohort_pathvar$RIT1)
# 0 1

table(EntireCohort_pathvar$RIT1)
#      0      1 
# 469792      7




### 'RRAS2'
EntireCohort_pathvar$RRAS2 <- NA
EntireCohort_pathvar$RRAS2[EntireCohort_pathvar$GENE1 %in% c("RRAS2")] <- 1
EntireCohort_pathvar$RRAS2[EntireCohort_pathvar$GENE2 %in% c("RRAS2")] <- 1
EntireCohort_pathvar$RRAS2[EntireCohort_pathvar$GENE3 %in% c("RRAS2")] <- 1
EntireCohort_pathvar$RRAS2[EntireCohort_pathvar$GENE4 %in% c("RRAS2")] <- 1
EntireCohort_pathvar$RRAS2[EntireCohort_pathvar$GENE5 %in% c("RRAS2")] <- 1
EntireCohort_pathvar$RRAS2[EntireCohort_pathvar$GENE6 %in% c("RRAS2")] <- 1
EntireCohort_pathvar$RRAS2[EntireCohort_pathvar$GENE7 %in% c("RRAS2")] <- 1
EntireCohort_pathvar$RRAS2[EntireCohort_pathvar$GENE8 %in% c("RRAS2")] <- 1
EntireCohort_pathvar$RRAS2[EntireCohort_pathvar$GENE9 %in% c("RRAS2")] <- 1

EntireCohort_pathvar$RRAS2[is.na(EntireCohort_pathvar$RRAS2)] <- 0
EntireCohort_pathvar$RRAS2 <- as.numeric(EntireCohort_pathvar$RRAS2)
unique(EntireCohort_pathvar$RRAS2)
# 0 1

table(EntireCohort_pathvar$RRAS2)
#      0      1 
# 469798      1




### 'SBDS'
EntireCohort_pathvar$SBDS <- NA
EntireCohort_pathvar$SBDS[EntireCohort_pathvar$GENE1 %in% c("SBDS")] <- 1
EntireCohort_pathvar$SBDS[EntireCohort_pathvar$GENE2 %in% c("SBDS")] <- 1
EntireCohort_pathvar$SBDS[EntireCohort_pathvar$GENE3 %in% c("SBDS")] <- 1
EntireCohort_pathvar$SBDS[EntireCohort_pathvar$GENE4 %in% c("SBDS")] <- 1
EntireCohort_pathvar$SBDS[EntireCohort_pathvar$GENE5 %in% c("SBDS")] <- 1
EntireCohort_pathvar$SBDS[EntireCohort_pathvar$GENE6 %in% c("SBDS")] <- 1
EntireCohort_pathvar$SBDS[EntireCohort_pathvar$GENE7 %in% c("SBDS")] <- 1
EntireCohort_pathvar$SBDS[EntireCohort_pathvar$GENE8 %in% c("SBDS")] <- 1
EntireCohort_pathvar$SBDS[EntireCohort_pathvar$GENE9 %in% c("SBDS")] <- 1

EntireCohort_pathvar$SBDS[is.na(EntireCohort_pathvar$SBDS)] <- 0
EntireCohort_pathvar$SBDS <- as.numeric(EntireCohort_pathvar$SBDS)
unique(EntireCohort_pathvar$SBDS)
# 0 1

table(EntireCohort_pathvar$SBDS)
#      0      1 
# 465479   4320





### 'SMAD4'   
EntireCohort_pathvar$SMAD4 <- NA
EntireCohort_pathvar$SMAD4[EntireCohort_pathvar$GENE1 %in% c("SMAD4")] <- 1
EntireCohort_pathvar$SMAD4[EntireCohort_pathvar$GENE2 %in% c("SMAD4")] <- 1
EntireCohort_pathvar$SMAD4[EntireCohort_pathvar$GENE3 %in% c("SMAD4")] <- 1
EntireCohort_pathvar$SMAD4[EntireCohort_pathvar$GENE4 %in% c("SMAD4")] <- 1
EntireCohort_pathvar$SMAD4[EntireCohort_pathvar$GENE5 %in% c("SMAD4")] <- 1
EntireCohort_pathvar$SMAD4[EntireCohort_pathvar$GENE6 %in% c("SMAD4")] <- 1
EntireCohort_pathvar$SMAD4[EntireCohort_pathvar$GENE7 %in% c("SMAD4")] <- 1
EntireCohort_pathvar$SMAD4[EntireCohort_pathvar$GENE8 %in% c("SMAD4")] <- 1
EntireCohort_pathvar$SMAD4[EntireCohort_pathvar$GENE9 %in% c("SMAD4")] <- 1

EntireCohort_pathvar$SMAD4[is.na(EntireCohort_pathvar$SMAD4)] <- 0
EntireCohort_pathvar$SMAD4 <- as.numeric(EntireCohort_pathvar$SMAD4)
unique(EntireCohort_pathvar$SMAD4)
# 0 1

table(EntireCohort_pathvar$SMAD4)
#      0      1 
# 469795      4




### 'SMARCAL1' 
EntireCohort_pathvar$SMARCAL1 <- NA
EntireCohort_pathvar$SMARCAL1[EntireCohort_pathvar$GENE1 %in% c("SMARCAL1")] <- 1
EntireCohort_pathvar$SMARCAL1[EntireCohort_pathvar$GENE2 %in% c("SMARCAL1")] <- 1
EntireCohort_pathvar$SMARCAL1[EntireCohort_pathvar$GENE3 %in% c("SMARCAL1")] <- 1
EntireCohort_pathvar$SMARCAL1[EntireCohort_pathvar$GENE4 %in% c("SMARCAL1")] <- 1
EntireCohort_pathvar$SMARCAL1[EntireCohort_pathvar$GENE5 %in% c("SMARCAL1")] <- 1
EntireCohort_pathvar$SMARCAL1[EntireCohort_pathvar$GENE6 %in% c("SMARCAL1")] <- 1
EntireCohort_pathvar$SMARCAL1[EntireCohort_pathvar$GENE7 %in% c("SMARCAL1")] <- 1
EntireCohort_pathvar$SMARCAL1[EntireCohort_pathvar$GENE8 %in% c("SMARCAL1")] <- 1
EntireCohort_pathvar$SMARCAL1[EntireCohort_pathvar$GENE9 %in% c("SMARCAL1")] <- 1

EntireCohort_pathvar$SMARCAL1[is.na(EntireCohort_pathvar$SMARCAL1)] <- 0
EntireCohort_pathvar$SMARCAL1 <- as.numeric(EntireCohort_pathvar$SMARCAL1)
unique(EntireCohort_pathvar$SMARCAL1)
# 0 1

table(EntireCohort_pathvar$SMARCAL1)
#      0      1 
# 468114   1685




### 'SOS1'
EntireCohort_pathvar$SOS1 <- NA
EntireCohort_pathvar$SOS1[EntireCohort_pathvar$GENE1 %in% c("SOS1")] <- 1
EntireCohort_pathvar$SOS1[EntireCohort_pathvar$GENE2 %in% c("SOS1")] <- 1
EntireCohort_pathvar$SOS1[EntireCohort_pathvar$GENE3 %in% c("SOS1")] <- 1
EntireCohort_pathvar$SOS1[EntireCohort_pathvar$GENE4 %in% c("SOS1")] <- 1
EntireCohort_pathvar$SOS1[EntireCohort_pathvar$GENE5 %in% c("SOS1")] <- 1
EntireCohort_pathvar$SOS1[EntireCohort_pathvar$GENE6 %in% c("SOS1")] <- 1
EntireCohort_pathvar$SOS1[EntireCohort_pathvar$GENE7 %in% c("SOS1")] <- 1
EntireCohort_pathvar$SOS1[EntireCohort_pathvar$GENE8 %in% c("SOS1")] <- 1
EntireCohort_pathvar$SOS1[EntireCohort_pathvar$GENE9 %in% c("SOS1")] <- 1

EntireCohort_pathvar$SOS1[is.na(EntireCohort_pathvar$SOS1)] <- 0
EntireCohort_pathvar$SOS1 <- as.numeric(EntireCohort_pathvar$SOS1)
unique(EntireCohort_pathvar$SOS1)
# 0 1

table(EntireCohort_pathvar$SOS1)
#      0      1 
# 469781     18




### 'TALDO1'
EntireCohort_pathvar$TALDO1 <- NA
EntireCohort_pathvar$TALDO1[EntireCohort_pathvar$GENE1 %in% c("TALDO1")] <- 1
EntireCohort_pathvar$TALDO1[EntireCohort_pathvar$GENE2 %in% c("TALDO1")] <- 1
EntireCohort_pathvar$TALDO1[EntireCohort_pathvar$GENE3 %in% c("TALDO1")] <- 1
EntireCohort_pathvar$TALDO1[EntireCohort_pathvar$GENE4 %in% c("TALDO1")] <- 1
EntireCohort_pathvar$TALDO1[EntireCohort_pathvar$GENE5 %in% c("TALDO1")] <- 1
EntireCohort_pathvar$TALDO1[EntireCohort_pathvar$GENE6 %in% c("TALDO1")] <- 1
EntireCohort_pathvar$TALDO1[EntireCohort_pathvar$GENE7 %in% c("TALDO1")] <- 1
EntireCohort_pathvar$TALDO1[EntireCohort_pathvar$GENE8 %in% c("TALDO1")] <- 1
EntireCohort_pathvar$TALDO1[EntireCohort_pathvar$GENE9 %in% c("TALDO1")] <- 1

EntireCohort_pathvar$TALDO1[is.na(EntireCohort_pathvar$TALDO1)] <- 0
EntireCohort_pathvar$TALDO1 <- as.numeric(EntireCohort_pathvar$TALDO1)
unique(EntireCohort_pathvar$TALDO1)
# 0 1

table(EntireCohort_pathvar$TALDO1)
#      0      1 
# 468277   1522






### 'TRIM37'
EntireCohort_pathvar$TRIM37 <- NA
EntireCohort_pathvar$TRIM37[EntireCohort_pathvar$GENE1 %in% c("TRIM37")] <- 1
EntireCohort_pathvar$TRIM37[EntireCohort_pathvar$GENE2 %in% c("TRIM37")] <- 1
EntireCohort_pathvar$TRIM37[EntireCohort_pathvar$GENE3 %in% c("TRIM37")] <- 1
EntireCohort_pathvar$TRIM37[EntireCohort_pathvar$GENE4 %in% c("TRIM37")] <- 1
EntireCohort_pathvar$TRIM37[EntireCohort_pathvar$GENE5 %in% c("TRIM37")] <- 1
EntireCohort_pathvar$TRIM37[EntireCohort_pathvar$GENE6 %in% c("TRIM37")] <- 1
EntireCohort_pathvar$TRIM37[EntireCohort_pathvar$GENE7 %in% c("TRIM37")] <- 1
EntireCohort_pathvar$TRIM37[EntireCohort_pathvar$GENE8 %in% c("TRIM37")] <- 1
EntireCohort_pathvar$TRIM37[EntireCohort_pathvar$GENE9 %in% c("TRIM37")] <- 1

EntireCohort_pathvar$TRIM37[is.na(EntireCohort_pathvar$TRIM37)] <- 0
EntireCohort_pathvar$TRIM37 <- as.numeric(EntireCohort_pathvar$TRIM37)
unique(EntireCohort_pathvar$TRIM37)
# 0 1

table(EntireCohort_pathvar$TRIM37)
#      0      1 
# 468767   1032





### 'XRCC4'
EntireCohort_pathvar$XRCC4 <- NA
EntireCohort_pathvar$XRCC4[EntireCohort_pathvar$GENE1 %in% c("XRCC4")] <- 1
EntireCohort_pathvar$XRCC4[EntireCohort_pathvar$GENE2 %in% c("XRCC4")] <- 1
EntireCohort_pathvar$XRCC4[EntireCohort_pathvar$GENE3 %in% c("XRCC4")] <- 1
EntireCohort_pathvar$XRCC4[EntireCohort_pathvar$GENE4 %in% c("XRCC4")] <- 1
EntireCohort_pathvar$XRCC4[EntireCohort_pathvar$GENE5 %in% c("XRCC4")] <- 1
EntireCohort_pathvar$XRCC4[EntireCohort_pathvar$GENE6 %in% c("XRCC4")] <- 1
EntireCohort_pathvar$XRCC4[EntireCohort_pathvar$GENE7 %in% c("XRCC4")] <- 1
EntireCohort_pathvar$XRCC4[EntireCohort_pathvar$GENE8 %in% c("XRCC4")] <- 1
EntireCohort_pathvar$XRCC4[EntireCohort_pathvar$GENE9 %in% c("XRCC4")] <- 1

EntireCohort_pathvar$XRCC4[is.na(EntireCohort_pathvar$XRCC4)] <- 0
EntireCohort_pathvar$XRCC4 <- as.numeric(EntireCohort_pathvar$XRCC4)
unique(EntireCohort_pathvar$XRCC4)
# 0 1

table(EntireCohort_pathvar$XRCC4)
#      0      1 
# 469547    252






# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)

















# CAgene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_CAgene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_CAgene_, "pathvar_CAgene_.csv", row.names=FALSE)
# system('dx upload pathvar_CAgene_.csv --path Emily-folder/results_2/pathvar_CAgene_.csv')


# CAgene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ CAgene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_CAgene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_incidence_1, "pathvar_CAgene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_incidence_1.csv --path Emily-folder/results_2/pathvar_CAgene_incidence_1.csv')


# GOIs
# removed RAF1(n=0), BRAF(n=1), RRAS2(n=1)
LR_2 <- glm(formula = cancer_ROOC_dummy ~ BUB1B+CBL+CREBBP+EP300+ERCC3+FANCA+KDM6A+KRAS+LIG4+LZTR1+NRAS+POLE+PTPN11+RECQL3+RECQL4+RIT1+SBDS+SMAD4+SMARCAL1+SOS1+TALDO1+ TRIM37+XRCC4, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_CAgene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_incidence_2, "pathvar_CAgene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_incidence_2.csv --path Emily-folder/results_2/pathvar_CAgene_incidence_2.csv')






### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "CAgene", 'BUB1B', 'CBL', 'CREBBP', 
                                     'EP300', 'ERCC3', 'FANCA', 'KDM6A', 'KRAS', 'LIG4', 
                                     'LZTR1', 'NRAS', 'POLE', 'PTPN11', 
                                     'RECQL3', 'RECQL4', 'RIT1', 'SBDS', 'SMAD4', 
                                     'SMARCAL1', 'SOS1', 'TALDO1',     
                                     'TRIM37', 'XRCC4')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     32

with(data_setA, table(CAgene))
# CAgene
# 0      1 
# 432972  32909
# 18/9/25

# CAgene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ CAgene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_CAgene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_incidence_3, "pathvar_CAgene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_incidence_3.csv --path Emily-folder/results_2/pathvar_CAgene_incidence_3.csv')


# age + ethnicity + smoking + sex + CAgene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_CAgene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_incidence_4, "pathvar_CAgene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_incidence_4.csv --path Emily-folder/results_2/pathvar_CAgene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ BUB1B+CBL+CREBBP+EP300+ERCC3+FANCA+KDM6A+KRAS+LIG4+LZTR1+NRAS+POLE+PTPN11+RECQL3+RECQL4+RIT1+SBDS+SMAD4+SMARCAL1+SOS1+TALDO1+ TRIM37+XRCC4, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_CAgene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_incidence_5, "pathvar_CAgene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_incidence_5.csv --path Emily-folder/results_2/pathvar_CAgene_incidence_5.csv')


# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B+CBL+CREBBP+EP300+ERCC3+FANCA+KDM6A+KRAS+LIG4+LZTR1+NRAS+POLE+PTPN11+RECQL3+RECQL4+RIT1+SBDS+SMAD4+SMARCAL1+SOS1+TALDO1+ TRIM37+XRCC4, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_CAgene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_incidence_6, "pathvar_CAgene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_incidence_6.csv --path Emily-folder/results_2/pathvar_CAgene_incidence_6.csv')







### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "CAgene", 'BUB1B', 'CBL', 'CREBBP', 'EP300', 
                                     'ERCC3', 'FANCA', 'KDM6A', 'KRAS', 'LIG4', 'LZTR1', 
                                     'NRAS', 'POLE', 'PTPN11', 'RECQL3', 'RECQL4', 
                                     'RIT1', 'SBDS', 'SMAD4', 'SMARCAL1', 'SOS1', 
                                     'TALDO1',     'TRIM37', 'XRCC4')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     37

with(data_setB, table(CAgene))
# CAgene
# 0      1 
# 354320  26839
# 18/9/25

# CAgene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ CAgene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_CAgene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_incidence_7, "pathvar_CAgene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_incidence_7.csv --path Emily-folder/results_2/pathvar_CAgene_incidence_7.csv')



# all variables + CAgene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+CAgene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_CAgene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_incidence_8, "pathvar_CAgene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_incidence_8.csv --path Emily-folder/results_2/pathvar_CAgene_incidence_8.csv')




# GOIs alone
LR_9 <- glm(formula = cancer_ROOC_dummy ~ BUB1B+CBL+CREBBP+EP300+ERCC3+FANCA+KDM6A+KRAS+LIG4+LZTR1+NRAS+POLE+PTPN11+RECQL3+RECQL4+RIT1+SBDS+SMAD4+SMARCAL1+SOS1+TALDO1+ TRIM37+XRCC4, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_CAgene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_incidence_9, "pathvar_CAgene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_incidence_9.csv --path Emily-folder/results_2/pathvar_CAgene_incidence_9.csv')




# all variables + GOIs
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+BUB1B+CBL+CREBBP+EP300+ERCC3+FANCA+KDM6A+KRAS+LIG4+LZTR1+NRAS+POLE+PTPN11+RECQL3+RECQL4+RIT1+SBDS+SMAD4+SMARCAL1+SOS1+TALDO1+ TRIM37+XRCC4, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_CAgene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_incidence_10, "pathvar_CAgene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_incidence_10.csv --path Emily-folder/results_2/pathvar_CAgene_incidence_10.csv')









### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(CAgene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(CAgene, cancer_ROOC_dummy)))
# p-value = 0.5286

with(EntireCohort_pathvar, table(CAgene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(CAgene, Ethnicity_White)))
# p-value = 2.704e-10
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.8171089 0.8984703
# sample estimates:
#   odds ratio 
# 0.8566349

with(EntireCohort_pathvar, table(CAgene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(CAgene, Ever_smoked_dummy)))
# p-value = 0.9119

with(EntireCohort_pathvar, table(CAgene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(CAgene, Sex_Male)))
#p-value = 0.008162
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.007836 1.054086
# sample estimates:
#   odds ratio 
# 1.030734

with(EntireCohort_pathvar, table(CAgene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(CAgene, University_education)))
# p-value = 0.2042

with(EntireCohort_pathvar, table(CAgene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(CAgene, Annual_income))	)
# p-value = 0.7188

with(EntireCohort_pathvar, table(CAgene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(CAgene, Smoker_previous)))
# p-value = 0.07401

with(EntireCohort_pathvar, table(CAgene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(CAgene, Smoker_current))	)
# p-value = 0.8093



### data_setA
with(data_setA, table(CAgene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(CAgene, cancer_ROOC_dummy)))
# p-value = 0.505

with(data_setA, table(CAgene, Ethnicity_White))			
fisher.test(with(data_setA, table(CAgene, Ethnicity_White)))
# p-value = 1.945e-10
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.8155052 0.8969327
# sample estimates:
#   odds ratio 
# 0.8550572 

with(data_setA, table(CAgene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(CAgene, Ever_smoked_dummy)))
# p-value = 0.852

with(data_setA, table(CAgene, Sex_Male))			
fisher.test(with(data_setA, table(CAgene, Sex_Male)))
# p-value = 0.007488
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.008196 1.054692
# sample estimates:
#   odds ratio 
# 1.031203



### data_setB
with(data_setB, table(CAgene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(CAgene, cancer_ROOC_dummy)))
# p-value = 0.2704

with(data_setB, table(CAgene, Ethnicity_White))			
fisher.test(with(data_setB, table(CAgene, Ethnicity_White)))
# p-value = 2.419e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.8106587 0.9092258
# sample estimates:
#   odds ratio 
# 0.858264 

with(data_setB, table(CAgene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(CAgene, Ever_smoked_dummy)))
# p-value = 0.509

with(data_setB, table(CAgene, Sex_Male))			
fisher.test(with(data_setB, table(CAgene, Sex_Male)))
# p-value = 0.09191

with(data_setB, table(CAgene, University_education))			
fisher.test(with(data_setB, table(CAgene, University_education)))
# p-value = 0.1802

with(data_setB, table(CAgene, Annual_income))			
fisher.test(with(data_setB, table(CAgene, Annual_income)))
# p-value = 0.7252






### can plot the continuous variables: 

# make CAgene into factor
df_to_plot <- EntireCohort_pathvar
df_to_plot$CAgene <- as.factor(df_to_plot$CAgene)

# try a violin plot

# Age_at_recruitment  
df_to_plot %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ CAgene, data = df_to_plot)
# W = 7303603142, p-value = 0.01811

# BMI  
df_to_plot %>% 
  ggplot(mapping = aes(x = CAgene, y = BMI)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = CAgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ CAgene, data = df_to_plot)
# W = 7189776912, p-value = 0.9538

# Standing_height 
df_to_plot %>% 
  ggplot(mapping = aes(x = CAgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = CAgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ CAgene, data = df_to_plot)
# W = 7175514608, p-value = 0.2207

# make this a bar plot
# Weekly_exercise
df_to_plot %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(CAgene), scales = "free_y") +
  theme_light()
# 



### data_setA

# make CAgene into factor
data_setA$CAgene <- as.factor(data_setA$CAgene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ CAgene, data = data_setA)
# W = 7178254577, p-value = 0.02179



### data_setB

# make CAgene into factor
data_setB$CAgene <- as.factor(data_setB$CAgene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ CAgene, data = data_setB)
# W = 4784633584, p-value = 0.08582

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = CAgene, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = CAgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ CAgene, data = data_setB)
# W = 4765408926, p-value = 0.5415

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = CAgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = CAgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ CAgene, data = data_setB)
# W = 4749320974, p-value = 0.7526


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(CAgene), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# CAgene Type ------------------------------------------------------------

# make sure the execute the code in the 'CAgenes' section

## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_CAgene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_type_1, "pathvar_CAgene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_type_1.csv --path Emily-folder/results_2/pathvar_CAgene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_CAgene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_type_2, "pathvar_CAgene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_type_2.csv --path Emily-folder/results_2/pathvar_CAgene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_CAgene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_type_3, "pathvar_CAgene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_type_3.csv --path Emily-folder/results_2/pathvar_CAgene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_CAgene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_type_4, "pathvar_CAgene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_type_4.csv --path Emily-folder/results_2/pathvar_CAgene_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_CAgene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_type_5, "pathvar_CAgene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_type_5.csv --path Emily-folder/results_2/pathvar_CAgene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_CAgene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_type_6, "pathvar_CAgene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_type_6.csv --path Emily-folder/results_2/pathvar_CAgene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_CAgene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_type_7, "pathvar_CAgene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_type_7.csv --path Emily-folder/results_2/pathvar_CAgene_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_CAgene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_type_8, "pathvar_CAgene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_type_8.csv --path Emily-folder/results_2/pathvar_CAgene_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_CAgene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_type_9, "pathvar_CAgene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_type_9.csv --path Emily-folder/results_2/pathvar_CAgene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_CAgene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_type_10, "pathvar_CAgene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_type_10.csv --path Emily-folder/results_2/pathvar_CAgene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_CAgene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_type_11, "pathvar_CAgene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_type_11.csv --path Emily-folder/results_2/pathvar_CAgene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_CAgene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_type_12, "pathvar_CAgene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_type_12.csv --path Emily-folder/results_2/pathvar_CAgene_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_CAgene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_type_13, "pathvar_CAgene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_type_13.csv --path Emily-folder/results_2/pathvar_CAgene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_CAgene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_type_14, "pathvar_CAgene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_type_14.csv --path Emily-folder/results_2/pathvar_CAgene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_CAgene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CAgene_type_15, "pathvar_CAgene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_type_15.csv --path Emily-folder/results_2/pathvar_CAgene_type_15.csv')






# RECQL3 Type -------------------------------------------------------------
# repeat the type analysis for just RECQL3


LR_oral_RECQL3 <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL3, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_RECQL3)
# exp(cbind(OR = coef(LR_oral_RECQL3), confint(LR_oral_RECQL3)))
pathvar_RECQL3_type_1 <- tidy(LR_oral_RECQL3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_type_1, "pathvar_RECQL3_type_1.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_type_1.csv --path Emily-folder/results_2/pathvar_RECQL3_type_1.csv')



LR_digestive_RECQL3 <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL3, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_RECQL3)
# exp(cbind(OR = coef(LR_digestive_RECQL3), confint(LR_digestive_RECQL3)))
pathvar_RECQL3_type_2 <- tidy(LR_digestive_RECQL3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_type_2, "pathvar_RECQL3_type_2.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_type_2.csv --path Emily-folder/results_2/pathvar_RECQL3_type_2.csv')



LR_respiratory_RECQL3 <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL3, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_RECQL3)
# exp(cbind(OR = coef(LR_respiratory_RECQL3), confint(LR_respiratory_RECQL3)))
pathvar_RECQL3_type_3 <- tidy(LR_respiratory_RECQL3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_type_3, "pathvar_RECQL3_type_3.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_type_3.csv --path Emily-folder/results_2/pathvar_RECQL3_type_3.csv')



LR_bone_RECQL3 <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL3, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_RECQL3)
# exp(cbind(OR = coef(LR_bone_RECQL3), confint(LR_bone_RECQL3)))
pathvar_RECQL3_type_4 <- tidy(LR_bone_RECQL3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_type_4, "pathvar_RECQL3_type_4.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_type_4.csv --path Emily-folder/results_2/pathvar_RECQL3_type_4.csv')



LR_skin_RECQL3 <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL3, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_RECQL3)
# exp(cbind(OR = coef(LR_skin_RECQL3), confint(LR_skin_RECQL3)))
pathvar_RECQL3_type_5 <- tidy(LR_skin_RECQL3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_type_5, "pathvar_RECQL3_type_5.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_type_5.csv --path Emily-folder/results_2/pathvar_RECQL3_type_5.csv')



LR_mesothelium_RECQL3 <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL3, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_RECQL3)
# exp(cbind(OR = coef(LR_mesothelium_RECQL3), confint(LR_mesothelium_RECQL3)))
pathvar_RECQL3_type_6 <- tidy(LR_mesothelium_RECQL3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_type_6, "pathvar_RECQL3_type_6.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_type_6.csv --path Emily-folder/results_2/pathvar_RECQL3_type_6.csv')



LR_breast_RECQL3 <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL3, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_RECQL3)
# exp(cbind(OR = coef(LR_breast_RECQL3), confint(LR_breast_RECQL3)))
pathvar_RECQL3_type_7 <- tidy(LR_breast_RECQL3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_type_7, "pathvar_RECQL3_type_7.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_type_7.csv --path Emily-folder/results_2/pathvar_RECQL3_type_7.csv')



LR_femalegenital_RECQL3 <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL3, 
                               data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_RECQL3)
# exp(cbind(OR = coef(LR_femalegenital_RECQL3), confint(LR_femalegenital_RECQL3)))
pathvar_RECQL3_type_8 <- tidy(LR_femalegenital_RECQL3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_type_8, "pathvar_RECQL3_type_8.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_type_8.csv --path Emily-folder/results_2/pathvar_RECQL3_type_8.csv')



LR_malegenital_RECQL3 <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL3, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_RECQL3)
# exp(cbind(OR = coef(LR_malegenital_RECQL3), confint(LR_malegenital_RECQL3)))
pathvar_RECQL3_type_9 <- tidy(LR_malegenital_RECQL3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_type_9, "pathvar_RECQL3_type_9.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_type_9.csv --path Emily-folder/results_2/pathvar_RECQL3_type_9.csv')



LR_urinary_RECQL3 <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL3, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_RECQL3)
# exp(cbind(OR = coef(LR_urinary_RECQL3), confint(LR_urinary_RECQL3)))
pathvar_RECQL3_type_10 <- tidy(LR_urinary_RECQL3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_type_10, "pathvar_RECQL3_type_10.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_type_10.csv --path Emily-folder/results_2/pathvar_RECQL3_type_10.csv')



LR_cns_RECQL3 <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL3, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_RECQL3)
# exp(cbind(OR = coef(LR_cns_RECQL3), confint(LR_cns_RECQL3)))
pathvar_RECQL3_type_11 <- tidy(LR_cns_RECQL3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_type_11, "pathvar_RECQL3_type_11.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_type_11.csv --path Emily-folder/results_2/pathvar_RECQL3_type_11.csv')



LR_endocrine_RECQL3 <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL3, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_RECQL3)
# exp(cbind(OR = coef(LR_endocrine_RECQL3), confint(LR_endocrine_RECQL3)))
pathvar_RECQL3_type_12 <- tidy(LR_endocrine_RECQL3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_type_12, "pathvar_RECQL3_type_12.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_type_12.csv --path Emily-folder/results_2/pathvar_RECQL3_type_12.csv')



LR_lymphatic_RECQL3 <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL3, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_RECQL3)
# exp(cbind(OR = coef(LR_lymphatic_RECQL3), confint(LR_lymphatic_RECQL3)))
pathvar_RECQL3_type_13 <- tidy(LR_lymphatic_RECQL3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_type_13, "pathvar_RECQL3_type_13.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_type_13.csv --path Emily-folder/results_2/pathvar_RECQL3_type_13.csv')



LR_secondary_RECQL3 <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL3, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_RECQL3)
# exp(cbind(OR = coef(LR_secondary_RECQL3), confint(LR_secondary_RECQL3)))
pathvar_RECQL3_type_14 <- tidy(LR_secondary_RECQL3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_type_14, "pathvar_RECQL3_type_14.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_type_14.csv --path Emily-folder/results_2/pathvar_RECQL3_type_14.csv')



LR_other_RECQL3 <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL3, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_RECQL3)
# exp(cbind(OR = coef(LR_other_RECQL3), confint(LR_other_RECQL3)))
pathvar_RECQL3_type_15 <- tidy(LR_other_RECQL3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_type_15, "pathvar_RECQL3_type_15.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_type_15.csv --path Emily-folder/results_2/pathvar_RECQL3_type_15.csv')












# CAgene Age -------------------------------------------------------------

# make sure the execute the code in the 'CAgenes' section

summary(EntireCohort_pathvar)

### linear regression
# "CAgene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ CAgene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_CAgene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_CAgene_age_1, "pathvar_CAgene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_age_1.csv --path Emily-folder/results_2/pathvar_CAgene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ BUB1B+CBL+CREBBP+EP300+ERCC3+FANCA+KDM6A+KRAS+LIG4+LZTR1+NRAS+POLE+PTPN11+RECQL3+RECQL4+RIT1+SBDS+SMAD4+SMARCAL1+SOS1+TALDO1+ TRIM37+XRCC4, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_CAgene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_CAgene_age_2, "pathvar_CAgene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_age_2.csv --path Emily-folder/results_2/pathvar_CAgene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "CAgene", 'BUB1B', 'CBL', 'CREBBP', 'EP300', 
                                     'ERCC3', 'FANCA', 'KDM6A', 'LIG4', 'LZTR1', 
                                     'NRAS', 'POLE', 'PTPN11', 'RECQL3', 'RECQL4', 
                                     'RIT1', 'SBDS', 'SMARCAL1', 'SOS1', 'TALDO1', 
                                     'TRIM37', 'XRCC4')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     35

# make sure genes with <2 are removed
with(data_setD, table(CAgene))
# CAgene
# 0      1 
# 106251   8131
# remving KRAS and SMAD4 bc only 1 participant


data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "CAgene", 'BUB1B', 'CBL', 
                                     'CREBBP', 'EP300', 'ERCC3', 'FANCA', 'KDM6A', 'LIG4', 
                                     'LZTR1', 'NRAS', 'POLE', 'PTPN11', 'RECQL3', 
                                     'RECQL4', 'RIT1', 'SBDS', 'SMARCAL1', 'SOS1', 
                                     'TALDO1',     'TRIM37', 'XRCC4')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    35

# make sure genes with <2 are removed
with(data_setE, table(CAgene))
#CAgene
# 0     1 
# 85552  6561
# didn't have to remove any more genes


### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ CAgene, data = data_setD)
# summary(LM_15)
pathvar_CAgene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_CAgene_age_3, "pathvar_CAgene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_age_3.csv --path Emily-folder/results_2/pathvar_CAgene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CAgene, data = data_setD)
# summary(LM_16)
pathvar_CAgene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_CAgene_age_4, "pathvar_CAgene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_age_4.csv --path Emily-folder/results_2/pathvar_CAgene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ BUB1B+CBL+CREBBP+EP300+ERCC3+FANCA+KDM6A+LIG4+LZTR1+NRAS+POLE+PTPN11+RECQL3+RECQL4+RIT1+SBDS+SMARCAL1+SOS1+TALDO1+ TRIM37+XRCC4, data = data_setD)
# summary(LM_17)
pathvar_CAgene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_CAgene_age_5, "pathvar_CAgene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_age_5.csv --path Emily-folder/results_2/pathvar_CAgene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B+CBL+CREBBP+EP300+ERCC3+FANCA+KDM6A+LIG4+LZTR1+NRAS+POLE+PTPN11+RECQL3+RECQL4+RIT1+SBDS+SMARCAL1+SOS1+TALDO1+ TRIM37+XRCC4, data = data_setD)
# summary(LM_18)
pathvar_CAgene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_CAgene_age_6, "pathvar_CAgene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_age_6.csv --path Emily-folder/results_2/pathvar_CAgene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ CAgene, data = data_setE)
# summary(LM_15)
pathvar_CAgene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_CAgene_age_7, "pathvar_CAgene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_age_7.csv --path Emily-folder/results_2/pathvar_CAgene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+CAgene, data = data_setE)
# summary(LM_16)
pathvar_CAgene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_CAgene_age_8, "pathvar_CAgene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_age_8.csv --path Emily-folder/results_2/pathvar_CAgene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ BUB1B+CBL+CREBBP+EP300+ERCC3+FANCA+KDM6A+LIG4+LZTR1+NRAS+POLE+PTPN11+RECQL3+RECQL4+RIT1+SBDS+SMARCAL1+SOS1+TALDO1+ TRIM37+XRCC4, data = data_setE)
# summary(LM_17)
pathvar_CAgene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_CAgene_age_9, "pathvar_CAgene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_age_9.csv --path Emily-folder/results_2/pathvar_CAgene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+BUB1B+CBL+CREBBP+EP300+ERCC3+FANCA+KDM6A+LIG4+LZTR1+NRAS+POLE+PTPN11+RECQL3+RECQL4+RIT1+SBDS+SMARCAL1+SOS1+TALDO1+ TRIM37+XRCC4, data = data_setE)
# summary(LM_18)
pathvar_CAgene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_CAgene_age_10, "pathvar_CAgene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_CAgene_age_10.csv --path Emily-folder/results_2/pathvar_CAgene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(df_cancer_diagnosis, table(CAgene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(CAgene, Ethnicity_White)))
# p-value = 0.1557

with(df_cancer_diagnosis, table(CAgene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(CAgene, Ever_smoked_dummy)))
# p-value = 0.2243

with(df_cancer_diagnosis, table(CAgene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(CAgene, Sex_Male)))
# p-value = 0.28

with(df_cancer_diagnosis, table(CAgene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(CAgene, University_education)))
# p-value = 0.2231

with(df_cancer_diagnosis, table(CAgene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(CAgene, Annual_income)))
# p-value = 0.3648


### data_setD
with(data_setD, table(CAgene, Ethnicity_White))			
fisher.test(with(data_setD, table(CAgene, Ethnicity_White)))
# p-value = 0.1754

with(data_setD, table(CAgene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(CAgene, Ever_smoked_dummy)))
# p-value = 0.219

with(data_setD, table(CAgene, Sex_Male))			
fisher.test(with(data_setD, table(CAgene, Sex_Male)))
# p-value = 0.2678


### data_setE
with(data_setE, table(CAgene, Ethnicity_White))			
fisher.test(with(data_setE, table(CAgene, Ethnicity_White)))
# p-value = 0.1188

with(data_setE, table(CAgene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(CAgene, Ever_smoked_dummy)))
# p-value = 0.1224

with(data_setE, table(CAgene, Sex_Male))			
fisher.test(with(data_setE, table(CAgene, Sex_Male)))
# p-value = 0.6262

with(data_setE, table(CAgene, University_education))			
fisher.test(with(data_setE, table(CAgene, University_education)))
# p-value = 0.1314

with(data_setE, table(CAgene, Annual_income))			
fisher.test(with(data_setE, table(CAgene, Annual_income)))
# p-value = 0.503



## make plots for the linear regression data


# make CAgene into factor
df_cancer_diagnosis$CAgene <- as.factor(df_cancer_diagnosis$CAgene)

# try a violin plot


# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ CAgene, data = df_cancer_diagnosis)
# W = 437387774, p-value = 0.4867



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ CAgene, data = df_cancer_diagnosis)
# W = 437981691, p-value = 0.6228



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CAgene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CAgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ CAgene, data = df_cancer_diagnosis)
#W = 439689156, p-value = 0.1967



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CAgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = CAgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ CAgene, data = df_cancer_diagnosis)
# W = 437042885, p-value = 0.9857



# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(CAgene), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$CAgene <- as.factor(data_setD$CAgene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ CAgene, data = data_setD)
# W = 430237028, p-value = 0.5474



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ CAgene, data = data_setD)
# W = 430838554, p-value = 0.6947




### data_setE

data_setE$CAgene <- as.factor(data_setE$CAgene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ CAgene, data = data_setE)
# W = 278423608, p-value = 0.2827



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = CAgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ CAgene, data = data_setE)
#W = 279452906, p-value = 0.5626



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = CAgene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = CAgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ CAgene, data = data_setE)
# W = 284770056, p-value = 0.04734



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = CAgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = CAgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ CAgene, data = data_setE)
# W = 283396739, p-value = 0.1861


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(CAgene), scales = "free_y") +
  theme_light()
# 









# RECQL3 extras -----------------------------------------------------------



# look at the demographics of the RECQL3 participants

### Fisher test on these datasets

### EntireCohort_pathvar
with(EntireCohort_pathvar, table(RECQL3, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(RECQL3, cancer_ROOC_dummy)))
# p-value = 0.002192
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.064794 1.335259
# sample estimates:
#   odds ratio 
# 1.193316

with(EntireCohort_pathvar, table(RECQL3, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(RECQL3, Ethnicity_White)))
# p-value = 0.07965

with(EntireCohort_pathvar, table(RECQL3, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(RECQL3, Ever_smoked_dummy)))
# p-value = 0.8568

with(EntireCohort_pathvar, table(RECQL3, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(RECQL3, Sex_Male)))
# p-value = 0.1351

with(EntireCohort_pathvar, table(RECQL3, University_education))			
fisher.test(with(EntireCohort_pathvar, table(RECQL3, University_education)))
# p-value = 0.8283

with(EntireCohort_pathvar, table(RECQL3, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(RECQL3, Annual_income))	)
# p-value = 0.6885

with(EntireCohort_pathvar, table(RECQL3, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(RECQL3, Smoker_previous)))
# p-value = 0.8319

with(EntireCohort_pathvar, table(RECQL3, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(RECQL3, Smoker_current))	)
# p-value = 0.5931



### data_setA
with(data_setA, table(RECQL3, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(RECQL3, cancer_ROOC_dummy)))
# p-value = 0.001935
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.067815 1.339556
# sample estimates:
#   odds ratio 
# 1.196916

with(data_setA, table(RECQL3, Ethnicity_White))			
fisher.test(with(data_setA, table(RECQL3, Ethnicity_White)))
# p-value = 0.06953

with(data_setA, table(RECQL3, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(RECQL3, Ever_smoked_dummy)))
# p-value = 0.877

with(data_setA, table(RECQL3, Sex_Male))			
fisher.test(with(data_setA, table(RECQL3, Sex_Male)))
# p-value = 0.1627



### data_setB
with(data_setB, table(RECQL3, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(RECQL3, cancer_ROOC_dummy)))
# p-value = 0.0006975
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.093177 1.401695
# sample estimates:
#   odds ratio 
# 1.238962

with(data_setB, table(RECQL3, Ethnicity_White))			
fisher.test(with(data_setB, table(RECQL3, Ethnicity_White)))
# p-value = 0.05778

with(data_setB, table(RECQL3, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(RECQL3, Ever_smoked_dummy)))
# p-value = 0.7333

with(data_setB, table(RECQL3, Sex_Male))			
fisher.test(with(data_setB, table(RECQL3, Sex_Male)))
# p-value = 0.5047

with(data_setB, table(RECQL3, University_education))			
fisher.test(with(data_setB, table(RECQL3, University_education)))
# p-value = 0.5239

with(data_setB, table(RECQL3, Annual_income))			
fisher.test(with(data_setB, table(RECQL3, Annual_income)))
# p-value = 0.8012




### can plot the continuous variables: 

# make RECQL3 into factor
df_to_plot <- EntireCohort_pathvar
df_to_plot$RECQL3 <- as.factor(df_to_plot$RECQL3)

# try a violin plot

# Age_at_recruitment  
df_to_plot %>% 
  ggplot(mapping = aes(x = RECQL3, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = RECQL3, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ RECQL3, data = df_to_plot)
# W = 369721387, p-value = 0.853

# BMI  
df_to_plot %>% 
  ggplot(mapping = aes(x = RECQL3, y = BMI)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = RECQL3, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ RECQL3, data = df_to_plot)
# W = 363637078, p-value = 0.7119

# Standing_height 
df_to_plot %>% 
  ggplot(mapping = aes(x = RECQL3, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = RECQL3, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ RECQL3, data = df_to_plot)
# W = 370543660, p-value = 0.4469

# make this a bar plot
# Weekly_exercise
df_to_plot %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(RECQL3), scales = "free_y") +
  theme_light()
# 



### data_setA

# make RECQL3 into factor
data_setA$RECQL3 <- as.factor(data_setA$RECQL3)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = RECQL3, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = RECQL3, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ RECQL3, data = data_setA)
# W = 364967998, p-value = 0.8586



### data_setB

# make RECQL3 into factor
data_setB$RECQL3 <- as.factor(data_setB$RECQL3)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = RECQL3, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = RECQL3, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ RECQL3, data = data_setB)
# W = 246214463, p-value = 0.8608

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = RECQL3, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = RECQL3, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ RECQL3, data = data_setB)
# W = 243750026, p-value = 0.4252

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = RECQL3, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = RECQL3, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ RECQL3, data = data_setB)
# W = 251590146, p-value = 0.237


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(RECQL3), scales = "free_y") +
  theme_light()
# follows the usual profile















# RECQL3 p.Q548X ----------------------------------------------------------


# -	Look at the p.Q548X variant (there have been quite a few studies on this variant and cancer risk)

unique(EntireCohort_pathvar$RECQL3)
# 0 1

RECQL3_df <- EntireCohort_pathvar[EntireCohort_pathvar$RECQL3 == 1,]
# 1575 obs. 

# p.Q548X = 15:90761015:C:T_T (GRCh38)

# separate the var names into 9 columns and then use these to make a new column
RECQL3_df <- separate(data = RECQL3_df, col = VAR, into = c("VARIANT1","VARIANT2","VARIANT3","VARIANT4","VARIANT5","VARIANT6","VARIANT7","VARIANT8","VARIANT9"), sep = ";")

# find participants with the 15:90761015:C:T_T variant
VARIANT1_df <- RECQL3_df[RECQL3_df$VARIANT1 == "15:90761015:C:T_T",]
VARIANT1_df$eid
# # [1] 2520015 3637195 4839300 5937961 2690742 2740419 1154132 1422480 3477890 2364041 3247763 3832791 2281401 2729683 5896647 1077655
# # [17] 4021672 3662351 5637750 3158925 2615302 1068710 3095718 4100571 1161396 5312744 1504704 5933176 3125046 3770061 4435376 2319402
# # [33] 5248709 1895330

VARIANT2_df <- RECQL3_df[RECQL3_df$VARIANT2 == "15:90761015:C:T_T",]
unique(VARIANT2_df$eid)
# [1]      NA 3287818 4008709 3122774 2182357 1792079 2998057 1937090 4083551 2088360 5932824
# [12] 2684750 4931468 3422053 1538819 2995926 2769336 3718607 4626623 4565898 2987048 1492191
# [23] 2482945 4346259 3286283 1278579 2805615 1005157 1586758 2821282 4230882

VARIANT3_df <- RECQL3_df[RECQL3_df$VARIANT3 == "15:90761015:C:T_T",]
unique(VARIANT3_df$eid)
# [1]      NA 1697034 3971758 5247705 4306266 3308631 1523633 4768106 4483533 2708719 3194129
# [12] 4406687 1000686 3416505 1790019 4084183 3759189 2107105 3552575 3467572

VARIANT4_df <- RECQL3_df[RECQL3_df$VARIANT4 == "15:90761015:C:T_T",]
unique(VARIANT4_df$eid)
# [1]      NA 2267832 2822347 5649697 4696204

VARIANT5_df <- RECQL3_df[RECQL3_df$VARIANT5 == "15:90761015:C:T_T",]
unique(VARIANT5_df$eid)
# [1]      NA 2702585 5860773

# none in the rest of the VAR columns



Q548X_eid <- c(2520015, 3637195, 4839300, 5937961, 2690742, 2740419, 1154132, 1422480, 3477890, 2364041, 3247763, 
               3832791, 2281401, 2729683, 5896647, 1077655, 4021672, 3662351, 5637750, 3158925, 2615302, 1068710, 
               3095718, 4100571, 1161396, 5312744, 1504704, 5933176, 3125046, 3770061, 4435376, 2319402, 5248709, 
               1895330, 3287818, 4008709, 3122774, 2182357, 1792079, 2998057, 1937090, 4083551, 2088360, 5932824,               
               2684750, 4931468, 3422053, 1538819, 2995926, 2769336, 3718607, 4626623, 4565898, 2987048, 1492191,
               2482945, 4346259, 3286283, 1278579, 2805615, 1005157, 1586758, 2821282, 4230882,
               1697034, 3971758, 5247705, 4306266, 3308631, 1523633, 4768106, 4483533, 2708719, 3194129,
               4406687, 1000686, 3416505, 1790019, 4084183, 3759189, 2107105, 3552575, 3467572,
               2267832, 2822347, 5649697, 4696204, 2702585, 5860773)

# 89


# annotate the EntireCohort_pathvar df with this info
Q548X_df <- data.frame(eid = Q548X_eid,
                       Q548X = 1)

# merge this with the EntireCohort_pathvar df
EntireCohort_pathvar_Q548X <- full_join(EntireCohort_pathvar, Q548X_df, by = 'eid')

EntireCohort_pathvar_Q548X$Q548X[is.na(EntireCohort_pathvar_Q548X$Q548X)] <- 0

unique(EntireCohort_pathvar_Q548X$Q548X)
# 0 1





# incidence:
LR_1_test <- glm(formula = cancer_ROOC_dummy ~ Q548X, data = EntireCohort_pathvar_Q548X, family = binomial)
# summary(LR_1_test)
# Coefficients:
#   Estimate Std. Error  z value Pr(>|z|)    
# (Intercept) -1.174754   0.003436 -341.872   <2e-16 ***
#   Q548X        0.061104   0.245747    0.249    0.804 
# exp(cbind(OR = coef(LR_1_test), confint(LR_1_test)))
pathvar_RECQL3_Q548X_incidence_1 <- tidy(LR_1_test, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_incidence_1, "pathvar_RECQL3_Q548X_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_incidence_1.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_incidence_1.csv')


data_setA <- EntireCohort_pathvar_Q548X[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                           "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                           "CAgene", 'BUB1B', 'CBL', 'CREBBP', 
                                           'EP300', 'ERCC3', 'FANCA', 'KDM6A', 'KRAS', 'LIG4', 
                                           'LZTR1', 'NRAS', 'POLE', 'PTPN11', 
                                           'RECQL3', 'RECQL4', 'RIT1', 'SBDS', 'SMAD4', 
                                           'SMARCAL1', 'SOS1', 'TALDO1',     
                                           'TRIM37', 'XRCC4', 'Q548X')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     33

LR_4_test <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, 
                 data = data_setA, family = binomial)
# summary(LR_4_test)
# # p-value = 0.813
# exp(cbind(OR = coef(LR_4_test), confint(LR_4_test)))
pathvar_RECQL3_Q548X_incidence_4 <- tidy(LR_4_test, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_incidence_4, "pathvar_RECQL3_Q548X_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_incidence_4.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_incidence_4.csv')


data_setB <- EntireCohort_pathvar_Q548X[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                           "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                           "BMI", "Standing_height", "University_education",
                                           "Annual_income", "Weekly_exercise",
                                           "CAgene", 'BUB1B', 'CBL', 'CREBBP', 'EP300', 
                                           'ERCC3', 'FANCA', 'KDM6A', 'KRAS', 'LIG4', 'LZTR1', 
                                           'NRAS', 'POLE', 'PTPN11', 'RECQL3', 'RECQL4', 
                                           'RIT1', 'SBDS', 'SMAD4', 'SMARCAL1', 'SOS1', 
                                           'TALDO1',     'TRIM37', 'XRCC4', 'Q548X')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     38

LR_8_test <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+Q548X, 
                 data = data_setB, family = binomial)
# summary(LR_8_test)
# # p = 0.5440
# exp(cbind(OR = coef(LR_8_test), confint(LR_8_test)))
pathvar_RECQL3_Q548X_incidence_8 <- tidy(LR_8_test, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_incidence_8, "pathvar_RECQL3_Q548X_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_incidence_8.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_incidence_8.csv')




# age:
LM_1_test <- lm(Age_at_cancer_diagnosis_earliest ~ Q548X, data = EntireCohort_pathvar_Q548X)
# summary(LM_1_test)
# 
pathvar_RECQL3_Q548X_age_1 <- tidy(LM_1_test, conf.int = TRUE)
write.csv(pathvar_RECQL3_Q548X_age_1, "pathvar_RECQL3_Q548X_age_1.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_age_1.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_age_1.csv')


data_setD <- EntireCohort_pathvar_Q548X[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                           "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                           "CAgene", 'BUB1B', 'CBL', 'CREBBP', 'EP300', 
                                           'ERCC3', 'FANCA', 'KDM6A', 'LIG4', 'LZTR1', 
                                           'NRAS', 'POLE', 'PTPN11', 'RECQL3', 'RECQL4', 
                                           'RIT1', 'SBDS', 'SMARCAL1', 'SOS1', 'TALDO1', 
                                           'TRIM37', 'XRCC4', 'Q548X')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     31

LM_4_test <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, data = data_setD)
# summary(LM_4_test)
# 
pathvar_RECQL3_Q548X_age_4 <- tidy(LM_4_test, conf.int = TRUE)
write.csv(pathvar_RECQL3_Q548X_age_4, "pathvar_RECQL3_Q548X_age_4.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_age_4.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_age_4.csv')

data_setE <- EntireCohort_pathvar_Q548X[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                           "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                           "Standing_height", "University_education", "Annual_income",
                                           "Weekly_exercise", "CAgene", 'BUB1B', 'CBL', 
                                           'CREBBP', 'EP300', 'ERCC3', 'FANCA', 'KDM6A', 'LIG4', 
                                           'LZTR1', 'NRAS', 'POLE', 'PTPN11', 'RECQL3', 
                                           'RECQL4', 'RIT1', 'SBDS', 'SMARCAL1', 'SOS1', 
                                           'TALDO1',     'TRIM37', 'XRCC4', 'Q548X')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    36

LM_16_test <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+Q548X, data = data_setE)
# summary(LM_16_test)
# 
pathvar_RECQL3_Q548X_age_16 <- tidy(LM_16_test, conf.int = TRUE)
write.csv(pathvar_RECQL3_Q548X_age_16, "pathvar_RECQL3_Q548X_age_16.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_age_16.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_age_16.csv')




## type 


data_setA2 <- EntireCohort_pathvar_Q548X[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                            "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                            "CAgene", 'BUB1B', 'CBL', 'CREBBP', 
                                            'EP300', 'ERCC3', 'FANCA', 'KDM6A', 'KRAS', 'LIG4', 
                                            'LZTR1', 'NRAS', 'POLE', 'PTPN11', 
                                            'RECQL3', 'RECQL4', 'RIT1', 'SBDS', 'SMAD4', 
                                            'SMARCAL1', 'SOS1', 'TALDO1',     
                                            'TRIM37', 'XRCC4', 'Q548X', 'ICD_code')]
data_setA2 <- data_setA2[complete.cases(data_setA2), ]
dim(data_setA2)
# 465881     34


# full type analysis for this variant


LR_oral_RECQL3_Q548X <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, 
                            data = EntireCohort_pathvar_Q548X, family = binomial)
# summary(LR_oral_RECQL3_Q548X)
# exp(cbind(OR = coef(LR_oral_RECQL3_Q548X), confint(LR_oral_RECQL3_Q548X)))
pathvar_RECQL3_Q548X_type_1 <- tidy(LR_oral_RECQL3_Q548X, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_type_1, "pathvar_RECQL3_Q548X_type_1.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_type_1.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_type_1.csv')



LR_digestive_RECQL3_Q548X <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, 
                                 data = EntireCohort_pathvar_Q548X, family = binomial)
# summary(LR_digestive_RECQL3_Q548X)
# exp(cbind(OR = coef(LR_digestive_RECQL3_Q548X), confint(LR_digestive_RECQL3_Q548X)))
pathvar_RECQL3_Q548X_type_2 <- tidy(LR_digestive_RECQL3_Q548X, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_type_2, "pathvar_RECQL3_Q548X_type_2.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_type_2.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_type_2.csv')



LR_respiratory_RECQL3_Q548X <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, 
                                   data = EntireCohort_pathvar_Q548X, family = binomial)
# summary(LR_respiratory_RECQL3_Q548X)
# exp(cbind(OR = coef(LR_respiratory_RECQL3_Q548X), confint(LR_respiratory_RECQL3_Q548X)))
pathvar_RECQL3_Q548X_type_3 <- tidy(LR_respiratory_RECQL3_Q548X, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_type_3, "pathvar_RECQL3_Q548X_type_3.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_type_3.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_type_3.csv')



LR_bone_RECQL3_Q548X <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, 
                            data = EntireCohort_pathvar_Q548X, family = binomial)
# summary(LR_bone_RECQL3_Q548X)
# exp(cbind(OR = coef(LR_bone_RECQL3_Q548X), confint(LR_bone_RECQL3_Q548X)))
pathvar_RECQL3_Q548X_type_4 <- tidy(LR_bone_RECQL3_Q548X, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_type_4, "pathvar_RECQL3_Q548X_type_4.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_type_4.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_type_4.csv')



LR_skin_RECQL3_Q548X <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, 
                            data = EntireCohort_pathvar_Q548X, family = binomial)
# summary(LR_skin_RECQL3_Q548X)
# exp(cbind(OR = coef(LR_skin_RECQL3_Q548X), confint(LR_skin_RECQL3_Q548X)))
pathvar_RECQL3_Q548X_type_5 <- tidy(LR_skin_RECQL3_Q548X, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_type_5, "pathvar_RECQL3_Q548X_type_5.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_type_5.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_type_5.csv')



LR_mesothelium_RECQL3_Q548X <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, 
                                   data = EntireCohort_pathvar_Q548X, family = binomial)
# summary(LR_mesothelium_RECQL3_Q548X)
# exp(cbind(OR = coef(LR_mesothelium_RECQL3_Q548X), confint(LR_mesothelium_RECQL3_Q548X)))
pathvar_RECQL3_Q548X_type_6 <- tidy(LR_mesothelium_RECQL3_Q548X, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_type_6, "pathvar_RECQL3_Q548X_type_6.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_type_6.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_type_6.csv')



LR_breast_RECQL3_Q548X <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, 
                              data = EntireCohort_pathvar_Q548X, family = binomial)
# summary(LR_breast_RECQL3_Q548X)
# exp(cbind(OR = coef(LR_breast_RECQL3_Q548X), confint(LR_breast_RECQL3_Q548X)))
pathvar_RECQL3_Q548X_type_7 <- tidy(LR_breast_RECQL3_Q548X, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_type_7, "pathvar_RECQL3_Q548X_type_7.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_type_7.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_type_7.csv')



LR_femalegenital_RECQL3_Q548X <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, 
                                     data = EntireCohort_pathvar_Q548X, family = binomial)
# summary(LR_femalegenital_RECQL3_Q548X)
# exp(cbind(OR = coef(LR_femalegenital_RECQL3_Q548X), confint(LR_femalegenital_RECQL3_Q548X)))
pathvar_RECQL3_Q548X_type_8 <- tidy(LR_femalegenital_RECQL3_Q548X, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_type_8, "pathvar_RECQL3_Q548X_type_8.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_type_8.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_type_8.csv')



LR_malegenital_RECQL3_Q548X <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, 
                                   data = EntireCohort_pathvar_Q548X, family = binomial)
# summary(LR_malegenital_RECQL3_Q548X)
# exp(cbind(OR = coef(LR_malegenital_RECQL3_Q548X), confint(LR_malegenital_RECQL3_Q548X)))
pathvar_RECQL3_Q548X_type_9 <- tidy(LR_malegenital_RECQL3_Q548X, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_type_9, "pathvar_RECQL3_Q548X_type_9.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_type_9.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_type_9.csv')



LR_urinary_RECQL3_Q548X <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, 
                               data = EntireCohort_pathvar_Q548X, family = binomial)
# summary(LR_urinary_RECQL3_Q548X)
# exp(cbind(OR = coef(LR_urinary_RECQL3_Q548X), confint(LR_urinary_RECQL3_Q548X)))
pathvar_RECQL3_Q548X_type_10 <- tidy(LR_urinary_RECQL3_Q548X, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_type_10, "pathvar_RECQL3_Q548X_type_10.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_type_10.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_type_10.csv')



LR_cns_RECQL3_Q548X <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, 
                           data = EntireCohort_pathvar_Q548X, family = binomial)
# summary(LR_cns_RECQL3_Q548X)
# exp(cbind(OR = coef(LR_cns_RECQL3_Q548X), confint(LR_cns_RECQL3_Q548X)))
pathvar_RECQL3_Q548X_type_11 <- tidy(LR_cns_RECQL3_Q548X, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_type_11, "pathvar_RECQL3_Q548X_type_11.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_type_11.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_type_11.csv')



LR_endocrine_RECQL3_Q548X <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, 
                                 data = EntireCohort_pathvar_Q548X, family = binomial)
# summary(LR_endocrine_RECQL3_Q548X)
# exp(cbind(OR = coef(LR_endocrine_RECQL3_Q548X), confint(LR_endocrine_RECQL3_Q548X)))
pathvar_RECQL3_Q548X_type_12 <- tidy(LR_endocrine_RECQL3_Q548X, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_type_12, "pathvar_RECQL3_Q548X_type_12.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_type_12.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_type_12.csv')



LR_lymphatic_RECQL3_Q548X <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, 
                                 data = EntireCohort_pathvar_Q548X, family = binomial)
# summary(LR_lymphatic_RECQL3_Q548X)
# exp(cbind(OR = coef(LR_lymphatic_RECQL3_Q548X), confint(LR_lymphatic_RECQL3_Q548X)))
pathvar_RECQL3_Q548X_type_13 <- tidy(LR_lymphatic_RECQL3_Q548X, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_type_13, "pathvar_RECQL3_Q548X_type_13.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_type_13.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_type_13.csv')



LR_secondary_RECQL3_Q548X <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, 
                                 data = EntireCohort_pathvar_Q548X, family = binomial)
# summary(LR_secondary_RECQL3_Q548X)
# exp(cbind(OR = coef(LR_secondary_RECQL3_Q548X), confint(LR_secondary_RECQL3_Q548X)))
pathvar_RECQL3_Q548X_type_14 <- tidy(LR_secondary_RECQL3_Q548X, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_type_14, "pathvar_RECQL3_Q548X_type_14.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_type_14.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_type_14.csv')



LR_other_RECQL3_Q548X <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Q548X, 
                             data = EntireCohort_pathvar_Q548X, family = binomial)
# summary(LR_other_RECQL3_Q548X)
# exp(cbind(OR = coef(LR_other_RECQL3_Q548X), confint(LR_other_RECQL3_Q548X)))
pathvar_RECQL3_Q548X_type_15 <- tidy(LR_other_RECQL3_Q548X, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_Q548X_type_15, "pathvar_RECQL3_Q548X_type_15.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_Q548X_type_15.csv --path Emily-folder/results_2/pathvar_RECQL3_Q548X_type_15.csv')









## looking at proportions of variants in the RECQL3 group

# all vars
RECQL3_df <- EntireCohort_pathvar[EntireCohort_pathvar$RECQL3 == 1,]
ALL <- RECQL3_df$VAR


# all setA vars
data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "CAgene", 'BUB1B', 'CBL', 'CREBBP', 
                                     'EP300', 'ERCC3', 'FANCA', 'KDM6A', 'KRAS', 'LIG4', 
                                     'LZTR1', 'NRAS', 'POLE', 'PTPN11', 
                                     'RECQL3', 'RECQL4', 'RIT1', 'SBDS', 'SMAD4', 
                                     'SMARCAL1', 'SOS1', 'TALDO1',     
                                     'TRIM37', 'XRCC4', 'VAR')]
data_setA <- data_setA[complete.cases(data_setA), ]

RECQL3_df_A <- data_setA[data_setA$RECQL3 == 1,]
setA <- RECQL3_df_A$VAR


# all setB vars
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "CAgene", 'BUB1B', 'CBL', 'CREBBP', 'EP300', 
                                     'ERCC3', 'FANCA', 'KDM6A', 'KRAS', 'LIG4', 'LZTR1', 
                                     'NRAS', 'POLE', 'PTPN11', 'RECQL3', 'RECQL4', 
                                     'RIT1', 'SBDS', 'SMAD4', 'SMARCAL1', 'SOS1', 
                                     'TALDO1',     'TRIM37', 'XRCC4', 'VAR')]
data_setB <- data_setB[complete.cases(data_setB), ]

RECQL3_df_B <- data_setB[data_setB$RECQL3 == 1,]
setB <- RECQL3_df_B$VAR


# all setD vars
data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "CAgene", 'BUB1B', 'CBL', 'CREBBP', 'EP300', 
                                     'ERCC3', 'FANCA', 'KDM6A', 'LIG4', 'LZTR1', 
                                     'NRAS', 'POLE', 'PTPN11', 'RECQL3', 'RECQL4', 
                                     'RIT1', 'SBDS', 'SMARCAL1', 'SOS1', 'TALDO1', 
                                     'TRIM37', 'XRCC4', 'VAR')]
data_setD <- data_setD[complete.cases(data_setD), ]

RECQL3_df_D <- data_setD[data_setD$RECQL3 == 1,]
setD <- RECQL3_df_D$VAR


# all setE vars
data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "CAgene", 'BUB1B', 'CBL', 
                                     'CREBBP', 'EP300', 'ERCC3', 'FANCA', 'KDM6A', 'LIG4', 
                                     'LZTR1', 'NRAS', 'POLE', 'PTPN11', 'RECQL3', 
                                     'RECQL4', 'RIT1', 'SBDS', 'SMARCAL1', 'SOS1', 
                                     'TALDO1',     'TRIM37', 'XRCC4', 'VAR')]
data_setE <- data_setE[complete.cases(data_setE), ]

RECQL3_df_E <- data_setE[data_setE$RECQL3 == 1,]
setE <- RECQL3_df_E$VAR


# save these as csv files to work with in R on my own computer
write.csv(ALL, "RECQL3_VAR_all.csv", row.names=FALSE)
system('dx upload RECQL3_VAR_all.csv --path Emily-folder/results_2/RECQL3_VAR_all.csv')

write.csv(setA, "RECQL3_VAR_setA.csv", row.names=FALSE)
system('dx upload RECQL3_VAR_setA.csv --path Emily-folder/results_2/RECQL3_VAR_setA.csv')

write.csv(setB, "RECQL3_VAR_setB.csv", row.names=FALSE)
system('dx upload RECQL3_VAR_setB.csv --path Emily-folder/results_2/RECQL3_VAR_setB.csv')

write.csv(setD, "RECQL3_VAR_setD.csv", row.names=FALSE)
system('dx upload RECQL3_VAR_setD.csv --path Emily-folder/results_2/RECQL3_VAR_setD.csv')

write.csv(setE, "RECQL3_VAR_setE.csv", row.names=FALSE)
system('dx upload RECQL3_VAR_setE.csv --path Emily-folder/results_2/RECQL3_VAR_setE.csv')









# RECQL3 15:90784953:C:T_T (pop_var)-----------------------------------------------------



## the most common RECQL3 variant = 15:90784953:C:T_T


## the most common RECQL3 variant = 15:90784953:C:T_T



unique(EntireCohort_pathvar$RECQL3)
# 0 1

RECQL3_df <- EntireCohort_pathvar[EntireCohort_pathvar$RECQL3 == 1,]
# 1575 obs. 


# separate the var names into 9 columns and then use these to make a new column
RECQL3_df <- separate(data = RECQL3_df, col = VAR, into = c("VARIANT1","VARIANT2","VARIANT3","VARIANT4","VARIANT5","VARIANT6","VARIANT7","VARIANT8","VARIANT9"), sep = ";")

# find participants with the 15:90784953:C:T_T variant
VARIANT1_df <- RECQL3_df[RECQL3_df$VARIANT1 == "15:90784953:C:T_T",]
VARIANT1_df$eid
# [1] 3012822 4662767 2650647 3570806 3149709 5860285 4773710 1675903 5058849 3874281 2819054
# [12] 1242919 1407119 5163310 1077634 2858144 1564649 5190584 4655904 3050136 4586382 1316942
# [23] 5610402 5684905 3261395 5480008 4877282 5749406 5860903 5192065 4728146 3387070 4675833
# [34] 5052981 1898254 3592515 5079241 2626723 1357984 4663806 2253890 4956302 5401099 2704287
# [45] 4977754 1255741 5050619 3780141 5053533 1091194 5715976 4086708 4694334 1129463 1859780
# [56] 5971706 1780482 3520999 3326213 1421632 3750393 2848805 2061361 2728005 3380714 1660812
# [67] 1841423 4336860 1604985 1706832 4849724 4111835 4626587 3946486 5567101 2839925 3708276
# [78] 1382118 1138462 2408219 1777028 4792265 2981748 2342771 1742226 4399259 1274137 2059963
# [89] 1523591

VARIANT2_df <- RECQL3_df[RECQL3_df$VARIANT2 == "15:90784953:C:T_T",]
unique(VARIANT2_df$eid)
# [1]      NA 5200356 2631333 5147643 4107968 1495066 1464388 1591403 3783881 4335869 1524171
# [12] 5491324 1176083 1887851 3732444 1508389 2659870 5463261 4961068 2772018 4282269 3658918
# [23] 3881355 3221014 4658222 2287519 5001194 3584891 3840011 2100807 6003836 1284388 1625591
# [34] 3776955 5854076 3511287 2103185 5581054 1865142 4223909 2315569 1685041 4808174 1952717
# [45] 3165794 2339134 3113867 1128801 3839064 5407984 4151644 4473869 3846160 2638501 2279572
# [56] 3793858 3427810 2040030 1004226 5824426 1889480 3459444 1472400 4936221 4296675 1204702
# [67] 4048518 2926212 4783437 2691026 4513883 3822652 4475665 2913284 5591299 3355802 5524352
# [78] 1942200 4721074 3017176 2501038 3824767 1553858

VARIANT3_df <- RECQL3_df[RECQL3_df$VARIANT3 == "15:90784953:C:T_T",]
unique(VARIANT3_df$eid)
# [1]      NA 2906951 5623647 1028020 4657948 5750466 2833511 3036372 2827008 2578703 1114910
# [12] 4525979 5413811 1722461 1979660 4861663 3845716 4991672 3339475 4136619 3959581 1547072
# [23] 5452784 3212218 1230309 3162886 4992095 2575263 1657031 3368070 1014760 3128078 2048876

VARIANT4_df <- RECQL3_df[RECQL3_df$VARIANT4 == "15:90784953:C:T_T",]
unique(VARIANT4_df$eid)
# [1]      NA 2306336 5811111 5715389 4089646 2749708 1819725 4290757 1260350 4163565 3720129


# none in the rest of the VAR columns



pop_var_eid <- c(3012822, 4662767, 2650647, 3570806, 3149709, 5860285, 4773710, 1675903, 5058849, 3874281, 2819054, 
                 1242919, 1407119, 5163310, 1077634, 2858144, 1564649, 5190584, 4655904, 3050136, 4586382, 1316942,
                 5610402, 5684905, 3261395, 5480008, 4877282, 5749406, 5860903, 5192065, 4728146, 3387070, 4675833,
                 5052981, 1898254, 3592515, 5079241, 2626723, 1357984, 4663806, 2253890, 4956302, 5401099, 2704287,
                 4977754, 1255741, 5050619, 3780141, 5053533, 1091194, 5715976, 4086708, 4694334, 1129463, 1859780,
                 5971706, 1780482, 3520999, 3326213, 1421632, 3750393, 2848805, 2061361, 2728005, 3380714, 1660812,
                 1841423, 4336860, 1604985, 1706832, 4849724, 4111835, 4626587, 3946486, 5567101, 2839925, 3708276,
                 1382118, 1138462, 2408219, 1777028, 4792265, 2981748, 2342771, 1742226, 4399259, 1274137, 2059963,
                 1523591, 5200356, 2631333, 5147643, 4107968, 1495066, 1464388, 1591403, 3783881, 4335869, 1524171,
                 5491324, 1176083, 1887851, 3732444, 1508389, 2659870, 5463261, 4961068, 2772018, 4282269, 3658918,
                 3881355, 3221014, 4658222, 2287519, 5001194, 3584891, 3840011, 2100807, 6003836, 1284388, 1625591,
                 3776955, 5854076, 3511287, 2103185, 5581054, 1865142, 4223909, 2315569, 1685041, 4808174, 1952717,
                 3165794, 2339134, 3113867, 1128801, 3839064, 5407984, 4151644, 4473869, 3846160, 2638501, 2279572,
                 3793858, 3427810, 2040030, 1004226, 5824426, 1889480, 3459444, 1472400, 4936221, 4296675, 1204702,
                 4048518, 2926212, 4783437, 2691026, 4513883, 3822652, 4475665, 2913284, 5591299, 3355802, 5524352,
                 1942200, 4721074, 3017176, 2501038, 3824767, 1553858, 2906951, 5623647, 1028020, 4657948, 5750466, 
                 2833511, 3036372, 2827008, 2578703, 1114910,
                 4525979, 5413811, 1722461, 1979660, 4861663, 3845716, 4991672, 3339475, 4136619, 3959581, 1547072,
                 5452784, 3212218, 1230309, 3162886, 4992095, 2575263, 1657031, 3368070, 1014760, 3128078, 2048876,
                 2306336, 5811111, 5715389, 4089646, 2749708, 1819725, 4290757, 1260350, 4163565, 3720129)
# 213


# annotate the EntireCohort_pathvar df with this info
pop_var_df <- data.frame(eid = pop_var_eid,
                         pop_var = 1)

# merge this with the EntireCohort_pathvar df
EntireCohort_pathvar_pop_var <- full_join(EntireCohort_pathvar, pop_var_df, by = 'eid')

EntireCohort_pathvar_pop_var$pop_var[is.na(EntireCohort_pathvar_pop_var$pop_var)] <- 0

unique(EntireCohort_pathvar_pop_var$pop_var)
# 0 1





# incidence:
LR_1_test <- glm(formula = cancer_ROOC_dummy ~ pop_var, data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_1_test)
# Coefficients:
#   Estimate Std. Error  z value Pr(>|z|)    
# (Intercept) -1.174893   0.003437 -341.854   <2e-16 ***
#   pop_var      0.307392   0.150171    2.047   0.0407 *
# exp(cbind(OR = coef(LR_1_test), confint(LR_1_test)))
pathvar_RECQL3_pop_var_incidence_1 <- tidy(LR_1_test, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_incidence_1, "pathvar_RECQL3_pop_var_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_incidence_1.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_incidence_1.csv')


data_setA <- EntireCohort_pathvar_pop_var[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                             "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                             "CAgene", 'BUB1B', 'CBL', 'CREBBP', 
                                             'EP300', 'ERCC3', 'FANCA', 'KDM6A', 'KRAS', 'LIG4', 
                                             'LZTR1', 'NRAS', 'POLE', 'PTPN11', 
                                             'RECQL3', 'RECQL4', 'RIT1', 'SBDS', 'SMAD4', 
                                             'SMARCAL1', 'SOS1', 'TALDO1',     
                                             'TRIM37', 'XRCC4', 'pop_var')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     33

LR_4_test <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, 
                 data = data_setA, family = binomial)
# summary(LR_4_test)
# # p-value = 0.145
# exp(cbind(OR = coef(LR_4_test), confint(LR_4_test)))
pathvar_RECQL3_pop_var_incidence_4 <- tidy(LR_4_test, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_incidence_4, "pathvar_RECQL3_pop_var_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_incidence_4.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_incidence_4.csv')


data_setB <- EntireCohort_pathvar_pop_var[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                             "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                             "BMI", "Standing_height", "University_education",
                                             "Annual_income", "Weekly_exercise",
                                             "CAgene", 'BUB1B', 'CBL', 'CREBBP', 'EP300', 
                                             'ERCC3', 'FANCA', 'KDM6A', 'KRAS', 'LIG4', 'LZTR1', 
                                             'NRAS', 'POLE', 'PTPN11', 'RECQL3', 'RECQL4', 
                                             'RIT1', 'SBDS', 'SMAD4', 'SMARCAL1', 'SOS1', 
                                             'TALDO1',     'TRIM37', 'XRCC4', 'pop_var')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     38

LR_8_test <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+pop_var, 
                 data = data_setB, family = binomial)
# summary(LR_8_test)
# # p = 0.1659
# exp(cbind(OR = coef(LR_8_test), confint(LR_8_test)))
pathvar_RECQL3_pop_var_incidence_8 <- tidy(LR_8_test, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_incidence_8, "pathvar_RECQL3_pop_var_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_incidence_8.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_incidence_8.csv')




# age:
LM_1_test <- lm(Age_at_cancer_diagnosis_earliest ~ pop_var, data = EntireCohort_pathvar_pop_var)
# summary(LM_1_test)
pathvar_RECQL3_pop_var_age_1 <- tidy(LM_1_test, conf.int = TRUE)
write.csv(pathvar_RECQL3_pop_var_age_1, "pathvar_RECQL3_pop_var_age_1.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_age_1.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_age_1.csv')


data_setD <- EntireCohort_pathvar_pop_var[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                             "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                             "CAgene", 'BUB1B', 'CBL', 'CREBBP', 'EP300', 
                                             'ERCC3', 'FANCA', 'KDM6A', 'LIG4', 'LZTR1', 
                                             'NRAS', 'POLE', 'PTPN11', 'RECQL3', 'RECQL4', 
                                             'RIT1', 'SBDS', 'SMARCAL1', 'SOS1', 'TALDO1', 
                                             'TRIM37', 'XRCC4', 'pop_var')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     31

LM_4_test <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, data = data_setD)
# summary(LM_4_test)
# 
pathvar_RECQL3_pop_var_age_4 <- tidy(LM_4_test, conf.int = TRUE)
write.csv(pathvar_RECQL3_pop_var_age_4, "pathvar_RECQL3_pop_var_age_4.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_age_4.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_age_4.csv')

data_setE <- EntireCohort_pathvar_pop_var[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                             "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                             "Standing_height", "University_education", "Annual_income",
                                             "Weekly_exercise", "CAgene", 'BUB1B', 'CBL', 
                                             'CREBBP', 'EP300', 'ERCC3', 'FANCA', 'KDM6A', 'LIG4', 
                                             'LZTR1', 'NRAS', 'POLE', 'PTPN11', 'RECQL3', 
                                             'RECQL4', 'RIT1', 'SBDS', 'SMARCAL1', 'SOS1', 
                                             'TALDO1',     'TRIM37', 'XRCC4', 'pop_var')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    36

LM_16_test <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+pop_var, data = data_setE)
# summary(LM_16_test)
# 
pathvar_RECQL3_pop_var_age_16 <- tidy(LM_16_test, conf.int = TRUE)
write.csv(pathvar_RECQL3_pop_var_age_16, "pathvar_RECQL3_pop_var_age_16.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_age_16.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_age_16.csv')




## type 


data_setA2 <- EntireCohort_pathvar_pop_var[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                              "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                              "CAgene", 'BUB1B', 'CBL', 'CREBBP', 
                                              'EP300', 'ERCC3', 'FANCA', 'KDM6A', 'KRAS', 'LIG4', 
                                              'LZTR1', 'NRAS', 'POLE', 'PTPN11', 
                                              'RECQL3', 'RECQL4', 'RIT1', 'SBDS', 'SMAD4', 
                                              'SMARCAL1', 'SOS1', 'TALDO1',     
                                              'TRIM37', 'XRCC4', 'pop_var', 'ICD_code')]
data_setA2 <- data_setA2[complete.cases(data_setA2), ]
dim(data_setA2)
# 465881     34

# full type analysis for this variant




LR_oral_RECQL3_pop_var <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, 
                              data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_oral_RECQL3_pop_var)
# exp(cbind(OR = coef(LR_oral_RECQL3_pop_var), confint(LR_oral_RECQL3_pop_var)))
pathvar_RECQL3_pop_var_type_1 <- tidy(LR_oral_RECQL3_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_type_1, "pathvar_RECQL3_pop_var_type_1.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_type_1.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_type_1.csv')



LR_digestive_RECQL3_pop_var <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, 
                                   data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_digestive_RECQL3_pop_var)
# exp(cbind(OR = coef(LR_digestive_RECQL3_pop_var), confint(LR_digestive_RECQL3_pop_var)))
pathvar_RECQL3_pop_var_type_2 <- tidy(LR_digestive_RECQL3_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_type_2, "pathvar_RECQL3_pop_var_type_2.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_type_2.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_type_2.csv')



LR_respiratory_RECQL3_pop_var <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, 
                                     data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_respiratory_RECQL3_pop_var)
# exp(cbind(OR = coef(LR_respiratory_RECQL3_pop_var), confint(LR_respiratory_RECQL3_pop_var)))
pathvar_RECQL3_pop_var_type_3 <- tidy(LR_respiratory_RECQL3_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_type_3, "pathvar_RECQL3_pop_var_type_3.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_type_3.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_type_3.csv')



LR_bone_RECQL3_pop_var <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, 
                              data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_bone_RECQL3_pop_var)
# exp(cbind(OR = coef(LR_bone_RECQL3_pop_var), confint(LR_bone_RECQL3_pop_var)))
pathvar_RECQL3_pop_var_type_4 <- tidy(LR_bone_RECQL3_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_type_4, "pathvar_RECQL3_pop_var_type_4.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_type_4.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_type_4.csv')



LR_skin_RECQL3_pop_var <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, 
                              data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_skin_RECQL3_pop_var)
# exp(cbind(OR = coef(LR_skin_RECQL3_pop_var), confint(LR_skin_RECQL3_pop_var)))
pathvar_RECQL3_pop_var_type_5 <- tidy(LR_skin_RECQL3_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_type_5, "pathvar_RECQL3_pop_var_type_5.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_type_5.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_type_5.csv')



LR_mesothelium_RECQL3_pop_var <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, 
                                     data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_mesothelium_RECQL3_pop_var)
# exp(cbind(OR = coef(LR_mesothelium_RECQL3_pop_var), confint(LR_mesothelium_RECQL3_pop_var)))
pathvar_RECQL3_pop_var_type_6 <- tidy(LR_mesothelium_RECQL3_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_type_6, "pathvar_RECQL3_pop_var_type_6.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_type_6.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_type_6.csv')



LR_breast_RECQL3_pop_var <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, 
                                data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_breast_RECQL3_pop_var)
# exp(cbind(OR = coef(LR_breast_RECQL3_pop_var), confint(LR_breast_RECQL3_pop_var)))
pathvar_RECQL3_pop_var_type_7 <- tidy(LR_breast_RECQL3_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_type_7, "pathvar_RECQL3_pop_var_type_7.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_type_7.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_type_7.csv')



LR_femalegenital_RECQL3_pop_var <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, 
                                       data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_femalegenital_RECQL3_pop_var)
# exp(cbind(OR = coef(LR_femalegenital_RECQL3_pop_var), confint(LR_femalegenital_RECQL3_pop_var)))
pathvar_RECQL3_pop_var_type_8 <- tidy(LR_femalegenital_RECQL3_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_type_8, "pathvar_RECQL3_pop_var_type_8.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_type_8.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_type_8.csv')



LR_malegenital_RECQL3_pop_var <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, 
                                     data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_malegenital_RECQL3_pop_var)
# exp(cbind(OR = coef(LR_malegenital_RECQL3_pop_var), confint(LR_malegenital_RECQL3_pop_var)))
pathvar_RECQL3_pop_var_type_9 <- tidy(LR_malegenital_RECQL3_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_type_9, "pathvar_RECQL3_pop_var_type_9.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_type_9.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_type_9.csv')



LR_urinary_RECQL3_pop_var <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, 
                                 data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_urinary_RECQL3_pop_var)
# exp(cbind(OR = coef(LR_urinary_RECQL3_pop_var), confint(LR_urinary_RECQL3_pop_var)))
pathvar_RECQL3_pop_var_type_10 <- tidy(LR_urinary_RECQL3_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_type_10, "pathvar_RECQL3_pop_var_type_10.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_type_10.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_type_10.csv')



LR_cns_RECQL3_pop_var <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, 
                             data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_cns_RECQL3_pop_var)
# exp(cbind(OR = coef(LR_cns_RECQL3_pop_var), confint(LR_cns_RECQL3_pop_var)))
pathvar_RECQL3_pop_var_type_11 <- tidy(LR_cns_RECQL3_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_type_11, "pathvar_RECQL3_pop_var_type_11.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_type_11.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_type_11.csv')



LR_endocrine_RECQL3_pop_var <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, 
                                   data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_endocrine_RECQL3_pop_var)
# exp(cbind(OR = coef(LR_endocrine_RECQL3_pop_var), confint(LR_endocrine_RECQL3_pop_var)))
pathvar_RECQL3_pop_var_type_12 <- tidy(LR_endocrine_RECQL3_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_type_12, "pathvar_RECQL3_pop_var_type_12.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_type_12.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_type_12.csv')



LR_lymphatic_RECQL3_pop_var <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, 
                                   data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_lymphatic_RECQL3_pop_var)
# exp(cbind(OR = coef(LR_lymphatic_RECQL3_pop_var), confint(LR_lymphatic_RECQL3_pop_var)))
pathvar_RECQL3_pop_var_type_13 <- tidy(LR_lymphatic_RECQL3_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_type_13, "pathvar_RECQL3_pop_var_type_13.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_type_13.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_type_13.csv')



LR_secondary_RECQL3_pop_var <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, 
                                   data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_secondary_RECQL3_pop_var)
# exp(cbind(OR = coef(LR_secondary_RECQL3_pop_var), confint(LR_secondary_RECQL3_pop_var)))
pathvar_RECQL3_pop_var_type_14 <- tidy(LR_secondary_RECQL3_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_type_14, "pathvar_RECQL3_pop_var_type_14.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_type_14.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_type_14.csv')



LR_other_RECQL3_pop_var <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var, 
                               data = EntireCohort_pathvar_pop_var, family = binomial)
# summary(LR_other_RECQL3_pop_var)
# exp(cbind(OR = coef(LR_other_RECQL3_pop_var), confint(LR_other_RECQL3_pop_var)))
pathvar_RECQL3_pop_var_type_15 <- tidy(LR_other_RECQL3_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var_type_15, "pathvar_RECQL3_pop_var_type_15.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var_type_15.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var_type_15.csv')












# RECQL3 15:90763016:C:T_T (pop_var2)-----------------------------------------------------



## the 2nd most common RECQL3 variant = 15:90763016:C:T_T





unique(EntireCohort_pathvar$RECQL3)
# 0 1

RECQL3_df <- EntireCohort_pathvar[EntireCohort_pathvar$RECQL3 == 1,]
# 1575 obs. 


# separate the var names into 9 columns and then use these to make a new column
RECQL3_df <- separate(data = RECQL3_df, col = VAR, into = c("VARIANT1","VARIANT2","VARIANT3","VARIANT4","VARIANT5","VARIANT6","VARIANT7","VARIANT8","VARIANT9"), sep = ";")

# find participants with the 15:90763016:C:T_T variant
VARIANT1_df <- RECQL3_df[RECQL3_df$VARIANT1 == "15:90763016:C:T_T",]
VARIANT1_df$eid
# 4522305 1377507 4817908 3995839 5094586 1120730 5938248 1383352 5224694 3116376 2529193 5892110 1831119 2380109 2504779 2430674 5908134 1255796
# 4329871 1671754 1869034 1065529 1650593 5984488 3002523 1040191 1165525 3659715 5162123 4679914 4008966 5422099 5880121 2032901 4409199 4150619
# 4674928 3536123 3598061 1979655 5980326 5739779 1948962 4789048 5384983

VARIANT2_df <- RECQL3_df[RECQL3_df$VARIANT2 == "15:90763016:C:T_T",]
unique(VARIANT2_df$eid)
# 4151661 4318066 3420717 5930623 4487350 4768813 1621984 4005055 3254818 2735731 5183398 5830767 4584780 4823324 5169342 5734539 5365761
# 5315482 5077369 2070083 4788136 4248390 6010494 3742826 3610311 4510219 1820352 5362437 3173665 2178754 2989931 4948763 5303759 4784220 5114937
# 2394427 1450974 4103468 5265025 4925684 5678598 4023680 1087675

VARIANT3_df <- RECQL3_df[RECQL3_df$VARIANT3 == "15:90763016:C:T_T",]
unique(VARIANT3_df$eid)
# 2040088 2130757 1957967 2568138 3676980 1379043 5520986 1128330 3119711 2371284 1251244 3024554 3678653 5456574 4790388 1485544 2081213
# 3608949 2383649 1626605 5634211 1274143 2012674 5602814 4954105 3769352

VARIANT4_df <- RECQL3_df[RECQL3_df$VARIANT4 == "15:90763016:C:T_T",]
unique(VARIANT4_df$eid)
# 1116969 2898293 4361984 5517578 2774910 3671055 2461180 1122543

VARIANT5_df <- RECQL3_df[RECQL3_df$VARIANT5 == "15:90763016:C:T_T",]
unique(VARIANT5_df$eid)
# 3635764

VARIANT6_df <- RECQL3_df[RECQL3_df$VARIANT6 == "15:90763016:C:T_T",]
unique(VARIANT6_df$eid)
# 2162821


# none in the rest of the VAR columns



pop_var2_eid <- c(4522305, 1377507, 4817908, 3995839, 5094586, 1120730, 5938248, 1383352, 5224694, 
                  3116376, 2529193, 5892110, 1831119, 2380109, 2504779, 2430674, 5908134, 1255796, 
                  4329871, 1671754, 1869034, 1065529, 1650593, 5984488, 3002523, 1040191, 1165525, 
                  3659715, 5162123, 4679914, 4008966, 5422099, 5880121, 2032901, 4409199, 4150619, 
                  4674928, 3536123, 3598061, 1979655, 5980326, 5739779, 1948962, 4789048, 5384983, 
                  4151661, 4318066, 3420717, 5930623, 4487350, 4768813, 1621984, 4005055, 3254818, 
                  2735731, 5183398, 5830767, 4584780, 4823324, 5169342, 5734539, 5365761, 5315482, 
                  5077369, 2070083, 4788136, 4248390, 6010494, 3742826, 3610311, 4510219, 1820352, 
                  5362437, 3173665, 2178754, 2989931, 4948763, 5303759, 4784220, 5114937, 2394427, 
                  1450974, 4103468, 5265025, 4925684, 5678598, 4023680, 1087675, 2040088, 2130757, 
                  1957967, 2568138, 3676980, 1379043, 5520986, 1128330, 3119711, 2371284, 1251244, 
                  3024554, 3678653, 5456574, 4790388, 1485544, 2081213, 3608949, 2383649, 1626605, 
                  5634211, 1274143, 2012674, 5602814, 4954105, 3769352, 1116969, 2898293, 4361984, 
                  5517578, 2774910, 3671055, 2461180, 1122543, 3635764, 2162821)
# 124 (matches :))


# annotate the EntireCohort_pathvar df with this info
pop_var2_df <- data.frame(eid = pop_var2_eid,
                          pop_var2 = 1)

# merge this with the EntireCohort_pathvar df
EntireCohort_pathvar_pop_var2 <- full_join(EntireCohort_pathvar, pop_var2_df, by = 'eid')

EntireCohort_pathvar_pop_var2$pop_var2[is.na(EntireCohort_pathvar_pop_var2$pop_var2)] <- 0

unique(EntireCohort_pathvar_pop_var2$pop_var2)
# 0 1





# incidence:
LR_1_test <- glm(formula = cancer_ROOC_dummy ~ pop_var2, data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_1_test)
# Coefficients:
# Estimate Std. Error  z value Pr(>|z|)    
# (Intercept) -1.174917   0.003437 -341.891  < 2e-16 ***
#   pop_var2     0.577080   0.187721    3.074  0.00211 ** 
pathvar_RECQL3_pop_var2_incidence_1 <- tidy(LR_1_test, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_incidence_1, "pathvar_RECQL3_pop_var2_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_incidence_1.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_incidence_1.csv')


data_setA <- EntireCohort_pathvar_pop_var2[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                              "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                              "CAgene", 'BUB1B', 'CBL', 'CREBBP', 
                                              'EP300', 'ERCC3', 'FANCA', 'KDM6A', 'KRAS', 'LIG4', 
                                              'LZTR1', 'NRAS', 'POLE', 'PTPN11', 
                                              'RECQL3', 'RECQL4', 'RIT1', 'SBDS', 'SMAD4', 
                                              'SMARCAL1', 'SOS1', 'TALDO1',     
                                              'TRIM37', 'XRCC4', 'pop_var2')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     33

LR_4_test <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, 
                 data = data_setA, family = binomial)
# summary(LR_4_test)
# p = 0.00265 ** 
# exp(cbind(OR = coef(LR_4_test), confint(LR_4_test)))
pathvar_RECQL3_pop_var2_incidence_4 <- tidy(LR_4_test, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_incidence_4, "pathvar_RECQL3_pop_var2_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_incidence_4.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_incidence_4.csv')


data_setB <- EntireCohort_pathvar_pop_var2[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                              "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                              "BMI", "Standing_height", "University_education",
                                              "Annual_income", "Weekly_exercise",
                                              "CAgene", 'BUB1B', 'CBL', 'CREBBP', 'EP300', 
                                              'ERCC3', 'FANCA', 'KDM6A', 'KRAS', 'LIG4', 'LZTR1', 
                                              'NRAS', 'POLE', 'PTPN11', 'RECQL3', 'RECQL4', 
                                              'RIT1', 'SBDS', 'SMAD4', 'SMARCAL1', 'SOS1', 
                                              'TALDO1',     'TRIM37', 'XRCC4', 'pop_var2')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     38

LR_8_test <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+pop_var2, 
                 data = data_setB, family = binomial)
# summary(LR_8_test)
# p=0.0109 *
# exp(cbind(OR = coef(LR_8_test), confint(LR_8_test)))
pathvar_RECQL3_pop_var2_incidence_8 <- tidy(LR_8_test, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_incidence_8, "pathvar_RECQL3_pop_var2_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_incidence_8.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_incidence_8.csv')




# age:
LM_1_test <- lm(Age_at_cancer_diagnosis_earliest ~ pop_var2, data = EntireCohort_pathvar_pop_var2)
# summary(LM_1_test)
#   
pathvar_RECQL3_pop_var2_age_1 <- tidy(LM_1_test, conf.int = TRUE)
write.csv(pathvar_RECQL3_pop_var2_age_1, "pathvar_RECQL3_pop_var2_age_1.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_age_1.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_age_1.csv')


data_setD <- EntireCohort_pathvar_pop_var2[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                              "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                              "CAgene", 'BUB1B', 'CBL', 'CREBBP', 'EP300', 
                                              'ERCC3', 'FANCA', 'KDM6A', 'LIG4', 'LZTR1', 
                                              'NRAS', 'POLE', 'PTPN11', 'RECQL3', 'RECQL4', 
                                              'RIT1', 'SBDS', 'SMARCAL1', 'SOS1', 'TALDO1', 
                                              'TRIM37', 'XRCC4', 'pop_var2')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     31

LM_4_test <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, data = data_setD)
# summary(LM_4_test)
#    
pathvar_RECQL3_pop_var2_age_4 <- tidy(LM_4_test, conf.int = TRUE)
write.csv(pathvar_RECQL3_pop_var2_age_4, "pathvar_RECQL3_pop_var2_age_4.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_age_4.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_age_4.csv')

data_setE <- EntireCohort_pathvar_pop_var2[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                              "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                              "Standing_height", "University_education", "Annual_income",
                                              "Weekly_exercise", "CAgene", 'BUB1B', 'CBL', 
                                              'CREBBP', 'EP300', 'ERCC3', 'FANCA', 'KDM6A', 'LIG4', 
                                              'LZTR1', 'NRAS', 'POLE', 'PTPN11', 'RECQL3', 
                                              'RECQL4', 'RIT1', 'SBDS', 'SMARCAL1', 'SOS1', 
                                              'TALDO1',     'TRIM37', 'XRCC4', 'pop_var2')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    36

LM_16_test <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+pop_var2, data = data_setE)
# summary(LM_16_test)
#     
pathvar_RECQL3_pop_var2_age_16 <- tidy(LM_16_test, conf.int = TRUE)
write.csv(pathvar_RECQL3_pop_var2_age_16, "pathvar_RECQL3_pop_var2_age_16.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_age_16.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_age_16.csv')




## type 


data_setA2 <- EntireCohort_pathvar_pop_var2[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                               "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                               "CAgene", 'BUB1B', 'CBL', 'CREBBP', 
                                               'EP300', 'ERCC3', 'FANCA', 'KDM6A', 'KRAS', 'LIG4', 
                                               'LZTR1', 'NRAS', 'POLE', 'PTPN11', 
                                               'RECQL3', 'RECQL4', 'RIT1', 'SBDS', 'SMAD4', 
                                               'SMARCAL1', 'SOS1', 'TALDO1',     
                                               'TRIM37', 'XRCC4', 'pop_var2', 'ICD_code')]
data_setA2 <- data_setA2[complete.cases(data_setA2), ]
dim(data_setA2)
# 465881     34

# full type analysis for this variant




LR_oral_RECQL3_pop_var2 <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, 
                               data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_oral_RECQL3_pop_var2)
# exp(cbind(OR = coef(LR_oral_RECQL3_pop_var2), confint(LR_oral_RECQL3_pop_var2)))
pathvar_RECQL3_pop_var2_type_1 <- tidy(LR_oral_RECQL3_pop_var2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_type_1, "pathvar_RECQL3_pop_var2_type_1.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_type_1.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_type_1.csv')



LR_digestive_RECQL3_pop_var2 <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, 
                                    data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_digestive_RECQL3_pop_var2)
# exp(cbind(OR = coef(LR_digestive_RECQL3_pop_var2), confint(LR_digestive_RECQL3_pop_var2)))
pathvar_RECQL3_pop_var2_type_2 <- tidy(LR_digestive_RECQL3_pop_var2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_type_2, "pathvar_RECQL3_pop_var2_type_2.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_type_2.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_type_2.csv')



LR_respiratory_RECQL3_pop_var2 <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, 
                                      data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_respiratory_RECQL3_pop_var2)
# exp(cbind(OR = coef(LR_respiratory_RECQL3_pop_var2), confint(LR_respiratory_RECQL3_pop_var2)))
pathvar_RECQL3_pop_var2_type_3 <- tidy(LR_respiratory_RECQL3_pop_var2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_type_3, "pathvar_RECQL3_pop_var2_type_3.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_type_3.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_type_3.csv')



LR_bone_RECQL3_pop_var2 <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, 
                               data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_bone_RECQL3_pop_var2)
# exp(cbind(OR = coef(LR_bone_RECQL3_pop_var2), confint(LR_bone_RECQL3_pop_var2)))
pathvar_RECQL3_pop_var2_type_4 <- tidy(LR_bone_RECQL3_pop_var2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_type_4, "pathvar_RECQL3_pop_var2_type_4.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_type_4.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_type_4.csv')



LR_skin_RECQL3_pop_var2 <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, 
                               data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_skin_RECQL3_pop_var2)
# exp(cbind(OR = coef(LR_skin_RECQL3_pop_var2), confint(LR_skin_RECQL3_pop_var2)))
pathvar_RECQL3_pop_var2_type_5 <- tidy(LR_skin_RECQL3_pop_var2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_type_5, "pathvar_RECQL3_pop_var2_type_5.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_type_5.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_type_5.csv')



LR_mesothelium_RECQL3_pop_var2 <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, 
                                      data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_mesothelium_RECQL3_pop_var2)
# exp(cbind(OR = coef(LR_mesothelium_RECQL3_pop_var2), confint(LR_mesothelium_RECQL3_pop_var2)))
pathvar_RECQL3_pop_var2_type_6 <- tidy(LR_mesothelium_RECQL3_pop_var2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_type_6, "pathvar_RECQL3_pop_var2_type_6.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_type_6.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_type_6.csv')



LR_breast_RECQL3_pop_var2 <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, 
                                 data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_breast_RECQL3_pop_var2)
# exp(cbind(OR = coef(LR_breast_RECQL3_pop_var2), confint(LR_breast_RECQL3_pop_var2)))
pathvar_RECQL3_pop_var2_type_7 <- tidy(LR_breast_RECQL3_pop_var2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_type_7, "pathvar_RECQL3_pop_var2_type_7.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_type_7.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_type_7.csv')



LR_femalegenital_RECQL3_pop_var2 <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, 
                                        data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_femalegenital_RECQL3_pop_var2)
# exp(cbind(OR = coef(LR_femalegenital_RECQL3_pop_var2), confint(LR_femalegenital_RECQL3_pop_var2)))
pathvar_RECQL3_pop_var2_type_8 <- tidy(LR_femalegenital_RECQL3_pop_var2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_type_8, "pathvar_RECQL3_pop_var2_type_8.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_type_8.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_type_8.csv')



LR_malegenital_RECQL3_pop_var2 <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, 
                                      data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_malegenital_RECQL3_pop_var2)
# exp(cbind(OR = coef(LR_malegenital_RECQL3_pop_var2), confint(LR_malegenital_RECQL3_pop_var2)))
pathvar_RECQL3_pop_var2_type_9 <- tidy(LR_malegenital_RECQL3_pop_var2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_type_9, "pathvar_RECQL3_pop_var2_type_9.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_type_9.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_type_9.csv')



LR_urinary_RECQL3_pop_var2 <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, 
                                  data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_urinary_RECQL3_pop_var2)
# exp(cbind(OR = coef(LR_urinary_RECQL3_pop_var2), confint(LR_urinary_RECQL3_pop_var2)))
pathvar_RECQL3_pop_var2_type_10 <- tidy(LR_urinary_RECQL3_pop_var2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_type_10, "pathvar_RECQL3_pop_var2_type_10.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_type_10.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_type_10.csv')



LR_cns_RECQL3_pop_var2 <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, 
                              data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_cns_RECQL3_pop_var2)
# exp(cbind(OR = coef(LR_cns_RECQL3_pop_var2), confint(LR_cns_RECQL3_pop_var2)))
pathvar_RECQL3_pop_var2_type_11 <- tidy(LR_cns_RECQL3_pop_var2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_type_11, "pathvar_RECQL3_pop_var2_type_11.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_type_11.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_type_11.csv')



LR_endocrine_RECQL3_pop_var2 <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, 
                                    data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_endocrine_RECQL3_pop_var2)
# exp(cbind(OR = coef(LR_endocrine_RECQL3_pop_var2), confint(LR_endocrine_RECQL3_pop_var2)))
pathvar_RECQL3_pop_var2_type_12 <- tidy(LR_endocrine_RECQL3_pop_var2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_type_12, "pathvar_RECQL3_pop_var2_type_12.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_type_12.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_type_12.csv')



LR_lymphatic_RECQL3_pop_var2 <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, 
                                    data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_lymphatic_RECQL3_pop_var2)
# exp(cbind(OR = coef(LR_lymphatic_RECQL3_pop_var2), confint(LR_lymphatic_RECQL3_pop_var2)))
pathvar_RECQL3_pop_var2_type_13 <- tidy(LR_lymphatic_RECQL3_pop_var2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_type_13, "pathvar_RECQL3_pop_var2_type_13.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_type_13.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_type_13.csv')



LR_secondary_RECQL3_pop_var2 <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, 
                                    data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_secondary_RECQL3_pop_var2)
# exp(cbind(OR = coef(LR_secondary_RECQL3_pop_var2), confint(LR_secondary_RECQL3_pop_var2)))
pathvar_RECQL3_pop_var2_type_14 <- tidy(LR_secondary_RECQL3_pop_var2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_type_14, "pathvar_RECQL3_pop_var2_type_14.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_type_14.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_type_14.csv')



LR_other_RECQL3_pop_var2 <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+pop_var2, 
                                data = EntireCohort_pathvar_pop_var2, family = binomial)
# summary(LR_other_RECQL3_pop_var2)
# exp(cbind(OR = coef(LR_other_RECQL3_pop_var2), confint(LR_other_RECQL3_pop_var2)))
pathvar_RECQL3_pop_var2_type_15 <- tidy(LR_other_RECQL3_pop_var2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_pop_var2_type_15, "pathvar_RECQL3_pop_var2_type_15.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_pop_var2_type_15.csv --path Emily-folder/results_2/pathvar_RECQL3_pop_var2_type_15.csv')






# RECQL3 not-pop-var ------------------------------------------------------



three_pop_vars_eid <- c(2520015, 3637195, 4839300, 5937961, 2690742, 2740419, 1154132, 1422480, 3477890, 2364041, 3247763, 
                        3832791, 2281401, 2729683, 5896647, 1077655, 4021672, 3662351, 5637750, 3158925, 2615302, 1068710, 
                        3095718, 4100571, 1161396, 5312744, 1504704, 5933176, 3125046, 3770061, 4435376, 2319402, 5248709, 
                        1895330, 3287818, 4008709, 3122774, 2182357, 1792079, 2998057, 1937090, 4083551, 2088360, 5932824,               
                        2684750, 4931468, 3422053, 1538819, 2995926, 2769336, 3718607, 4626623, 4565898, 2987048, 1492191,
                        2482945, 4346259, 3286283, 1278579, 2805615, 1005157, 1586758, 2821282, 4230882,
                        1697034, 3971758, 5247705, 4306266, 3308631, 1523633, 4768106, 4483533, 2708719, 3194129,
                        4406687, 1000686, 3416505, 1790019, 4084183, 3759189, 2107105, 3552575, 3467572,
                        2267832, 2822347, 5649697, 4696204, 2702585, 5860773, 3012822, 4662767, 2650647, 3570806, 3149709, 
                        5860285, 4773710, 1675903, 5058849, 3874281, 2819054, 
                        1242919, 1407119, 5163310, 1077634, 2858144, 1564649, 5190584, 4655904, 3050136, 4586382, 1316942,
                        5610402, 5684905, 3261395, 5480008, 4877282, 5749406, 5860903, 5192065, 4728146, 3387070, 4675833,
                        5052981, 1898254, 3592515, 5079241, 2626723, 1357984, 4663806, 2253890, 4956302, 5401099, 2704287,
                        4977754, 1255741, 5050619, 3780141, 5053533, 1091194, 5715976, 4086708, 4694334, 1129463, 1859780,
                        5971706, 1780482, 3520999, 3326213, 1421632, 3750393, 2848805, 2061361, 2728005, 3380714, 1660812,
                        1841423, 4336860, 1604985, 1706832, 4849724, 4111835, 4626587, 3946486, 5567101, 2839925, 3708276,
                        1382118, 1138462, 2408219, 1777028, 4792265, 2981748, 2342771, 1742226, 4399259, 1274137, 2059963,
                        1523591, 5200356, 2631333, 5147643, 4107968, 1495066, 1464388, 1591403, 3783881, 4335869, 1524171,
                        5491324, 1176083, 1887851, 3732444, 1508389, 2659870, 5463261, 4961068, 2772018, 4282269, 3658918,
                        3881355, 3221014, 4658222, 2287519, 5001194, 3584891, 3840011, 2100807, 6003836, 1284388, 1625591,
                        3776955, 5854076, 3511287, 2103185, 5581054, 1865142, 4223909, 2315569, 1685041, 4808174, 1952717,
                        3165794, 2339134, 3113867, 1128801, 3839064, 5407984, 4151644, 4473869, 3846160, 2638501, 2279572,
                        3793858, 3427810, 2040030, 1004226, 5824426, 1889480, 3459444, 1472400, 4936221, 4296675, 1204702,
                        4048518, 2926212, 4783437, 2691026, 4513883, 3822652, 4475665, 2913284, 5591299, 3355802, 5524352,
                        1942200, 4721074, 3017176, 2501038, 3824767, 1553858, 2906951, 5623647, 1028020, 4657948, 5750466, 
                        2833511, 3036372, 2827008, 2578703, 1114910,
                        4525979, 5413811, 1722461, 1979660, 4861663, 3845716, 4991672, 3339475, 4136619, 3959581, 1547072,
                        5452784, 3212218, 1230309, 3162886, 4992095, 2575263, 1657031, 3368070, 1014760, 3128078, 2048876,
                        2306336, 5811111, 5715389, 4089646, 2749708, 1819725, 4290757, 1260350, 4163565, 3720129, 
                        4522305, 1377507, 4817908, 3995839, 5094586, 1120730, 5938248, 1383352, 5224694, 
                        3116376, 2529193, 5892110, 1831119, 2380109, 2504779, 2430674, 5908134, 1255796, 
                        4329871, 1671754, 1869034, 1065529, 1650593, 5984488, 3002523, 1040191, 1165525, 
                        3659715, 5162123, 4679914, 4008966, 5422099, 5880121, 2032901, 4409199, 4150619, 
                        4674928, 3536123, 3598061, 1979655, 5980326, 5739779, 1948962, 4789048, 5384983, 
                        4151661, 4318066, 3420717, 5930623, 4487350, 4768813, 1621984, 4005055, 3254818, 
                        2735731, 5183398, 5830767, 4584780, 4823324, 5169342, 5734539, 5365761, 5315482, 
                        5077369, 2070083, 4788136, 4248390, 6010494, 3742826, 3610311, 4510219, 1820352, 
                        5362437, 3173665, 2178754, 2989931, 4948763, 5303759, 4784220, 5114937, 2394427, 
                        1450974, 4103468, 5265025, 4925684, 5678598, 4023680, 1087675, 2040088, 2130757, 
                        1957967, 2568138, 3676980, 1379043, 5520986, 1128330, 3119711, 2371284, 1251244, 
                        3024554, 3678653, 5456574, 4790388, 1485544, 2081213, 3608949, 2383649, 1626605, 
                        5634211, 1274143, 2012674, 5602814, 4954105, 3769352, 1116969, 2898293, 4361984, 
                        5517578, 2774910, 3671055, 2461180, 1122543, 3635764, 2162821)





# annotate the EntireCohort_pathvar df with this info
not_pop_var_df <- data.frame(eid = three_pop_vars_eid,
                             not_pop_var = 0)

# merge this with the EntireCohort_pathvar df
EntireCohort_pathvar_not_pop_var <- full_join(EntireCohort_pathvar, not_pop_var_df, by = 'eid')

EntireCohort_pathvar_not_pop_var$not_pop_var[is.na(EntireCohort_pathvar_not_pop_var$not_pop_var)] <- 1

# NA = 1
# pop_var = 0

EntireCohort_pathvar_not_pop_var$not_pop_var[EntireCohort_pathvar_not_pop_var$RECQL3 == 0 & EntireCohort_pathvar_not_pop_var$not_pop_var != 0] <- 0

unique(EntireCohort_pathvar_not_pop_var$not_pop_var)
# 1 0

with(EntireCohort_pathvar_not_pop_var, table(RECQL3, not_pop_var))
#          not_pop_var
# RECQL3      0      1
# 0      468224      0
# 1         426   1149











# incidence:
LR_1_test <- glm(formula = cancer_ROOC_dummy ~ not_pop_var, data = EntireCohort_pathvar_not_pop_var, family = binomial)
# summary(LR_1_test)
# 0.0975 .
pathvar_RECQL3_not_pop_var_incidence_1 <- tidy(LR_1_test, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_incidence_1, "pathvar_RECQL3_not_pop_var_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_incidence_1.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_incidence_1.csv')


data_setA <- EntireCohort_pathvar_not_pop_var[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                                 "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                                 "CAgene", 'BUB1B', 'CBL', 'CREBBP', 
                                                 'EP300', 'ERCC3', 'FANCA', 'KDM6A', 'KRAS', 'LIG4', 
                                                 'LZTR1', 'NRAS', 'POLE', 'PTPN11', 
                                                 'RECQL3', 'RECQL4', 'RIT1', 'SBDS', 'SMAD4', 
                                                 'SMARCAL1', 'SOS1', 'TALDO1',     
                                                 'TRIM37', 'XRCC4', 'not_pop_var')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     33

LR_4_test <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, 
                 data = data_setA, family = binomial)
# summary(LR_4_test)
#  0.024 *
# exp(cbind(OR = coef(LR_4_test), confint(LR_4_test)))
pathvar_RECQL3_not_pop_var_incidence_4 <- tidy(LR_4_test, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_incidence_4, "pathvar_RECQL3_not_pop_var_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_incidence_4.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_incidence_4.csv')


data_setB <- EntireCohort_pathvar_not_pop_var[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                                 "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                                 "BMI", "Standing_height", "University_education",
                                                 "Annual_income", "Weekly_exercise",
                                                 "CAgene", 'BUB1B', 'CBL', 'CREBBP', 'EP300', 
                                                 'ERCC3', 'FANCA', 'KDM6A', 'KRAS', 'LIG4', 'LZTR1', 
                                                 'NRAS', 'POLE', 'PTPN11', 'RECQL3', 'RECQL4', 
                                                 'RIT1', 'SBDS', 'SMAD4', 'SMARCAL1', 'SOS1', 
                                                 'TALDO1',     'TRIM37', 'XRCC4', 'not_pop_var')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     38

LR_8_test <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+not_pop_var, 
                 data = data_setB, family = binomial)
# summary(LR_8_test)
# 0.00834 **
# exp(cbind(OR = coef(LR_8_test), confint(LR_8_test)))
pathvar_RECQL3_not_pop_var_incidence_8 <- tidy(LR_8_test, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_incidence_8, "pathvar_RECQL3_not_pop_var_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_incidence_8.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_incidence_8.csv')




# age:
LM_1_test <- lm(Age_at_cancer_diagnosis_earliest ~ not_pop_var, data = EntireCohort_pathvar_not_pop_var)
# summary(LM_1_test)
pathvar_RECQL3_not_pop_var_age_1 <- tidy(LM_1_test, conf.int = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_age_1, "pathvar_RECQL3_not_pop_var_age_1.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_age_1.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_age_1.csv')


data_setD <- EntireCohort_pathvar_not_pop_var[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                                 "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                                 "CAgene", 'BUB1B', 'CBL', 'CREBBP', 'EP300', 
                                                 'ERCC3', 'FANCA', 'KDM6A', 'LIG4', 'LZTR1', 
                                                 'NRAS', 'POLE', 'PTPN11', 'RECQL3', 'RECQL4', 
                                                 'RIT1', 'SBDS', 'SMARCAL1', 'SOS1', 'TALDO1', 
                                                 'TRIM37', 'XRCC4', 'not_pop_var')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     31

LM_4_test <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, data = data_setD)
# summary(LM_4_test)
#     
pathvar_RECQL3_not_pop_var_age_4 <- tidy(LM_4_test, conf.int = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_age_4, "pathvar_RECQL3_not_pop_var_age_4.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_age_4.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_age_4.csv')

data_setE <- EntireCohort_pathvar_not_pop_var[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                                 "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                                 "Standing_height", "University_education", "Annual_income",
                                                 "Weekly_exercise", "CAgene", 'BUB1B', 'CBL', 
                                                 'CREBBP', 'EP300', 'ERCC3', 'FANCA', 'KDM6A', 'LIG4', 
                                                 'LZTR1', 'NRAS', 'POLE', 'PTPN11', 'RECQL3', 
                                                 'RECQL4', 'RIT1', 'SBDS', 'SMARCAL1', 'SOS1', 
                                                 'TALDO1',     'TRIM37', 'XRCC4', 'not_pop_var')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    36

LM_16_test <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+not_pop_var, data = data_setE)
# summary(LM_16_test)
#     
pathvar_RECQL3_not_pop_var_age_16 <- tidy(LM_16_test, conf.int = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_age_16, "pathvar_RECQL3_not_pop_var_age_16.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_age_16.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_age_16.csv')




## type 


data_setA2 <- EntireCohort_pathvar_not_pop_var[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                                  "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                                  "CAgene", 'BUB1B', 'CBL', 'CREBBP', 
                                                  'EP300', 'ERCC3', 'FANCA', 'KDM6A', 'KRAS', 'LIG4', 
                                                  'LZTR1', 'NRAS', 'POLE', 'PTPN11', 
                                                  'RECQL3', 'RECQL4', 'RIT1', 'SBDS', 'SMAD4', 
                                                  'SMARCAL1', 'SOS1', 'TALDO1',     
                                                  'TRIM37', 'XRCC4', 'not_pop_var', 'ICD_code')]
data_setA2 <- data_setA2[complete.cases(data_setA2), ]
dim(data_setA2)
# 465881     34


# full type analysis for this variant




LR_oral_RECQL3_not_pop_var <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, 
                                  data = EntireCohort_pathvar_not_pop_var, family = binomial)
# summary(LR_oral_RECQL3_not_pop_var)
# exp(cbind(OR = coef(LR_oral_RECQL3_not_pop_var), confint(LR_oral_RECQL3_not_pop_var)))
pathvar_RECQL3_not_pop_var_type_1 <- tidy(LR_oral_RECQL3_not_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_type_1, "pathvar_RECQL3_not_pop_var_type_1.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_type_1.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_type_1.csv')



LR_digestive_RECQL3_not_pop_var <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, 
                                       data = EntireCohort_pathvar_not_pop_var, family = binomial)
# summary(LR_digestive_RECQL3_not_pop_var)
# exp(cbind(OR = coef(LR_digestive_RECQL3_not_pop_var), confint(LR_digestive_RECQL3_not_pop_var)))
pathvar_RECQL3_not_pop_var_type_2 <- tidy(LR_digestive_RECQL3_not_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_type_2, "pathvar_RECQL3_not_pop_var_type_2.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_type_2.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_type_2.csv')



LR_respiratory_RECQL3_not_pop_var <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, 
                                         data = EntireCohort_pathvar_not_pop_var, family = binomial)
# summary(LR_respiratory_RECQL3_not_pop_var)
# exp(cbind(OR = coef(LR_respiratory_RECQL3_not_pop_var), confint(LR_respiratory_RECQL3_not_pop_var)))
pathvar_RECQL3_not_pop_var_type_3 <- tidy(LR_respiratory_RECQL3_not_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_type_3, "pathvar_RECQL3_not_pop_var_type_3.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_type_3.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_type_3.csv')



LR_bone_RECQL3_not_pop_var <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, 
                                  data = EntireCohort_pathvar_not_pop_var, family = binomial)
# summary(LR_bone_RECQL3_not_pop_var)
# exp(cbind(OR = coef(LR_bone_RECQL3_not_pop_var), confint(LR_bone_RECQL3_not_pop_var)))
pathvar_RECQL3_not_pop_var_type_4 <- tidy(LR_bone_RECQL3_not_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_type_4, "pathvar_RECQL3_not_pop_var_type_4.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_type_4.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_type_4.csv')



LR_skin_RECQL3_not_pop_var <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, 
                                  data = EntireCohort_pathvar_not_pop_var, family = binomial)
# summary(LR_skin_RECQL3_not_pop_var)
# exp(cbind(OR = coef(LR_skin_RECQL3_not_pop_var), confint(LR_skin_RECQL3_not_pop_var)))
pathvar_RECQL3_not_pop_var_type_5 <- tidy(LR_skin_RECQL3_not_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_type_5, "pathvar_RECQL3_not_pop_var_type_5.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_type_5.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_type_5.csv')



LR_mesothelium_RECQL3_not_pop_var <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, 
                                         data = EntireCohort_pathvar_not_pop_var, family = binomial)
# summary(LR_mesothelium_RECQL3_not_pop_var)
# exp(cbind(OR = coef(LR_mesothelium_RECQL3_not_pop_var), confint(LR_mesothelium_RECQL3_not_pop_var)))
pathvar_RECQL3_not_pop_var_type_6 <- tidy(LR_mesothelium_RECQL3_not_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_type_6, "pathvar_RECQL3_not_pop_var_type_6.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_type_6.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_type_6.csv')



LR_breast_RECQL3_not_pop_var <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, 
                                    data = EntireCohort_pathvar_not_pop_var, family = binomial)
# summary(LR_breast_RECQL3_not_pop_var)
# exp(cbind(OR = coef(LR_breast_RECQL3_not_pop_var), confint(LR_breast_RECQL3_not_pop_var)))
pathvar_RECQL3_not_pop_var_type_7 <- tidy(LR_breast_RECQL3_not_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_type_7, "pathvar_RECQL3_not_pop_var_type_7.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_type_7.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_type_7.csv')



LR_femalegenital_RECQL3_not_pop_var <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, 
                                           data = EntireCohort_pathvar_not_pop_var, family = binomial)
# summary(LR_femalegenital_RECQL3_not_pop_var)
# exp(cbind(OR = coef(LR_femalegenital_RECQL3_not_pop_var), confint(LR_femalegenital_RECQL3_not_pop_var)))
pathvar_RECQL3_not_pop_var_type_8 <- tidy(LR_femalegenital_RECQL3_not_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_type_8, "pathvar_RECQL3_not_pop_var_type_8.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_type_8.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_type_8.csv')



LR_malegenital_RECQL3_not_pop_var <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, 
                                         data = EntireCohort_pathvar_not_pop_var, family = binomial)
# summary(LR_malegenital_RECQL3_not_pop_var)
# exp(cbind(OR = coef(LR_malegenital_RECQL3_not_pop_var), confint(LR_malegenital_RECQL3_not_pop_var)))
pathvar_RECQL3_not_pop_var_type_9 <- tidy(LR_malegenital_RECQL3_not_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_type_9, "pathvar_RECQL3_not_pop_var_type_9.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_type_9.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_type_9.csv')



LR_urinary_RECQL3_not_pop_var <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, 
                                     data = EntireCohort_pathvar_not_pop_var, family = binomial)
# summary(LR_urinary_RECQL3_not_pop_var)
# exp(cbind(OR = coef(LR_urinary_RECQL3_not_pop_var), confint(LR_urinary_RECQL3_not_pop_var)))
pathvar_RECQL3_not_pop_var_type_10 <- tidy(LR_urinary_RECQL3_not_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_type_10, "pathvar_RECQL3_not_pop_var_type_10.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_type_10.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_type_10.csv')



LR_cns_RECQL3_not_pop_var <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, 
                                 data = EntireCohort_pathvar_not_pop_var, family = binomial)
# summary(LR_cns_RECQL3_not_pop_var)
# exp(cbind(OR = coef(LR_cns_RECQL3_not_pop_var), confint(LR_cns_RECQL3_not_pop_var)))
pathvar_RECQL3_not_pop_var_type_11 <- tidy(LR_cns_RECQL3_not_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_type_11, "pathvar_RECQL3_not_pop_var_type_11.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_type_11.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_type_11.csv')



LR_endocrine_RECQL3_not_pop_var <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, 
                                       data = EntireCohort_pathvar_not_pop_var, family = binomial)
# summary(LR_endocrine_RECQL3_not_pop_var)
# exp(cbind(OR = coef(LR_endocrine_RECQL3_not_pop_var), confint(LR_endocrine_RECQL3_not_pop_var)))
pathvar_RECQL3_not_pop_var_type_12 <- tidy(LR_endocrine_RECQL3_not_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_type_12, "pathvar_RECQL3_not_pop_var_type_12.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_type_12.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_type_12.csv')



LR_lymphatic_RECQL3_not_pop_var <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, 
                                       data = EntireCohort_pathvar_not_pop_var, family = binomial)
# summary(LR_lymphatic_RECQL3_not_pop_var)
# exp(cbind(OR = coef(LR_lymphatic_RECQL3_not_pop_var), confint(LR_lymphatic_RECQL3_not_pop_var)))
pathvar_RECQL3_not_pop_var_type_13 <- tidy(LR_lymphatic_RECQL3_not_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_type_13, "pathvar_RECQL3_not_pop_var_type_13.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_type_13.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_type_13.csv')



LR_secondary_RECQL3_not_pop_var <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, 
                                       data = EntireCohort_pathvar_not_pop_var, family = binomial)
# summary(LR_secondary_RECQL3_not_pop_var)
# exp(cbind(OR = coef(LR_secondary_RECQL3_not_pop_var), confint(LR_secondary_RECQL3_not_pop_var)))
pathvar_RECQL3_not_pop_var_type_14 <- tidy(LR_secondary_RECQL3_not_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_type_14, "pathvar_RECQL3_not_pop_var_type_14.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_type_14.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_type_14.csv')



LR_other_RECQL3_not_pop_var <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+not_pop_var, 
                                   data = EntireCohort_pathvar_not_pop_var, family = binomial)
# summary(LR_other_RECQL3_not_pop_var)
# exp(cbind(OR = coef(LR_other_RECQL3_not_pop_var), confint(LR_other_RECQL3_not_pop_var)))
pathvar_RECQL3_not_pop_var_type_15 <- tidy(LR_other_RECQL3_not_pop_var, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL3_not_pop_var_type_15, "pathvar_RECQL3_not_pop_var_type_15.csv", row.names=FALSE)
system('dx upload pathvar_RECQL3_not_pop_var_type_15.csv --path Emily-folder/results_2/pathvar_RECQL3_not_pop_var_type_15.csv')

















# NRAS Type ------------------------------------------------------------------ 



# too few to do this analysis



# RIT1 Type ----------------------------------------------------------------- 

# too few to do this analysis





# ERCC3 Type ----------------------------------------------------------------- 



LR_oral_ERCC3 <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC3, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_ERCC3)
# exp(cbind(OR = coef(LR_oral_ERCC3), confint(LR_oral_ERCC3)))
pathvar_ERCC3_type_1 <- tidy(LR_oral_ERCC3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_ERCC3_type_1, "pathvar_ERCC3_type_1.csv", row.names=FALSE)
system('dx upload pathvar_ERCC3_type_1.csv --path Emily-folder/results_2/pathvar_ERCC3_type_1.csv')



LR_digestive_ERCC3 <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC3, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_ERCC3)
# exp(cbind(OR = coef(LR_digestive_ERCC3), confint(LR_digestive_ERCC3)))
pathvar_ERCC3_type_2 <- tidy(LR_digestive_ERCC3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_ERCC3_type_2, "pathvar_ERCC3_type_2.csv", row.names=FALSE)
system('dx upload pathvar_ERCC3_type_2.csv --path Emily-folder/results_2/pathvar_ERCC3_type_2.csv')



LR_respiratory_ERCC3 <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC3, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_ERCC3)
# exp(cbind(OR = coef(LR_respiratory_ERCC3), confint(LR_respiratory_ERCC3)))
pathvar_ERCC3_type_3 <- tidy(LR_respiratory_ERCC3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_ERCC3_type_3, "pathvar_ERCC3_type_3.csv", row.names=FALSE)
system('dx upload pathvar_ERCC3_type_3.csv --path Emily-folder/results_2/pathvar_ERCC3_type_3.csv')



LR_bone_ERCC3 <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC3, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_ERCC3)
# exp(cbind(OR = coef(LR_bone_ERCC3), confint(LR_bone_ERCC3)))
pathvar_ERCC3_type_4 <- tidy(LR_bone_ERCC3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_ERCC3_type_4, "pathvar_ERCC3_type_4.csv", row.names=FALSE)
system('dx upload pathvar_ERCC3_type_4.csv --path Emily-folder/results_2/pathvar_ERCC3_type_4.csv')



LR_skin_ERCC3 <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC3, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_ERCC3)
# exp(cbind(OR = coef(LR_skin_ERCC3), confint(LR_skin_ERCC3)))
pathvar_ERCC3_type_5 <- tidy(LR_skin_ERCC3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_ERCC3_type_5, "pathvar_ERCC3_type_5.csv", row.names=FALSE)
system('dx upload pathvar_ERCC3_type_5.csv --path Emily-folder/results_2/pathvar_ERCC3_type_5.csv')



LR_mesothelium_ERCC3 <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC3, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_ERCC3)
# exp(cbind(OR = coef(LR_mesothelium_ERCC3), confint(LR_mesothelium_ERCC3)))
pathvar_ERCC3_type_6 <- tidy(LR_mesothelium_ERCC3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_ERCC3_type_6, "pathvar_ERCC3_type_6.csv", row.names=FALSE)
system('dx upload pathvar_ERCC3_type_6.csv --path Emily-folder/results_2/pathvar_ERCC3_type_6.csv')



LR_breast_ERCC3 <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC3, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_ERCC3)
# exp(cbind(OR = coef(LR_breast_ERCC3), confint(LR_breast_ERCC3)))
pathvar_ERCC3_type_7 <- tidy(LR_breast_ERCC3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_ERCC3_type_7, "pathvar_ERCC3_type_7.csv", row.names=FALSE)
system('dx upload pathvar_ERCC3_type_7.csv --path Emily-folder/results_2/pathvar_ERCC3_type_7.csv')



LR_femalegenital_ERCC3 <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC3, 
                              data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_ERCC3)
# exp(cbind(OR = coef(LR_femalegenital_ERCC3), confint(LR_femalegenital_ERCC3)))
pathvar_ERCC3_type_8 <- tidy(LR_femalegenital_ERCC3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_ERCC3_type_8, "pathvar_ERCC3_type_8.csv", row.names=FALSE)
system('dx upload pathvar_ERCC3_type_8.csv --path Emily-folder/results_2/pathvar_ERCC3_type_8.csv')



LR_malegenital_ERCC3 <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC3, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_ERCC3)
# exp(cbind(OR = coef(LR_malegenital_ERCC3), confint(LR_malegenital_ERCC3)))
pathvar_ERCC3_type_9 <- tidy(LR_malegenital_ERCC3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_ERCC3_type_9, "pathvar_ERCC3_type_9.csv", row.names=FALSE)
system('dx upload pathvar_ERCC3_type_9.csv --path Emily-folder/results_2/pathvar_ERCC3_type_9.csv')



LR_urinary_ERCC3 <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC3, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_ERCC3)
# exp(cbind(OR = coef(LR_urinary_ERCC3), confint(LR_urinary_ERCC3)))
pathvar_ERCC3_type_10 <- tidy(LR_urinary_ERCC3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_ERCC3_type_10, "pathvar_ERCC3_type_10.csv", row.names=FALSE)
system('dx upload pathvar_ERCC3_type_10.csv --path Emily-folder/results_2/pathvar_ERCC3_type_10.csv')



LR_cns_ERCC3 <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC3, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_ERCC3)
# exp(cbind(OR = coef(LR_cns_ERCC3), confint(LR_cns_ERCC3)))
pathvar_ERCC3_type_11 <- tidy(LR_cns_ERCC3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_ERCC3_type_11, "pathvar_ERCC3_type_11.csv", row.names=FALSE)
system('dx upload pathvar_ERCC3_type_11.csv --path Emily-folder/results_2/pathvar_ERCC3_type_11.csv')



LR_endocrine_ERCC3 <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC3, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_ERCC3)
# exp(cbind(OR = coef(LR_endocrine_ERCC3), confint(LR_endocrine_ERCC3)))
pathvar_ERCC3_type_12 <- tidy(LR_endocrine_ERCC3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_ERCC3_type_12, "pathvar_ERCC3_type_12.csv", row.names=FALSE)
system('dx upload pathvar_ERCC3_type_12.csv --path Emily-folder/results_2/pathvar_ERCC3_type_12.csv')



LR_lymphatic_ERCC3 <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC3, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_ERCC3)
# exp(cbind(OR = coef(LR_lymphatic_ERCC3), confint(LR_lymphatic_ERCC3)))
pathvar_ERCC3_type_13 <- tidy(LR_lymphatic_ERCC3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_ERCC3_type_13, "pathvar_ERCC3_type_13.csv", row.names=FALSE)
system('dx upload pathvar_ERCC3_type_13.csv --path Emily-folder/results_2/pathvar_ERCC3_type_13.csv')



LR_secondary_ERCC3 <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC3, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_ERCC3)
# exp(cbind(OR = coef(LR_secondary_ERCC3), confint(LR_secondary_ERCC3)))
pathvar_ERCC3_type_14 <- tidy(LR_secondary_ERCC3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_ERCC3_type_14, "pathvar_ERCC3_type_14.csv", row.names=FALSE)
system('dx upload pathvar_ERCC3_type_14.csv --path Emily-folder/results_2/pathvar_ERCC3_type_14.csv')



LR_other_ERCC3 <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ERCC3, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_ERCC3)
# exp(cbind(OR = coef(LR_other_ERCC3), confint(LR_other_ERCC3)))
pathvar_ERCC3_type_15 <- tidy(LR_other_ERCC3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_ERCC3_type_15, "pathvar_ERCC3_type_15.csv", row.names=FALSE)
system('dx upload pathvar_ERCC3_type_15.csv --path Emily-folder/results_2/pathvar_ERCC3_type_15.csv')









# SMARCAL1 Type ----------------------------------------------------------------- 





LR_oral_SMARCAL1 <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SMARCAL1, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_SMARCAL1)
# exp(cbind(OR = coef(LR_oral_SMARCAL1), confint(LR_oral_SMARCAL1)))
pathvar_SMARCAL1_type_1 <- tidy(LR_oral_SMARCAL1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SMARCAL1_type_1, "pathvar_SMARCAL1_type_1.csv", row.names=FALSE)
system('dx upload pathvar_SMARCAL1_type_1.csv --path Emily-folder/results_2/pathvar_SMARCAL1_type_1.csv')



LR_digestive_SMARCAL1 <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SMARCAL1, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_SMARCAL1)
# exp(cbind(OR = coef(LR_digestive_SMARCAL1), confint(LR_digestive_SMARCAL1)))
pathvar_SMARCAL1_type_2 <- tidy(LR_digestive_SMARCAL1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SMARCAL1_type_2, "pathvar_SMARCAL1_type_2.csv", row.names=FALSE)
system('dx upload pathvar_SMARCAL1_type_2.csv --path Emily-folder/results_2/pathvar_SMARCAL1_type_2.csv')



LR_respiratory_SMARCAL1 <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SMARCAL1, 
                               data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_SMARCAL1)
# exp(cbind(OR = coef(LR_respiratory_SMARCAL1), confint(LR_respiratory_SMARCAL1)))
pathvar_SMARCAL1_type_3 <- tidy(LR_respiratory_SMARCAL1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SMARCAL1_type_3, "pathvar_SMARCAL1_type_3.csv", row.names=FALSE)
system('dx upload pathvar_SMARCAL1_type_3.csv --path Emily-folder/results_2/pathvar_SMARCAL1_type_3.csv')



LR_bone_SMARCAL1 <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SMARCAL1, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_SMARCAL1)
# exp(cbind(OR = coef(LR_bone_SMARCAL1), confint(LR_bone_SMARCAL1)))
pathvar_SMARCAL1_type_4 <- tidy(LR_bone_SMARCAL1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SMARCAL1_type_4, "pathvar_SMARCAL1_type_4.csv", row.names=FALSE)
system('dx upload pathvar_SMARCAL1_type_4.csv --path Emily-folder/results_2/pathvar_SMARCAL1_type_4.csv')



LR_skin_SMARCAL1 <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SMARCAL1, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_SMARCAL1)
# exp(cbind(OR = coef(LR_skin_SMARCAL1), confint(LR_skin_SMARCAL1)))
pathvar_SMARCAL1_type_5 <- tidy(LR_skin_SMARCAL1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SMARCAL1_type_5, "pathvar_SMARCAL1_type_5.csv", row.names=FALSE)
system('dx upload pathvar_SMARCAL1_type_5.csv --path Emily-folder/results_2/pathvar_SMARCAL1_type_5.csv')



LR_mesothelium_SMARCAL1 <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SMARCAL1, 
                               data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_SMARCAL1)
# exp(cbind(OR = coef(LR_mesothelium_SMARCAL1), confint(LR_mesothelium_SMARCAL1)))
pathvar_SMARCAL1_type_6 <- tidy(LR_mesothelium_SMARCAL1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SMARCAL1_type_6, "pathvar_SMARCAL1_type_6.csv", row.names=FALSE)
system('dx upload pathvar_SMARCAL1_type_6.csv --path Emily-folder/results_2/pathvar_SMARCAL1_type_6.csv')



LR_breast_SMARCAL1 <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SMARCAL1, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_SMARCAL1)
# exp(cbind(OR = coef(LR_breast_SMARCAL1), confint(LR_breast_SMARCAL1)))
pathvar_SMARCAL1_type_7 <- tidy(LR_breast_SMARCAL1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SMARCAL1_type_7, "pathvar_SMARCAL1_type_7.csv", row.names=FALSE)
system('dx upload pathvar_SMARCAL1_type_7.csv --path Emily-folder/results_2/pathvar_SMARCAL1_type_7.csv')



LR_femalegenital_SMARCAL1 <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SMARCAL1, 
                                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_SMARCAL1)
# exp(cbind(OR = coef(LR_femalegenital_SMARCAL1), confint(LR_femalegenital_SMARCAL1)))
pathvar_SMARCAL1_type_8 <- tidy(LR_femalegenital_SMARCAL1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SMARCAL1_type_8, "pathvar_SMARCAL1_type_8.csv", row.names=FALSE)
system('dx upload pathvar_SMARCAL1_type_8.csv --path Emily-folder/results_2/pathvar_SMARCAL1_type_8.csv')



LR_malegenital_SMARCAL1 <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SMARCAL1, 
                               data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_SMARCAL1)
# exp(cbind(OR = coef(LR_malegenital_SMARCAL1), confint(LR_malegenital_SMARCAL1)))
pathvar_SMARCAL1_type_9 <- tidy(LR_malegenital_SMARCAL1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SMARCAL1_type_9, "pathvar_SMARCAL1_type_9.csv", row.names=FALSE)
system('dx upload pathvar_SMARCAL1_type_9.csv --path Emily-folder/results_2/pathvar_SMARCAL1_type_9.csv')



LR_urinary_SMARCAL1 <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SMARCAL1, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_SMARCAL1)
# exp(cbind(OR = coef(LR_urinary_SMARCAL1), confint(LR_urinary_SMARCAL1)))
pathvar_SMARCAL1_type_10 <- tidy(LR_urinary_SMARCAL1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SMARCAL1_type_10, "pathvar_SMARCAL1_type_10.csv", row.names=FALSE)
system('dx upload pathvar_SMARCAL1_type_10.csv --path Emily-folder/results_2/pathvar_SMARCAL1_type_10.csv')



LR_cns_SMARCAL1 <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SMARCAL1, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_SMARCAL1)
# exp(cbind(OR = coef(LR_cns_SMARCAL1), confint(LR_cns_SMARCAL1)))
pathvar_SMARCAL1_type_11 <- tidy(LR_cns_SMARCAL1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SMARCAL1_type_11, "pathvar_SMARCAL1_type_11.csv", row.names=FALSE)
system('dx upload pathvar_SMARCAL1_type_11.csv --path Emily-folder/results_2/pathvar_SMARCAL1_type_11.csv')



LR_endocrine_SMARCAL1 <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SMARCAL1, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_SMARCAL1)
# exp(cbind(OR = coef(LR_endocrine_SMARCAL1), confint(LR_endocrine_SMARCAL1)))
pathvar_SMARCAL1_type_12 <- tidy(LR_endocrine_SMARCAL1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SMARCAL1_type_12, "pathvar_SMARCAL1_type_12.csv", row.names=FALSE)
system('dx upload pathvar_SMARCAL1_type_12.csv --path Emily-folder/results_2/pathvar_SMARCAL1_type_12.csv')



LR_lymphatic_SMARCAL1 <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SMARCAL1, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_SMARCAL1)
# exp(cbind(OR = coef(LR_lymphatic_SMARCAL1), confint(LR_lymphatic_SMARCAL1)))
pathvar_SMARCAL1_type_13 <- tidy(LR_lymphatic_SMARCAL1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SMARCAL1_type_13, "pathvar_SMARCAL1_type_13.csv", row.names=FALSE)
system('dx upload pathvar_SMARCAL1_type_13.csv --path Emily-folder/results_2/pathvar_SMARCAL1_type_13.csv')



LR_secondary_SMARCAL1 <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SMARCAL1, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_SMARCAL1)
# exp(cbind(OR = coef(LR_secondary_SMARCAL1), confint(LR_secondary_SMARCAL1)))
pathvar_SMARCAL1_type_14 <- tidy(LR_secondary_SMARCAL1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SMARCAL1_type_14, "pathvar_SMARCAL1_type_14.csv", row.names=FALSE)
system('dx upload pathvar_SMARCAL1_type_14.csv --path Emily-folder/results_2/pathvar_SMARCAL1_type_14.csv')



LR_other_SMARCAL1 <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SMARCAL1, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_SMARCAL1)
# exp(cbind(OR = coef(LR_other_SMARCAL1), confint(LR_other_SMARCAL1)))
pathvar_SMARCAL1_type_15 <- tidy(LR_other_SMARCAL1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SMARCAL1_type_15, "pathvar_SMARCAL1_type_15.csv", row.names=FALSE)
system('dx upload pathvar_SMARCAL1_type_15.csv --path Emily-folder/results_2/pathvar_SMARCAL1_type_15.csv')








# SBDS Type ----------------------------------------------------------------- 




LR_oral_SBDS <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SBDS, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_SBDS)
# exp(cbind(OR = coef(LR_oral_SBDS), confint(LR_oral_SBDS)))
pathvar_SBDS_type_1 <- tidy(LR_oral_SBDS, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SBDS_type_1, "pathvar_SBDS_type_1.csv", row.names=FALSE)
system('dx upload pathvar_SBDS_type_1.csv --path Emily-folder/results_2/pathvar_SBDS_type_1.csv')



LR_digestive_SBDS <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SBDS, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_SBDS)
# exp(cbind(OR = coef(LR_digestive_SBDS), confint(LR_digestive_SBDS)))
pathvar_SBDS_type_2 <- tidy(LR_digestive_SBDS, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SBDS_type_2, "pathvar_SBDS_type_2.csv", row.names=FALSE)
system('dx upload pathvar_SBDS_type_2.csv --path Emily-folder/results_2/pathvar_SBDS_type_2.csv')



LR_respiratory_SBDS <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SBDS, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_SBDS)
# exp(cbind(OR = coef(LR_respiratory_SBDS), confint(LR_respiratory_SBDS)))
pathvar_SBDS_type_3 <- tidy(LR_respiratory_SBDS, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SBDS_type_3, "pathvar_SBDS_type_3.csv", row.names=FALSE)
system('dx upload pathvar_SBDS_type_3.csv --path Emily-folder/results_2/pathvar_SBDS_type_3.csv')



LR_bone_SBDS <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SBDS, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_SBDS)
# exp(cbind(OR = coef(LR_bone_SBDS), confint(LR_bone_SBDS)))
pathvar_SBDS_type_4 <- tidy(LR_bone_SBDS, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SBDS_type_4, "pathvar_SBDS_type_4.csv", row.names=FALSE)
system('dx upload pathvar_SBDS_type_4.csv --path Emily-folder/results_2/pathvar_SBDS_type_4.csv')



LR_skin_SBDS <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SBDS, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_SBDS)
# exp(cbind(OR = coef(LR_skin_SBDS), confint(LR_skin_SBDS)))
pathvar_SBDS_type_5 <- tidy(LR_skin_SBDS, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SBDS_type_5, "pathvar_SBDS_type_5.csv", row.names=FALSE)
system('dx upload pathvar_SBDS_type_5.csv --path Emily-folder/results_2/pathvar_SBDS_type_5.csv')



LR_mesothelium_SBDS <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SBDS, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_SBDS)
# exp(cbind(OR = coef(LR_mesothelium_SBDS), confint(LR_mesothelium_SBDS)))
pathvar_SBDS_type_6 <- tidy(LR_mesothelium_SBDS, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SBDS_type_6, "pathvar_SBDS_type_6.csv", row.names=FALSE)
system('dx upload pathvar_SBDS_type_6.csv --path Emily-folder/results_2/pathvar_SBDS_type_6.csv')



LR_breast_SBDS <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SBDS, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_SBDS)
# exp(cbind(OR = coef(LR_breast_SBDS), confint(LR_breast_SBDS)))
pathvar_SBDS_type_7 <- tidy(LR_breast_SBDS, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SBDS_type_7, "pathvar_SBDS_type_7.csv", row.names=FALSE)
system('dx upload pathvar_SBDS_type_7.csv --path Emily-folder/results_2/pathvar_SBDS_type_7.csv')



LR_femalegenital_SBDS <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SBDS, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_SBDS)
# exp(cbind(OR = coef(LR_femalegenital_SBDS), confint(LR_femalegenital_SBDS)))
pathvar_SBDS_type_8 <- tidy(LR_femalegenital_SBDS, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SBDS_type_8, "pathvar_SBDS_type_8.csv", row.names=FALSE)
system('dx upload pathvar_SBDS_type_8.csv --path Emily-folder/results_2/pathvar_SBDS_type_8.csv')



LR_malegenital_SBDS <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SBDS, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_SBDS)
# exp(cbind(OR = coef(LR_malegenital_SBDS), confint(LR_malegenital_SBDS)))
pathvar_SBDS_type_9 <- tidy(LR_malegenital_SBDS, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SBDS_type_9, "pathvar_SBDS_type_9.csv", row.names=FALSE)
system('dx upload pathvar_SBDS_type_9.csv --path Emily-folder/results_2/pathvar_SBDS_type_9.csv')



LR_urinary_SBDS <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SBDS, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_SBDS)
# exp(cbind(OR = coef(LR_urinary_SBDS), confint(LR_urinary_SBDS)))
pathvar_SBDS_type_10 <- tidy(LR_urinary_SBDS, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SBDS_type_10, "pathvar_SBDS_type_10.csv", row.names=FALSE)
system('dx upload pathvar_SBDS_type_10.csv --path Emily-folder/results_2/pathvar_SBDS_type_10.csv')



LR_cns_SBDS <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SBDS, 
                   data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_SBDS)
# exp(cbind(OR = coef(LR_cns_SBDS), confint(LR_cns_SBDS)))
pathvar_SBDS_type_11 <- tidy(LR_cns_SBDS, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SBDS_type_11, "pathvar_SBDS_type_11.csv", row.names=FALSE)
system('dx upload pathvar_SBDS_type_11.csv --path Emily-folder/results_2/pathvar_SBDS_type_11.csv')



LR_endocrine_SBDS <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SBDS, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_SBDS)
# exp(cbind(OR = coef(LR_endocrine_SBDS), confint(LR_endocrine_SBDS)))
pathvar_SBDS_type_12 <- tidy(LR_endocrine_SBDS, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SBDS_type_12, "pathvar_SBDS_type_12.csv", row.names=FALSE)
system('dx upload pathvar_SBDS_type_12.csv --path Emily-folder/results_2/pathvar_SBDS_type_12.csv')



LR_lymphatic_SBDS <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SBDS, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_SBDS)
# exp(cbind(OR = coef(LR_lymphatic_SBDS), confint(LR_lymphatic_SBDS)))
pathvar_SBDS_type_13 <- tidy(LR_lymphatic_SBDS, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SBDS_type_13, "pathvar_SBDS_type_13.csv", row.names=FALSE)
system('dx upload pathvar_SBDS_type_13.csv --path Emily-folder/results_2/pathvar_SBDS_type_13.csv')



LR_secondary_SBDS <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SBDS, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_SBDS)
# exp(cbind(OR = coef(LR_secondary_SBDS), confint(LR_secondary_SBDS)))
pathvar_SBDS_type_14 <- tidy(LR_secondary_SBDS, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SBDS_type_14, "pathvar_SBDS_type_14.csv", row.names=FALSE)
system('dx upload pathvar_SBDS_type_14.csv --path Emily-folder/results_2/pathvar_SBDS_type_14.csv')



LR_other_SBDS <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SBDS, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_SBDS)
# exp(cbind(OR = coef(LR_other_SBDS), confint(LR_other_SBDS)))
pathvar_SBDS_type_15 <- tidy(LR_other_SBDS, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_SBDS_type_15, "pathvar_SBDS_type_15.csv", row.names=FALSE)
system('dx upload pathvar_SBDS_type_15.csv --path Emily-folder/results_2/pathvar_SBDS_type_15.csv')








# RECQL4 Type ----------------------------------------------------------------- 





LR_oral_RECQL4 <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL4, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_RECQL4)
# exp(cbind(OR = coef(LR_oral_RECQL4), confint(LR_oral_RECQL4)))
pathvar_RECQL4_type_1 <- tidy(LR_oral_RECQL4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL4_type_1, "pathvar_RECQL4_type_1.csv", row.names=FALSE)
system('dx upload pathvar_RECQL4_type_1.csv --path Emily-folder/results_2/pathvar_RECQL4_type_1.csv')



LR_digestive_RECQL4 <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL4, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_RECQL4)
# exp(cbind(OR = coef(LR_digestive_RECQL4), confint(LR_digestive_RECQL4)))
pathvar_RECQL4_type_2 <- tidy(LR_digestive_RECQL4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL4_type_2, "pathvar_RECQL4_type_2.csv", row.names=FALSE)
system('dx upload pathvar_RECQL4_type_2.csv --path Emily-folder/results_2/pathvar_RECQL4_type_2.csv')



LR_respiratory_RECQL4 <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL4, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_RECQL4)
# exp(cbind(OR = coef(LR_respiratory_RECQL4), confint(LR_respiratory_RECQL4)))
pathvar_RECQL4_type_3 <- tidy(LR_respiratory_RECQL4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL4_type_3, "pathvar_RECQL4_type_3.csv", row.names=FALSE)
system('dx upload pathvar_RECQL4_type_3.csv --path Emily-folder/results_2/pathvar_RECQL4_type_3.csv')



LR_bone_RECQL4 <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL4, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_RECQL4)
# exp(cbind(OR = coef(LR_bone_RECQL4), confint(LR_bone_RECQL4)))
pathvar_RECQL4_type_4 <- tidy(LR_bone_RECQL4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL4_type_4, "pathvar_RECQL4_type_4.csv", row.names=FALSE)
system('dx upload pathvar_RECQL4_type_4.csv --path Emily-folder/results_2/pathvar_RECQL4_type_4.csv')



LR_skin_RECQL4 <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL4, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_RECQL4)
# exp(cbind(OR = coef(LR_skin_RECQL4), confint(LR_skin_RECQL4)))
pathvar_RECQL4_type_5 <- tidy(LR_skin_RECQL4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL4_type_5, "pathvar_RECQL4_type_5.csv", row.names=FALSE)
system('dx upload pathvar_RECQL4_type_5.csv --path Emily-folder/results_2/pathvar_RECQL4_type_5.csv')



LR_mesothelium_RECQL4 <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL4, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_RECQL4)
# exp(cbind(OR = coef(LR_mesothelium_RECQL4), confint(LR_mesothelium_RECQL4)))
pathvar_RECQL4_type_6 <- tidy(LR_mesothelium_RECQL4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL4_type_6, "pathvar_RECQL4_type_6.csv", row.names=FALSE)
system('dx upload pathvar_RECQL4_type_6.csv --path Emily-folder/results_2/pathvar_RECQL4_type_6.csv')



LR_breast_RECQL4 <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL4, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_RECQL4)
# exp(cbind(OR = coef(LR_breast_RECQL4), confint(LR_breast_RECQL4)))
pathvar_RECQL4_type_7 <- tidy(LR_breast_RECQL4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL4_type_7, "pathvar_RECQL4_type_7.csv", row.names=FALSE)
system('dx upload pathvar_RECQL4_type_7.csv --path Emily-folder/results_2/pathvar_RECQL4_type_7.csv')



LR_femalegenital_RECQL4 <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL4, 
                               data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_RECQL4)
# exp(cbind(OR = coef(LR_femalegenital_RECQL4), confint(LR_femalegenital_RECQL4)))
pathvar_RECQL4_type_8 <- tidy(LR_femalegenital_RECQL4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL4_type_8, "pathvar_RECQL4_type_8.csv", row.names=FALSE)
system('dx upload pathvar_RECQL4_type_8.csv --path Emily-folder/results_2/pathvar_RECQL4_type_8.csv')



LR_malegenital_RECQL4 <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL4, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_RECQL4)
# exp(cbind(OR = coef(LR_malegenital_RECQL4), confint(LR_malegenital_RECQL4)))
pathvar_RECQL4_type_9 <- tidy(LR_malegenital_RECQL4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL4_type_9, "pathvar_RECQL4_type_9.csv", row.names=FALSE)
system('dx upload pathvar_RECQL4_type_9.csv --path Emily-folder/results_2/pathvar_RECQL4_type_9.csv')



LR_urinary_RECQL4 <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL4, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_RECQL4)
# exp(cbind(OR = coef(LR_urinary_RECQL4), confint(LR_urinary_RECQL4)))
pathvar_RECQL4_type_10 <- tidy(LR_urinary_RECQL4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL4_type_10, "pathvar_RECQL4_type_10.csv", row.names=FALSE)
system('dx upload pathvar_RECQL4_type_10.csv --path Emily-folder/results_2/pathvar_RECQL4_type_10.csv')



LR_cns_RECQL4 <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL4, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_RECQL4)
# exp(cbind(OR = coef(LR_cns_RECQL4), confint(LR_cns_RECQL4)))
pathvar_RECQL4_type_11 <- tidy(LR_cns_RECQL4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL4_type_11, "pathvar_RECQL4_type_11.csv", row.names=FALSE)
system('dx upload pathvar_RECQL4_type_11.csv --path Emily-folder/results_2/pathvar_RECQL4_type_11.csv')



LR_endocrine_RECQL4 <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL4, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_RECQL4)
# exp(cbind(OR = coef(LR_endocrine_RECQL4), confint(LR_endocrine_RECQL4)))
pathvar_RECQL4_type_12 <- tidy(LR_endocrine_RECQL4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL4_type_12, "pathvar_RECQL4_type_12.csv", row.names=FALSE)
system('dx upload pathvar_RECQL4_type_12.csv --path Emily-folder/results_2/pathvar_RECQL4_type_12.csv')



LR_lymphatic_RECQL4 <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL4, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_RECQL4)
# exp(cbind(OR = coef(LR_lymphatic_RECQL4), confint(LR_lymphatic_RECQL4)))
pathvar_RECQL4_type_13 <- tidy(LR_lymphatic_RECQL4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL4_type_13, "pathvar_RECQL4_type_13.csv", row.names=FALSE)
system('dx upload pathvar_RECQL4_type_13.csv --path Emily-folder/results_2/pathvar_RECQL4_type_13.csv')



LR_secondary_RECQL4 <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL4, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_RECQL4)
# exp(cbind(OR = coef(LR_secondary_RECQL4), confint(LR_secondary_RECQL4)))
pathvar_RECQL4_type_14 <- tidy(LR_secondary_RECQL4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL4_type_14, "pathvar_RECQL4_type_14.csv", row.names=FALSE)
system('dx upload pathvar_RECQL4_type_14.csv --path Emily-folder/results_2/pathvar_RECQL4_type_14.csv')



LR_other_RECQL4 <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RECQL4, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_RECQL4)
# exp(cbind(OR = coef(LR_other_RECQL4), confint(LR_other_RECQL4)))
pathvar_RECQL4_type_15 <- tidy(LR_other_RECQL4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RECQL4_type_15, "pathvar_RECQL4_type_15.csv", row.names=FALSE)
system('dx upload pathvar_RECQL4_type_15.csv --path Emily-folder/results_2/pathvar_RECQL4_type_15.csv')








# CBL Type ----------------------------------------------------------------- 

# too few participants!






# PTPN11 Type ----------------------------------------------------------------- 



# too few participants





# KRAS Type ----------------------------------------------------------------- 

# too few participants







# POLE Type ----------------------------------------------------------------- 



LR_oral_POLE <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POLE, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_POLE)
# exp(cbind(OR = coef(LR_oral_POLE), confint(LR_oral_POLE)))
pathvar_POLE_type_1 <- tidy(LR_oral_POLE, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_POLE_type_1, "pathvar_POLE_type_1.csv", row.names=FALSE)
system('dx upload pathvar_POLE_type_1.csv --path Emily-folder/results_2/pathvar_POLE_type_1.csv')



LR_digestive_POLE <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POLE, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_POLE)
# exp(cbind(OR = coef(LR_digestive_POLE), confint(LR_digestive_POLE)))
pathvar_POLE_type_2 <- tidy(LR_digestive_POLE, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_POLE_type_2, "pathvar_POLE_type_2.csv", row.names=FALSE)
system('dx upload pathvar_POLE_type_2.csv --path Emily-folder/results_2/pathvar_POLE_type_2.csv')



LR_respiratory_POLE <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POLE, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_POLE)
# exp(cbind(OR = coef(LR_respiratory_POLE), confint(LR_respiratory_POLE)))
pathvar_POLE_type_3 <- tidy(LR_respiratory_POLE, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_POLE_type_3, "pathvar_POLE_type_3.csv", row.names=FALSE)
system('dx upload pathvar_POLE_type_3.csv --path Emily-folder/results_2/pathvar_POLE_type_3.csv')



LR_bone_POLE <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POLE, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_POLE)
# exp(cbind(OR = coef(LR_bone_POLE), confint(LR_bone_POLE)))
pathvar_POLE_type_4 <- tidy(LR_bone_POLE, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_POLE_type_4, "pathvar_POLE_type_4.csv", row.names=FALSE)
system('dx upload pathvar_POLE_type_4.csv --path Emily-folder/results_2/pathvar_POLE_type_4.csv')



LR_skin_POLE <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POLE, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_POLE)
# exp(cbind(OR = coef(LR_skin_POLE), confint(LR_skin_POLE)))
pathvar_POLE_type_5 <- tidy(LR_skin_POLE, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_POLE_type_5, "pathvar_POLE_type_5.csv", row.names=FALSE)
system('dx upload pathvar_POLE_type_5.csv --path Emily-folder/results_2/pathvar_POLE_type_5.csv')



LR_mesothelium_POLE <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POLE, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_POLE)
# exp(cbind(OR = coef(LR_mesothelium_POLE), confint(LR_mesothelium_POLE)))
pathvar_POLE_type_6 <- tidy(LR_mesothelium_POLE, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_POLE_type_6, "pathvar_POLE_type_6.csv", row.names=FALSE)
system('dx upload pathvar_POLE_type_6.csv --path Emily-folder/results_2/pathvar_POLE_type_6.csv')



LR_breast_POLE <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POLE, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_POLE)
# exp(cbind(OR = coef(LR_breast_POLE), confint(LR_breast_POLE)))
pathvar_POLE_type_7 <- tidy(LR_breast_POLE, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_POLE_type_7, "pathvar_POLE_type_7.csv", row.names=FALSE)
system('dx upload pathvar_POLE_type_7.csv --path Emily-folder/results_2/pathvar_POLE_type_7.csv')



LR_femalegenital_POLE <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POLE, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_POLE)
# exp(cbind(OR = coef(LR_femalegenital_POLE), confint(LR_femalegenital_POLE)))
pathvar_POLE_type_8 <- tidy(LR_femalegenital_POLE, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_POLE_type_8, "pathvar_POLE_type_8.csv", row.names=FALSE)
system('dx upload pathvar_POLE_type_8.csv --path Emily-folder/results_2/pathvar_POLE_type_8.csv')



LR_malegenital_POLE <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POLE, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_POLE)
# exp(cbind(OR = coef(LR_malegenital_POLE), confint(LR_malegenital_POLE)))
pathvar_POLE_type_9 <- tidy(LR_malegenital_POLE, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_POLE_type_9, "pathvar_POLE_type_9.csv", row.names=FALSE)
system('dx upload pathvar_POLE_type_9.csv --path Emily-folder/results_2/pathvar_POLE_type_9.csv')



LR_urinary_POLE <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POLE, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_POLE)
# exp(cbind(OR = coef(LR_urinary_POLE), confint(LR_urinary_POLE)))
pathvar_POLE_type_10 <- tidy(LR_urinary_POLE, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_POLE_type_10, "pathvar_POLE_type_10.csv", row.names=FALSE)
system('dx upload pathvar_POLE_type_10.csv --path Emily-folder/results_2/pathvar_POLE_type_10.csv')



LR_cns_POLE <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POLE, 
                   data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_POLE)
# exp(cbind(OR = coef(LR_cns_POLE), confint(LR_cns_POLE)))
pathvar_POLE_type_11 <- tidy(LR_cns_POLE, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_POLE_type_11, "pathvar_POLE_type_11.csv", row.names=FALSE)
system('dx upload pathvar_POLE_type_11.csv --path Emily-folder/results_2/pathvar_POLE_type_11.csv')



LR_endocrine_POLE <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POLE, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_POLE)
# exp(cbind(OR = coef(LR_endocrine_POLE), confint(LR_endocrine_POLE)))
pathvar_POLE_type_12 <- tidy(LR_endocrine_POLE, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_POLE_type_12, "pathvar_POLE_type_12.csv", row.names=FALSE)
system('dx upload pathvar_POLE_type_12.csv --path Emily-folder/results_2/pathvar_POLE_type_12.csv')



LR_lymphatic_POLE <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POLE, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_POLE)
# exp(cbind(OR = coef(LR_lymphatic_POLE), confint(LR_lymphatic_POLE)))
pathvar_POLE_type_13 <- tidy(LR_lymphatic_POLE, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_POLE_type_13, "pathvar_POLE_type_13.csv", row.names=FALSE)
system('dx upload pathvar_POLE_type_13.csv --path Emily-folder/results_2/pathvar_POLE_type_13.csv')



LR_secondary_POLE <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POLE, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_POLE)
# exp(cbind(OR = coef(LR_secondary_POLE), confint(LR_secondary_POLE)))
pathvar_POLE_type_14 <- tidy(LR_secondary_POLE, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_POLE_type_14, "pathvar_POLE_type_14.csv", row.names=FALSE)
system('dx upload pathvar_POLE_type_14.csv --path Emily-folder/results_2/pathvar_POLE_type_14.csv')



LR_other_POLE <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+POLE, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_POLE)
# exp(cbind(OR = coef(LR_other_POLE), confint(LR_other_POLE)))
pathvar_POLE_type_15 <- tidy(LR_other_POLE, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_POLE_type_15, "pathvar_POLE_type_15.csv", row.names=FALSE)
system('dx upload pathvar_POLE_type_15.csv --path Emily-folder/results_2/pathvar_POLE_type_15.csv')










# LIG4 Type ----------------------------------------------------------------- 





LR_oral_LIG4 <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LIG4, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_LIG4)
# exp(cbind(OR = coef(LR_oral_LIG4), confint(LR_oral_LIG4)))
pathvar_LIG4_type_1 <- tidy(LR_oral_LIG4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LIG4_type_1, "pathvar_LIG4_type_1.csv", row.names=FALSE)
system('dx upload pathvar_LIG4_type_1.csv --path Emily-folder/results_2/pathvar_LIG4_type_1.csv')



LR_digestive_LIG4 <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LIG4, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_LIG4)
# exp(cbind(OR = coef(LR_digestive_LIG4), confint(LR_digestive_LIG4)))
pathvar_LIG4_type_2 <- tidy(LR_digestive_LIG4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LIG4_type_2, "pathvar_LIG4_type_2.csv", row.names=FALSE)
system('dx upload pathvar_LIG4_type_2.csv --path Emily-folder/results_2/pathvar_LIG4_type_2.csv')



LR_respiratory_LIG4 <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LIG4, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_LIG4)
# exp(cbind(OR = coef(LR_respiratory_LIG4), confint(LR_respiratory_LIG4)))
pathvar_LIG4_type_3 <- tidy(LR_respiratory_LIG4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LIG4_type_3, "pathvar_LIG4_type_3.csv", row.names=FALSE)
system('dx upload pathvar_LIG4_type_3.csv --path Emily-folder/results_2/pathvar_LIG4_type_3.csv')



LR_bone_LIG4 <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LIG4, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_LIG4)
# exp(cbind(OR = coef(LR_bone_LIG4), confint(LR_bone_LIG4)))
pathvar_LIG4_type_4 <- tidy(LR_bone_LIG4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LIG4_type_4, "pathvar_LIG4_type_4.csv", row.names=FALSE)
system('dx upload pathvar_LIG4_type_4.csv --path Emily-folder/results_2/pathvar_LIG4_type_4.csv')



LR_skin_LIG4 <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LIG4, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_LIG4)
# exp(cbind(OR = coef(LR_skin_LIG4), confint(LR_skin_LIG4)))
pathvar_LIG4_type_5 <- tidy(LR_skin_LIG4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LIG4_type_5, "pathvar_LIG4_type_5.csv", row.names=FALSE)
system('dx upload pathvar_LIG4_type_5.csv --path Emily-folder/results_2/pathvar_LIG4_type_5.csv')



LR_mesothelium_LIG4 <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LIG4, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_LIG4)
# exp(cbind(OR = coef(LR_mesothelium_LIG4), confint(LR_mesothelium_LIG4)))
pathvar_LIG4_type_6 <- tidy(LR_mesothelium_LIG4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LIG4_type_6, "pathvar_LIG4_type_6.csv", row.names=FALSE)
system('dx upload pathvar_LIG4_type_6.csv --path Emily-folder/results_2/pathvar_LIG4_type_6.csv')



LR_breast_LIG4 <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LIG4, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_LIG4)
# exp(cbind(OR = coef(LR_breast_LIG4), confint(LR_breast_LIG4)))
pathvar_LIG4_type_7 <- tidy(LR_breast_LIG4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LIG4_type_7, "pathvar_LIG4_type_7.csv", row.names=FALSE)
system('dx upload pathvar_LIG4_type_7.csv --path Emily-folder/results_2/pathvar_LIG4_type_7.csv')



LR_femalegenital_LIG4 <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LIG4, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_LIG4)
# exp(cbind(OR = coef(LR_femalegenital_LIG4), confint(LR_femalegenital_LIG4)))
pathvar_LIG4_type_8 <- tidy(LR_femalegenital_LIG4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LIG4_type_8, "pathvar_LIG4_type_8.csv", row.names=FALSE)
system('dx upload pathvar_LIG4_type_8.csv --path Emily-folder/results_2/pathvar_LIG4_type_8.csv')



LR_malegenital_LIG4 <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LIG4, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_LIG4)
# exp(cbind(OR = coef(LR_malegenital_LIG4), confint(LR_malegenital_LIG4)))
pathvar_LIG4_type_9 <- tidy(LR_malegenital_LIG4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LIG4_type_9, "pathvar_LIG4_type_9.csv", row.names=FALSE)
system('dx upload pathvar_LIG4_type_9.csv --path Emily-folder/results_2/pathvar_LIG4_type_9.csv')



LR_urinary_LIG4 <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LIG4, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_LIG4)
# exp(cbind(OR = coef(LR_urinary_LIG4), confint(LR_urinary_LIG4)))
pathvar_LIG4_type_10 <- tidy(LR_urinary_LIG4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LIG4_type_10, "pathvar_LIG4_type_10.csv", row.names=FALSE)
system('dx upload pathvar_LIG4_type_10.csv --path Emily-folder/results_2/pathvar_LIG4_type_10.csv')



LR_cns_LIG4 <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LIG4, 
                   data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_LIG4)
# exp(cbind(OR = coef(LR_cns_LIG4), confint(LR_cns_LIG4)))
pathvar_LIG4_type_11 <- tidy(LR_cns_LIG4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LIG4_type_11, "pathvar_LIG4_type_11.csv", row.names=FALSE)
system('dx upload pathvar_LIG4_type_11.csv --path Emily-folder/results_2/pathvar_LIG4_type_11.csv')



LR_endocrine_LIG4 <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LIG4, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_LIG4)
# exp(cbind(OR = coef(LR_endocrine_LIG4), confint(LR_endocrine_LIG4)))
pathvar_LIG4_type_12 <- tidy(LR_endocrine_LIG4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LIG4_type_12, "pathvar_LIG4_type_12.csv", row.names=FALSE)
system('dx upload pathvar_LIG4_type_12.csv --path Emily-folder/results_2/pathvar_LIG4_type_12.csv')



LR_lymphatic_LIG4 <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LIG4, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_LIG4)
# exp(cbind(OR = coef(LR_lymphatic_LIG4), confint(LR_lymphatic_LIG4)))
pathvar_LIG4_type_13 <- tidy(LR_lymphatic_LIG4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LIG4_type_13, "pathvar_LIG4_type_13.csv", row.names=FALSE)
system('dx upload pathvar_LIG4_type_13.csv --path Emily-folder/results_2/pathvar_LIG4_type_13.csv')



LR_secondary_LIG4 <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LIG4, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_LIG4)
# exp(cbind(OR = coef(LR_secondary_LIG4), confint(LR_secondary_LIG4)))
pathvar_LIG4_type_14 <- tidy(LR_secondary_LIG4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LIG4_type_14, "pathvar_LIG4_type_14.csv", row.names=FALSE)
system('dx upload pathvar_LIG4_type_14.csv --path Emily-folder/results_2/pathvar_LIG4_type_14.csv')



LR_other_LIG4 <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LIG4, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_LIG4)
# exp(cbind(OR = coef(LR_other_LIG4), confint(LR_other_LIG4)))
pathvar_LIG4_type_15 <- tidy(LR_other_LIG4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LIG4_type_15, "pathvar_LIG4_type_15.csv", row.names=FALSE)
system('dx upload pathvar_LIG4_type_15.csv --path Emily-folder/results_2/pathvar_LIG4_type_15.csv')









# BUB1B Type ----------------------------------------------------------------- 



LR_oral_BUB1B <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_BUB1B)
# exp(cbind(OR = coef(LR_oral_BUB1B), confint(LR_oral_BUB1B)))
pathvar_BUB1B_type_1 <- tidy(LR_oral_BUB1B, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BUB1B_type_1, "pathvar_BUB1B_type_1.csv", row.names=FALSE)
system('dx upload pathvar_BUB1B_type_1.csv --path Emily-folder/results_2/pathvar_BUB1B_type_1.csv')



LR_digestive_BUB1B <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_BUB1B)
# exp(cbind(OR = coef(LR_digestive_BUB1B), confint(LR_digestive_BUB1B)))
pathvar_BUB1B_type_2 <- tidy(LR_digestive_BUB1B, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BUB1B_type_2, "pathvar_BUB1B_type_2.csv", row.names=FALSE)
system('dx upload pathvar_BUB1B_type_2.csv --path Emily-folder/results_2/pathvar_BUB1B_type_2.csv')



LR_respiratory_BUB1B <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_BUB1B)
# exp(cbind(OR = coef(LR_respiratory_BUB1B), confint(LR_respiratory_BUB1B)))
pathvar_BUB1B_type_3 <- tidy(LR_respiratory_BUB1B, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BUB1B_type_3, "pathvar_BUB1B_type_3.csv", row.names=FALSE)
system('dx upload pathvar_BUB1B_type_3.csv --path Emily-folder/results_2/pathvar_BUB1B_type_3.csv')



LR_bone_BUB1B <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_BUB1B)
# exp(cbind(OR = coef(LR_bone_BUB1B), confint(LR_bone_BUB1B)))
pathvar_BUB1B_type_4 <- tidy(LR_bone_BUB1B, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BUB1B_type_4, "pathvar_BUB1B_type_4.csv", row.names=FALSE)
system('dx upload pathvar_BUB1B_type_4.csv --path Emily-folder/results_2/pathvar_BUB1B_type_4.csv')



LR_skin_BUB1B <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_BUB1B)
# exp(cbind(OR = coef(LR_skin_BUB1B), confint(LR_skin_BUB1B)))
pathvar_BUB1B_type_5 <- tidy(LR_skin_BUB1B, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BUB1B_type_5, "pathvar_BUB1B_type_5.csv", row.names=FALSE)
system('dx upload pathvar_BUB1B_type_5.csv --path Emily-folder/results_2/pathvar_BUB1B_type_5.csv')



LR_mesothelium_BUB1B <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_BUB1B)
# exp(cbind(OR = coef(LR_mesothelium_BUB1B), confint(LR_mesothelium_BUB1B)))
pathvar_BUB1B_type_6 <- tidy(LR_mesothelium_BUB1B, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BUB1B_type_6, "pathvar_BUB1B_type_6.csv", row.names=FALSE)
system('dx upload pathvar_BUB1B_type_6.csv --path Emily-folder/results_2/pathvar_BUB1B_type_6.csv')



LR_breast_BUB1B <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_BUB1B)
# exp(cbind(OR = coef(LR_breast_BUB1B), confint(LR_breast_BUB1B)))
pathvar_BUB1B_type_7 <- tidy(LR_breast_BUB1B, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BUB1B_type_7, "pathvar_BUB1B_type_7.csv", row.names=FALSE)
system('dx upload pathvar_BUB1B_type_7.csv --path Emily-folder/results_2/pathvar_BUB1B_type_7.csv')



LR_femalegenital_BUB1B <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B, 
                              data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_BUB1B)
# exp(cbind(OR = coef(LR_femalegenital_BUB1B), confint(LR_femalegenital_BUB1B)))
pathvar_BUB1B_type_8 <- tidy(LR_femalegenital_BUB1B, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BUB1B_type_8, "pathvar_BUB1B_type_8.csv", row.names=FALSE)
system('dx upload pathvar_BUB1B_type_8.csv --path Emily-folder/results_2/pathvar_BUB1B_type_8.csv')



LR_malegenital_BUB1B <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_BUB1B)
# exp(cbind(OR = coef(LR_malegenital_BUB1B), confint(LR_malegenital_BUB1B)))
pathvar_BUB1B_type_9 <- tidy(LR_malegenital_BUB1B, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BUB1B_type_9, "pathvar_BUB1B_type_9.csv", row.names=FALSE)
system('dx upload pathvar_BUB1B_type_9.csv --path Emily-folder/results_2/pathvar_BUB1B_type_9.csv')



LR_urinary_BUB1B <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_BUB1B)
# exp(cbind(OR = coef(LR_urinary_BUB1B), confint(LR_urinary_BUB1B)))
pathvar_BUB1B_type_10 <- tidy(LR_urinary_BUB1B, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BUB1B_type_10, "pathvar_BUB1B_type_10.csv", row.names=FALSE)
system('dx upload pathvar_BUB1B_type_10.csv --path Emily-folder/results_2/pathvar_BUB1B_type_10.csv')



LR_cns_BUB1B <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_BUB1B)
# exp(cbind(OR = coef(LR_cns_BUB1B), confint(LR_cns_BUB1B)))
pathvar_BUB1B_type_11 <- tidy(LR_cns_BUB1B, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BUB1B_type_11, "pathvar_BUB1B_type_11.csv", row.names=FALSE)
system('dx upload pathvar_BUB1B_type_11.csv --path Emily-folder/results_2/pathvar_BUB1B_type_11.csv')



LR_endocrine_BUB1B <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_BUB1B)
# exp(cbind(OR = coef(LR_endocrine_BUB1B), confint(LR_endocrine_BUB1B)))
pathvar_BUB1B_type_12 <- tidy(LR_endocrine_BUB1B, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BUB1B_type_12, "pathvar_BUB1B_type_12.csv", row.names=FALSE)
system('dx upload pathvar_BUB1B_type_12.csv --path Emily-folder/results_2/pathvar_BUB1B_type_12.csv')



LR_lymphatic_BUB1B <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_BUB1B)
# exp(cbind(OR = coef(LR_lymphatic_BUB1B), confint(LR_lymphatic_BUB1B)))
pathvar_BUB1B_type_13 <- tidy(LR_lymphatic_BUB1B, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BUB1B_type_13, "pathvar_BUB1B_type_13.csv", row.names=FALSE)
system('dx upload pathvar_BUB1B_type_13.csv --path Emily-folder/results_2/pathvar_BUB1B_type_13.csv')



LR_secondary_BUB1B <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_BUB1B)
# exp(cbind(OR = coef(LR_secondary_BUB1B), confint(LR_secondary_BUB1B)))
pathvar_BUB1B_type_14 <- tidy(LR_secondary_BUB1B, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BUB1B_type_14, "pathvar_BUB1B_type_14.csv", row.names=FALSE)
system('dx upload pathvar_BUB1B_type_14.csv --path Emily-folder/results_2/pathvar_BUB1B_type_14.csv')



LR_other_BUB1B <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BUB1B, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_BUB1B)
# exp(cbind(OR = coef(LR_other_BUB1B), confint(LR_other_BUB1B)))
pathvar_BUB1B_type_15 <- tidy(LR_other_BUB1B, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BUB1B_type_15, "pathvar_BUB1B_type_15.csv", row.names=FALSE)
system('dx upload pathvar_BUB1B_type_15.csv --path Emily-folder/results_2/pathvar_BUB1B_type_15.csv')









# FANCA Type ----------------------------------------------------------------- 



LR_oral_FANCA <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FANCA, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_FANCA)
# exp(cbind(OR = coef(LR_oral_FANCA), confint(LR_oral_FANCA)))
pathvar_FANCA_type_1 <- tidy(LR_oral_FANCA, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_FANCA_type_1, "pathvar_FANCA_type_1.csv", row.names=FALSE)
system('dx upload pathvar_FANCA_type_1.csv --path Emily-folder/results_2/pathvar_FANCA_type_1.csv')



LR_digestive_FANCA <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FANCA, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_FANCA)
# exp(cbind(OR = coef(LR_digestive_FANCA), confint(LR_digestive_FANCA)))
pathvar_FANCA_type_2 <- tidy(LR_digestive_FANCA, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_FANCA_type_2, "pathvar_FANCA_type_2.csv", row.names=FALSE)
system('dx upload pathvar_FANCA_type_2.csv --path Emily-folder/results_2/pathvar_FANCA_type_2.csv')



LR_respiratory_FANCA <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FANCA, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_FANCA)
# exp(cbind(OR = coef(LR_respiratory_FANCA), confint(LR_respiratory_FANCA)))
pathvar_FANCA_type_3 <- tidy(LR_respiratory_FANCA, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_FANCA_type_3, "pathvar_FANCA_type_3.csv", row.names=FALSE)
system('dx upload pathvar_FANCA_type_3.csv --path Emily-folder/results_2/pathvar_FANCA_type_3.csv')



LR_bone_FANCA <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FANCA, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_FANCA)
# exp(cbind(OR = coef(LR_bone_FANCA), confint(LR_bone_FANCA)))
pathvar_FANCA_type_4 <- tidy(LR_bone_FANCA, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_FANCA_type_4, "pathvar_FANCA_type_4.csv", row.names=FALSE)
system('dx upload pathvar_FANCA_type_4.csv --path Emily-folder/results_2/pathvar_FANCA_type_4.csv')



LR_skin_FANCA <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FANCA, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_FANCA)
# exp(cbind(OR = coef(LR_skin_FANCA), confint(LR_skin_FANCA)))
pathvar_FANCA_type_5 <- tidy(LR_skin_FANCA, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_FANCA_type_5, "pathvar_FANCA_type_5.csv", row.names=FALSE)
system('dx upload pathvar_FANCA_type_5.csv --path Emily-folder/results_2/pathvar_FANCA_type_5.csv')



LR_mesothelium_FANCA <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FANCA, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_FANCA)
# exp(cbind(OR = coef(LR_mesothelium_FANCA), confint(LR_mesothelium_FANCA)))
pathvar_FANCA_type_6 <- tidy(LR_mesothelium_FANCA, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_FANCA_type_6, "pathvar_FANCA_type_6.csv", row.names=FALSE)
system('dx upload pathvar_FANCA_type_6.csv --path Emily-folder/results_2/pathvar_FANCA_type_6.csv')



LR_breast_FANCA <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FANCA, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_FANCA)
# exp(cbind(OR = coef(LR_breast_FANCA), confint(LR_breast_FANCA)))
pathvar_FANCA_type_7 <- tidy(LR_breast_FANCA, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_FANCA_type_7, "pathvar_FANCA_type_7.csv", row.names=FALSE)
system('dx upload pathvar_FANCA_type_7.csv --path Emily-folder/results_2/pathvar_FANCA_type_7.csv')



LR_femalegenital_FANCA <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FANCA, 
                              data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_FANCA)
# exp(cbind(OR = coef(LR_femalegenital_FANCA), confint(LR_femalegenital_FANCA)))
pathvar_FANCA_type_8 <- tidy(LR_femalegenital_FANCA, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_FANCA_type_8, "pathvar_FANCA_type_8.csv", row.names=FALSE)
system('dx upload pathvar_FANCA_type_8.csv --path Emily-folder/results_2/pathvar_FANCA_type_8.csv')



LR_malegenital_FANCA <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FANCA, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_FANCA)
# exp(cbind(OR = coef(LR_malegenital_FANCA), confint(LR_malegenital_FANCA)))
pathvar_FANCA_type_9 <- tidy(LR_malegenital_FANCA, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_FANCA_type_9, "pathvar_FANCA_type_9.csv", row.names=FALSE)
system('dx upload pathvar_FANCA_type_9.csv --path Emily-folder/results_2/pathvar_FANCA_type_9.csv')



LR_urinary_FANCA <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FANCA, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_FANCA)
# exp(cbind(OR = coef(LR_urinary_FANCA), confint(LR_urinary_FANCA)))
pathvar_FANCA_type_10 <- tidy(LR_urinary_FANCA, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_FANCA_type_10, "pathvar_FANCA_type_10.csv", row.names=FALSE)
system('dx upload pathvar_FANCA_type_10.csv --path Emily-folder/results_2/pathvar_FANCA_type_10.csv')



LR_cns_FANCA <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FANCA, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_FANCA)
# exp(cbind(OR = coef(LR_cns_FANCA), confint(LR_cns_FANCA)))
pathvar_FANCA_type_11 <- tidy(LR_cns_FANCA, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_FANCA_type_11, "pathvar_FANCA_type_11.csv", row.names=FALSE)
system('dx upload pathvar_FANCA_type_11.csv --path Emily-folder/results_2/pathvar_FANCA_type_11.csv')



LR_endocrine_FANCA <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FANCA, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_FANCA)
# exp(cbind(OR = coef(LR_endocrine_FANCA), confint(LR_endocrine_FANCA)))
pathvar_FANCA_type_12 <- tidy(LR_endocrine_FANCA, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_FANCA_type_12, "pathvar_FANCA_type_12.csv", row.names=FALSE)
system('dx upload pathvar_FANCA_type_12.csv --path Emily-folder/results_2/pathvar_FANCA_type_12.csv')



LR_lymphatic_FANCA <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FANCA, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_FANCA)
# exp(cbind(OR = coef(LR_lymphatic_FANCA), confint(LR_lymphatic_FANCA)))
pathvar_FANCA_type_13 <- tidy(LR_lymphatic_FANCA, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_FANCA_type_13, "pathvar_FANCA_type_13.csv", row.names=FALSE)
system('dx upload pathvar_FANCA_type_13.csv --path Emily-folder/results_2/pathvar_FANCA_type_13.csv')



LR_secondary_FANCA <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FANCA, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_FANCA)
# exp(cbind(OR = coef(LR_secondary_FANCA), confint(LR_secondary_FANCA)))
pathvar_FANCA_type_14 <- tidy(LR_secondary_FANCA, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_FANCA_type_14, "pathvar_FANCA_type_14.csv", row.names=FALSE)
system('dx upload pathvar_FANCA_type_14.csv --path Emily-folder/results_2/pathvar_FANCA_type_14.csv')



LR_other_FANCA <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FANCA, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_FANCA)
# exp(cbind(OR = coef(LR_other_FANCA), confint(LR_other_FANCA)))
pathvar_FANCA_type_15 <- tidy(LR_other_FANCA, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_FANCA_type_15, "pathvar_FANCA_type_15.csv", row.names=FALSE)
system('dx upload pathvar_FANCA_type_15.csv --path Emily-folder/results_2/pathvar_FANCA_type_15.csv')








# CREBBP Type ----------------------------------------------------------------- 





LR_oral_CREBBP <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CREBBP, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_CREBBP)
# exp(cbind(OR = coef(LR_oral_CREBBP), confint(LR_oral_CREBBP)))
pathvar_CREBBP_type_1 <- tidy(LR_oral_CREBBP, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CREBBP_type_1, "pathvar_CREBBP_type_1.csv", row.names=FALSE)
system('dx upload pathvar_CREBBP_type_1.csv --path Emily-folder/results_2/pathvar_CREBBP_type_1.csv')



LR_digestive_CREBBP <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CREBBP, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_CREBBP)
# exp(cbind(OR = coef(LR_digestive_CREBBP), confint(LR_digestive_CREBBP)))
pathvar_CREBBP_type_2 <- tidy(LR_digestive_CREBBP, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CREBBP_type_2, "pathvar_CREBBP_type_2.csv", row.names=FALSE)
system('dx upload pathvar_CREBBP_type_2.csv --path Emily-folder/results_2/pathvar_CREBBP_type_2.csv')



LR_respiratory_CREBBP <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CREBBP, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_CREBBP)
# exp(cbind(OR = coef(LR_respiratory_CREBBP), confint(LR_respiratory_CREBBP)))
pathvar_CREBBP_type_3 <- tidy(LR_respiratory_CREBBP, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CREBBP_type_3, "pathvar_CREBBP_type_3.csv", row.names=FALSE)
system('dx upload pathvar_CREBBP_type_3.csv --path Emily-folder/results_2/pathvar_CREBBP_type_3.csv')



LR_bone_CREBBP <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CREBBP, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_CREBBP)
# exp(cbind(OR = coef(LR_bone_CREBBP), confint(LR_bone_CREBBP)))
pathvar_CREBBP_type_4 <- tidy(LR_bone_CREBBP, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CREBBP_type_4, "pathvar_CREBBP_type_4.csv", row.names=FALSE)
system('dx upload pathvar_CREBBP_type_4.csv --path Emily-folder/results_2/pathvar_CREBBP_type_4.csv')



LR_skin_CREBBP <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CREBBP, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_CREBBP)
# exp(cbind(OR = coef(LR_skin_CREBBP), confint(LR_skin_CREBBP)))
pathvar_CREBBP_type_5 <- tidy(LR_skin_CREBBP, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CREBBP_type_5, "pathvar_CREBBP_type_5.csv", row.names=FALSE)
system('dx upload pathvar_CREBBP_type_5.csv --path Emily-folder/results_2/pathvar_CREBBP_type_5.csv')



LR_mesothelium_CREBBP <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CREBBP, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_CREBBP)
# exp(cbind(OR = coef(LR_mesothelium_CREBBP), confint(LR_mesothelium_CREBBP)))
pathvar_CREBBP_type_6 <- tidy(LR_mesothelium_CREBBP, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CREBBP_type_6, "pathvar_CREBBP_type_6.csv", row.names=FALSE)
system('dx upload pathvar_CREBBP_type_6.csv --path Emily-folder/results_2/pathvar_CREBBP_type_6.csv')



LR_breast_CREBBP <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CREBBP, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_CREBBP)
# exp(cbind(OR = coef(LR_breast_CREBBP), confint(LR_breast_CREBBP)))
pathvar_CREBBP_type_7 <- tidy(LR_breast_CREBBP, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CREBBP_type_7, "pathvar_CREBBP_type_7.csv", row.names=FALSE)
system('dx upload pathvar_CREBBP_type_7.csv --path Emily-folder/results_2/pathvar_CREBBP_type_7.csv')



LR_femalegenital_CREBBP <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CREBBP, 
                               data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_CREBBP)
# exp(cbind(OR = coef(LR_femalegenital_CREBBP), confint(LR_femalegenital_CREBBP)))
pathvar_CREBBP_type_8 <- tidy(LR_femalegenital_CREBBP, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CREBBP_type_8, "pathvar_CREBBP_type_8.csv", row.names=FALSE)
system('dx upload pathvar_CREBBP_type_8.csv --path Emily-folder/results_2/pathvar_CREBBP_type_8.csv')



LR_malegenital_CREBBP <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CREBBP, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_CREBBP)
# exp(cbind(OR = coef(LR_malegenital_CREBBP), confint(LR_malegenital_CREBBP)))
pathvar_CREBBP_type_9 <- tidy(LR_malegenital_CREBBP, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CREBBP_type_9, "pathvar_CREBBP_type_9.csv", row.names=FALSE)
system('dx upload pathvar_CREBBP_type_9.csv --path Emily-folder/results_2/pathvar_CREBBP_type_9.csv')



LR_urinary_CREBBP <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CREBBP, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_CREBBP)
# exp(cbind(OR = coef(LR_urinary_CREBBP), confint(LR_urinary_CREBBP)))
pathvar_CREBBP_type_10 <- tidy(LR_urinary_CREBBP, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CREBBP_type_10, "pathvar_CREBBP_type_10.csv", row.names=FALSE)
system('dx upload pathvar_CREBBP_type_10.csv --path Emily-folder/results_2/pathvar_CREBBP_type_10.csv')



LR_cns_CREBBP <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CREBBP, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_CREBBP)
# exp(cbind(OR = coef(LR_cns_CREBBP), confint(LR_cns_CREBBP)))
pathvar_CREBBP_type_11 <- tidy(LR_cns_CREBBP, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CREBBP_type_11, "pathvar_CREBBP_type_11.csv", row.names=FALSE)
system('dx upload pathvar_CREBBP_type_11.csv --path Emily-folder/results_2/pathvar_CREBBP_type_11.csv')



LR_endocrine_CREBBP <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CREBBP, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_CREBBP)
# exp(cbind(OR = coef(LR_endocrine_CREBBP), confint(LR_endocrine_CREBBP)))
pathvar_CREBBP_type_12 <- tidy(LR_endocrine_CREBBP, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CREBBP_type_12, "pathvar_CREBBP_type_12.csv", row.names=FALSE)
system('dx upload pathvar_CREBBP_type_12.csv --path Emily-folder/results_2/pathvar_CREBBP_type_12.csv')



LR_lymphatic_CREBBP <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CREBBP, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_CREBBP)
# exp(cbind(OR = coef(LR_lymphatic_CREBBP), confint(LR_lymphatic_CREBBP)))
pathvar_CREBBP_type_13 <- tidy(LR_lymphatic_CREBBP, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CREBBP_type_13, "pathvar_CREBBP_type_13.csv", row.names=FALSE)
system('dx upload pathvar_CREBBP_type_13.csv --path Emily-folder/results_2/pathvar_CREBBP_type_13.csv')



LR_secondary_CREBBP <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CREBBP, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_CREBBP)
# exp(cbind(OR = coef(LR_secondary_CREBBP), confint(LR_secondary_CREBBP)))
pathvar_CREBBP_type_14 <- tidy(LR_secondary_CREBBP, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CREBBP_type_14, "pathvar_CREBBP_type_14.csv", row.names=FALSE)
system('dx upload pathvar_CREBBP_type_14.csv --path Emily-folder/results_2/pathvar_CREBBP_type_14.csv')



LR_other_CREBBP <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CREBBP, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_CREBBP)
# exp(cbind(OR = coef(LR_other_CREBBP), confint(LR_other_CREBBP)))
pathvar_CREBBP_type_15 <- tidy(LR_other_CREBBP, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CREBBP_type_15, "pathvar_CREBBP_type_15.csv", row.names=FALSE)
system('dx upload pathvar_CREBBP_type_15.csv --path Emily-folder/results_2/pathvar_CREBBP_type_15.csv')










# TRIM37 Type ----------------------------------------------------------------- 





LR_oral_TRIM37 <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TRIM37, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_TRIM37)
# exp(cbind(OR = coef(LR_oral_TRIM37), confint(LR_oral_TRIM37)))
pathvar_TRIM37_type_1 <- tidy(LR_oral_TRIM37, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TRIM37_type_1, "pathvar_TRIM37_type_1.csv", row.names=FALSE)
system('dx upload pathvar_TRIM37_type_1.csv --path Emily-folder/results_2/pathvar_TRIM37_type_1.csv')



LR_digestive_TRIM37 <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TRIM37, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_TRIM37)
# exp(cbind(OR = coef(LR_digestive_TRIM37), confint(LR_digestive_TRIM37)))
pathvar_TRIM37_type_2 <- tidy(LR_digestive_TRIM37, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TRIM37_type_2, "pathvar_TRIM37_type_2.csv", row.names=FALSE)
system('dx upload pathvar_TRIM37_type_2.csv --path Emily-folder/results_2/pathvar_TRIM37_type_2.csv')



LR_respiratory_TRIM37 <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TRIM37, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_TRIM37)
# exp(cbind(OR = coef(LR_respiratory_TRIM37), confint(LR_respiratory_TRIM37)))
pathvar_TRIM37_type_3 <- tidy(LR_respiratory_TRIM37, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TRIM37_type_3, "pathvar_TRIM37_type_3.csv", row.names=FALSE)
system('dx upload pathvar_TRIM37_type_3.csv --path Emily-folder/results_2/pathvar_TRIM37_type_3.csv')



LR_bone_TRIM37 <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TRIM37, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_TRIM37)
# exp(cbind(OR = coef(LR_bone_TRIM37), confint(LR_bone_TRIM37)))
pathvar_TRIM37_type_4 <- tidy(LR_bone_TRIM37, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TRIM37_type_4, "pathvar_TRIM37_type_4.csv", row.names=FALSE)
system('dx upload pathvar_TRIM37_type_4.csv --path Emily-folder/results_2/pathvar_TRIM37_type_4.csv')



LR_skin_TRIM37 <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TRIM37, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_TRIM37)
# exp(cbind(OR = coef(LR_skin_TRIM37), confint(LR_skin_TRIM37)))
pathvar_TRIM37_type_5 <- tidy(LR_skin_TRIM37, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TRIM37_type_5, "pathvar_TRIM37_type_5.csv", row.names=FALSE)
system('dx upload pathvar_TRIM37_type_5.csv --path Emily-folder/results_2/pathvar_TRIM37_type_5.csv')



LR_mesothelium_TRIM37 <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TRIM37, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_TRIM37)
# exp(cbind(OR = coef(LR_mesothelium_TRIM37), confint(LR_mesothelium_TRIM37)))
pathvar_TRIM37_type_6 <- tidy(LR_mesothelium_TRIM37, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TRIM37_type_6, "pathvar_TRIM37_type_6.csv", row.names=FALSE)
system('dx upload pathvar_TRIM37_type_6.csv --path Emily-folder/results_2/pathvar_TRIM37_type_6.csv')



LR_breast_TRIM37 <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TRIM37, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_TRIM37)
# exp(cbind(OR = coef(LR_breast_TRIM37), confint(LR_breast_TRIM37)))
pathvar_TRIM37_type_7 <- tidy(LR_breast_TRIM37, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TRIM37_type_7, "pathvar_TRIM37_type_7.csv", row.names=FALSE)
system('dx upload pathvar_TRIM37_type_7.csv --path Emily-folder/results_2/pathvar_TRIM37_type_7.csv')



LR_femalegenital_TRIM37 <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TRIM37, 
                               data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_TRIM37)
# exp(cbind(OR = coef(LR_femalegenital_TRIM37), confint(LR_femalegenital_TRIM37)))
pathvar_TRIM37_type_8 <- tidy(LR_femalegenital_TRIM37, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TRIM37_type_8, "pathvar_TRIM37_type_8.csv", row.names=FALSE)
system('dx upload pathvar_TRIM37_type_8.csv --path Emily-folder/results_2/pathvar_TRIM37_type_8.csv')



LR_malegenital_TRIM37 <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TRIM37, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_TRIM37)
# exp(cbind(OR = coef(LR_malegenital_TRIM37), confint(LR_malegenital_TRIM37)))
pathvar_TRIM37_type_9 <- tidy(LR_malegenital_TRIM37, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TRIM37_type_9, "pathvar_TRIM37_type_9.csv", row.names=FALSE)
system('dx upload pathvar_TRIM37_type_9.csv --path Emily-folder/results_2/pathvar_TRIM37_type_9.csv')



LR_urinary_TRIM37 <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TRIM37, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_TRIM37)
# exp(cbind(OR = coef(LR_urinary_TRIM37), confint(LR_urinary_TRIM37)))
pathvar_TRIM37_type_10 <- tidy(LR_urinary_TRIM37, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TRIM37_type_10, "pathvar_TRIM37_type_10.csv", row.names=FALSE)
system('dx upload pathvar_TRIM37_type_10.csv --path Emily-folder/results_2/pathvar_TRIM37_type_10.csv')



LR_cns_TRIM37 <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TRIM37, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_TRIM37)
# exp(cbind(OR = coef(LR_cns_TRIM37), confint(LR_cns_TRIM37)))
pathvar_TRIM37_type_11 <- tidy(LR_cns_TRIM37, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TRIM37_type_11, "pathvar_TRIM37_type_11.csv", row.names=FALSE)
system('dx upload pathvar_TRIM37_type_11.csv --path Emily-folder/results_2/pathvar_TRIM37_type_11.csv')



LR_endocrine_TRIM37 <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TRIM37, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_TRIM37)
# exp(cbind(OR = coef(LR_endocrine_TRIM37), confint(LR_endocrine_TRIM37)))
pathvar_TRIM37_type_12 <- tidy(LR_endocrine_TRIM37, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TRIM37_type_12, "pathvar_TRIM37_type_12.csv", row.names=FALSE)
system('dx upload pathvar_TRIM37_type_12.csv --path Emily-folder/results_2/pathvar_TRIM37_type_12.csv')



LR_lymphatic_TRIM37 <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TRIM37, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_TRIM37)
# exp(cbind(OR = coef(LR_lymphatic_TRIM37), confint(LR_lymphatic_TRIM37)))
pathvar_TRIM37_type_13 <- tidy(LR_lymphatic_TRIM37, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TRIM37_type_13, "pathvar_TRIM37_type_13.csv", row.names=FALSE)
system('dx upload pathvar_TRIM37_type_13.csv --path Emily-folder/results_2/pathvar_TRIM37_type_13.csv')



LR_secondary_TRIM37 <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TRIM37, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_TRIM37)
# exp(cbind(OR = coef(LR_secondary_TRIM37), confint(LR_secondary_TRIM37)))
pathvar_TRIM37_type_14 <- tidy(LR_secondary_TRIM37, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TRIM37_type_14, "pathvar_TRIM37_type_14.csv", row.names=FALSE)
system('dx upload pathvar_TRIM37_type_14.csv --path Emily-folder/results_2/pathvar_TRIM37_type_14.csv')



LR_other_TRIM37 <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TRIM37, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_TRIM37)
# exp(cbind(OR = coef(LR_other_TRIM37), confint(LR_other_TRIM37)))
pathvar_TRIM37_type_15 <- tidy(LR_other_TRIM37, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TRIM37_type_15, "pathvar_TRIM37_type_15.csv", row.names=FALSE)
system('dx upload pathvar_TRIM37_type_15.csv --path Emily-folder/results_2/pathvar_TRIM37_type_15.csv')









# SMAD4 Type ----------------------------------------------------------------- 

# too few participants











# EP300 Type ----------------------------------------------------------------- 



LR_oral_EP300 <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+EP300, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_EP300)
# exp(cbind(OR = coef(LR_oral_EP300), confint(LR_oral_EP300)))
pathvar_EP300_type_1 <- tidy(LR_oral_EP300, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_EP300_type_1, "pathvar_EP300_type_1.csv", row.names=FALSE)
system('dx upload pathvar_EP300_type_1.csv --path Emily-folder/results_2/pathvar_EP300_type_1.csv')



LR_digestive_EP300 <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+EP300, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_EP300)
# exp(cbind(OR = coef(LR_digestive_EP300), confint(LR_digestive_EP300)))
pathvar_EP300_type_2 <- tidy(LR_digestive_EP300, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_EP300_type_2, "pathvar_EP300_type_2.csv", row.names=FALSE)
system('dx upload pathvar_EP300_type_2.csv --path Emily-folder/results_2/pathvar_EP300_type_2.csv')



LR_respiratory_EP300 <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+EP300, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_EP300)
# exp(cbind(OR = coef(LR_respiratory_EP300), confint(LR_respiratory_EP300)))
pathvar_EP300_type_3 <- tidy(LR_respiratory_EP300, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_EP300_type_3, "pathvar_EP300_type_3.csv", row.names=FALSE)
system('dx upload pathvar_EP300_type_3.csv --path Emily-folder/results_2/pathvar_EP300_type_3.csv')



LR_bone_EP300 <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+EP300, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_EP300)
# exp(cbind(OR = coef(LR_bone_EP300), confint(LR_bone_EP300)))
pathvar_EP300_type_4 <- tidy(LR_bone_EP300, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_EP300_type_4, "pathvar_EP300_type_4.csv", row.names=FALSE)
system('dx upload pathvar_EP300_type_4.csv --path Emily-folder/results_2/pathvar_EP300_type_4.csv')



LR_skin_EP300 <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+EP300, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_EP300)
# exp(cbind(OR = coef(LR_skin_EP300), confint(LR_skin_EP300)))
pathvar_EP300_type_5 <- tidy(LR_skin_EP300, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_EP300_type_5, "pathvar_EP300_type_5.csv", row.names=FALSE)
system('dx upload pathvar_EP300_type_5.csv --path Emily-folder/results_2/pathvar_EP300_type_5.csv')



LR_mesothelium_EP300 <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+EP300, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_EP300)
# exp(cbind(OR = coef(LR_mesothelium_EP300), confint(LR_mesothelium_EP300)))
pathvar_EP300_type_6 <- tidy(LR_mesothelium_EP300, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_EP300_type_6, "pathvar_EP300_type_6.csv", row.names=FALSE)
system('dx upload pathvar_EP300_type_6.csv --path Emily-folder/results_2/pathvar_EP300_type_6.csv')



LR_breast_EP300 <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+EP300, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_EP300)
# exp(cbind(OR = coef(LR_breast_EP300), confint(LR_breast_EP300)))
pathvar_EP300_type_7 <- tidy(LR_breast_EP300, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_EP300_type_7, "pathvar_EP300_type_7.csv", row.names=FALSE)
system('dx upload pathvar_EP300_type_7.csv --path Emily-folder/results_2/pathvar_EP300_type_7.csv')



LR_femalegenital_EP300 <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+EP300, 
                              data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_EP300)
# exp(cbind(OR = coef(LR_femalegenital_EP300), confint(LR_femalegenital_EP300)))
pathvar_EP300_type_8 <- tidy(LR_femalegenital_EP300, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_EP300_type_8, "pathvar_EP300_type_8.csv", row.names=FALSE)
system('dx upload pathvar_EP300_type_8.csv --path Emily-folder/results_2/pathvar_EP300_type_8.csv')



LR_malegenital_EP300 <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+EP300, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_EP300)
# exp(cbind(OR = coef(LR_malegenital_EP300), confint(LR_malegenital_EP300)))
pathvar_EP300_type_9 <- tidy(LR_malegenital_EP300, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_EP300_type_9, "pathvar_EP300_type_9.csv", row.names=FALSE)
system('dx upload pathvar_EP300_type_9.csv --path Emily-folder/results_2/pathvar_EP300_type_9.csv')



LR_urinary_EP300 <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+EP300, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_EP300)
# exp(cbind(OR = coef(LR_urinary_EP300), confint(LR_urinary_EP300)))
pathvar_EP300_type_10 <- tidy(LR_urinary_EP300, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_EP300_type_10, "pathvar_EP300_type_10.csv", row.names=FALSE)
system('dx upload pathvar_EP300_type_10.csv --path Emily-folder/results_2/pathvar_EP300_type_10.csv')



LR_cns_EP300 <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+EP300, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_EP300)
# exp(cbind(OR = coef(LR_cns_EP300), confint(LR_cns_EP300)))
pathvar_EP300_type_11 <- tidy(LR_cns_EP300, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_EP300_type_11, "pathvar_EP300_type_11.csv", row.names=FALSE)
system('dx upload pathvar_EP300_type_11.csv --path Emily-folder/results_2/pathvar_EP300_type_11.csv')



LR_endocrine_EP300 <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+EP300, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_EP300)
# exp(cbind(OR = coef(LR_endocrine_EP300), confint(LR_endocrine_EP300)))
pathvar_EP300_type_12 <- tidy(LR_endocrine_EP300, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_EP300_type_12, "pathvar_EP300_type_12.csv", row.names=FALSE)
system('dx upload pathvar_EP300_type_12.csv --path Emily-folder/results_2/pathvar_EP300_type_12.csv')



LR_lymphatic_EP300 <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+EP300, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_EP300)
# exp(cbind(OR = coef(LR_lymphatic_EP300), confint(LR_lymphatic_EP300)))
pathvar_EP300_type_13 <- tidy(LR_lymphatic_EP300, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_EP300_type_13, "pathvar_EP300_type_13.csv", row.names=FALSE)
system('dx upload pathvar_EP300_type_13.csv --path Emily-folder/results_2/pathvar_EP300_type_13.csv')



LR_secondary_EP300 <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+EP300, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_EP300)
# exp(cbind(OR = coef(LR_secondary_EP300), confint(LR_secondary_EP300)))
pathvar_EP300_type_14 <- tidy(LR_secondary_EP300, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_EP300_type_14, "pathvar_EP300_type_14.csv", row.names=FALSE)
system('dx upload pathvar_EP300_type_14.csv --path Emily-folder/results_2/pathvar_EP300_type_14.csv')



LR_other_EP300 <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+EP300, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_EP300)
# exp(cbind(OR = coef(LR_other_EP300), confint(LR_other_EP300)))
pathvar_EP300_type_15 <- tidy(LR_other_EP300, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_EP300_type_15, "pathvar_EP300_type_15.csv", row.names=FALSE)
system('dx upload pathvar_EP300_type_15.csv --path Emily-folder/results_2/pathvar_EP300_type_15.csv')










# LZTR1 Type ----------------------------------------------------------------- 



LR_oral_LZTR1 <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LZTR1, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_LZTR1)
# exp(cbind(OR = coef(LR_oral_LZTR1), confint(LR_oral_LZTR1)))
pathvar_LZTR1_type_1 <- tidy(LR_oral_LZTR1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LZTR1_type_1, "pathvar_LZTR1_type_1.csv", row.names=FALSE)
system('dx upload pathvar_LZTR1_type_1.csv --path Emily-folder/results_2/pathvar_LZTR1_type_1.csv')



LR_digestive_LZTR1 <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LZTR1, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_LZTR1)
# exp(cbind(OR = coef(LR_digestive_LZTR1), confint(LR_digestive_LZTR1)))
pathvar_LZTR1_type_2 <- tidy(LR_digestive_LZTR1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LZTR1_type_2, "pathvar_LZTR1_type_2.csv", row.names=FALSE)
system('dx upload pathvar_LZTR1_type_2.csv --path Emily-folder/results_2/pathvar_LZTR1_type_2.csv')



LR_respiratory_LZTR1 <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LZTR1, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_LZTR1)
# exp(cbind(OR = coef(LR_respiratory_LZTR1), confint(LR_respiratory_LZTR1)))
pathvar_LZTR1_type_3 <- tidy(LR_respiratory_LZTR1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LZTR1_type_3, "pathvar_LZTR1_type_3.csv", row.names=FALSE)
system('dx upload pathvar_LZTR1_type_3.csv --path Emily-folder/results_2/pathvar_LZTR1_type_3.csv')



LR_bone_LZTR1 <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LZTR1, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_LZTR1)
# exp(cbind(OR = coef(LR_bone_LZTR1), confint(LR_bone_LZTR1)))
pathvar_LZTR1_type_4 <- tidy(LR_bone_LZTR1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LZTR1_type_4, "pathvar_LZTR1_type_4.csv", row.names=FALSE)
system('dx upload pathvar_LZTR1_type_4.csv --path Emily-folder/results_2/pathvar_LZTR1_type_4.csv')



LR_skin_LZTR1 <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LZTR1, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_LZTR1)
# exp(cbind(OR = coef(LR_skin_LZTR1), confint(LR_skin_LZTR1)))
pathvar_LZTR1_type_5 <- tidy(LR_skin_LZTR1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LZTR1_type_5, "pathvar_LZTR1_type_5.csv", row.names=FALSE)
system('dx upload pathvar_LZTR1_type_5.csv --path Emily-folder/results_2/pathvar_LZTR1_type_5.csv')



LR_mesothelium_LZTR1 <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LZTR1, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_LZTR1)
# exp(cbind(OR = coef(LR_mesothelium_LZTR1), confint(LR_mesothelium_LZTR1)))
pathvar_LZTR1_type_6 <- tidy(LR_mesothelium_LZTR1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LZTR1_type_6, "pathvar_LZTR1_type_6.csv", row.names=FALSE)
system('dx upload pathvar_LZTR1_type_6.csv --path Emily-folder/results_2/pathvar_LZTR1_type_6.csv')



LR_breast_LZTR1 <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LZTR1, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_LZTR1)
# exp(cbind(OR = coef(LR_breast_LZTR1), confint(LR_breast_LZTR1)))
pathvar_LZTR1_type_7 <- tidy(LR_breast_LZTR1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LZTR1_type_7, "pathvar_LZTR1_type_7.csv", row.names=FALSE)
system('dx upload pathvar_LZTR1_type_7.csv --path Emily-folder/results_2/pathvar_LZTR1_type_7.csv')



LR_femalegenital_LZTR1 <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LZTR1, 
                              data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_LZTR1)
# exp(cbind(OR = coef(LR_femalegenital_LZTR1), confint(LR_femalegenital_LZTR1)))
pathvar_LZTR1_type_8 <- tidy(LR_femalegenital_LZTR1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LZTR1_type_8, "pathvar_LZTR1_type_8.csv", row.names=FALSE)
system('dx upload pathvar_LZTR1_type_8.csv --path Emily-folder/results_2/pathvar_LZTR1_type_8.csv')



LR_malegenital_LZTR1 <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LZTR1, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_LZTR1)
# exp(cbind(OR = coef(LR_malegenital_LZTR1), confint(LR_malegenital_LZTR1)))
pathvar_LZTR1_type_9 <- tidy(LR_malegenital_LZTR1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LZTR1_type_9, "pathvar_LZTR1_type_9.csv", row.names=FALSE)
system('dx upload pathvar_LZTR1_type_9.csv --path Emily-folder/results_2/pathvar_LZTR1_type_9.csv')



LR_urinary_LZTR1 <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LZTR1, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_LZTR1)
# exp(cbind(OR = coef(LR_urinary_LZTR1), confint(LR_urinary_LZTR1)))
pathvar_LZTR1_type_10 <- tidy(LR_urinary_LZTR1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LZTR1_type_10, "pathvar_LZTR1_type_10.csv", row.names=FALSE)
system('dx upload pathvar_LZTR1_type_10.csv --path Emily-folder/results_2/pathvar_LZTR1_type_10.csv')



LR_cns_LZTR1 <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LZTR1, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_LZTR1)
# exp(cbind(OR = coef(LR_cns_LZTR1), confint(LR_cns_LZTR1)))
pathvar_LZTR1_type_11 <- tidy(LR_cns_LZTR1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LZTR1_type_11, "pathvar_LZTR1_type_11.csv", row.names=FALSE)
system('dx upload pathvar_LZTR1_type_11.csv --path Emily-folder/results_2/pathvar_LZTR1_type_11.csv')



LR_endocrine_LZTR1 <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LZTR1, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_LZTR1)
# exp(cbind(OR = coef(LR_endocrine_LZTR1), confint(LR_endocrine_LZTR1)))
pathvar_LZTR1_type_12 <- tidy(LR_endocrine_LZTR1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LZTR1_type_12, "pathvar_LZTR1_type_12.csv", row.names=FALSE)
system('dx upload pathvar_LZTR1_type_12.csv --path Emily-folder/results_2/pathvar_LZTR1_type_12.csv')



LR_lymphatic_LZTR1 <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LZTR1, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_LZTR1)
# exp(cbind(OR = coef(LR_lymphatic_LZTR1), confint(LR_lymphatic_LZTR1)))
pathvar_LZTR1_type_13 <- tidy(LR_lymphatic_LZTR1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LZTR1_type_13, "pathvar_LZTR1_type_13.csv", row.names=FALSE)
system('dx upload pathvar_LZTR1_type_13.csv --path Emily-folder/results_2/pathvar_LZTR1_type_13.csv')



LR_secondary_LZTR1 <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LZTR1, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_LZTR1)
# exp(cbind(OR = coef(LR_secondary_LZTR1), confint(LR_secondary_LZTR1)))
pathvar_LZTR1_type_14 <- tidy(LR_secondary_LZTR1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LZTR1_type_14, "pathvar_LZTR1_type_14.csv", row.names=FALSE)
system('dx upload pathvar_LZTR1_type_14.csv --path Emily-folder/results_2/pathvar_LZTR1_type_14.csv')



LR_other_LZTR1 <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LZTR1, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_LZTR1)
# exp(cbind(OR = coef(LR_other_LZTR1), confint(LR_other_LZTR1)))
pathvar_LZTR1_type_15 <- tidy(LR_other_LZTR1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LZTR1_type_15, "pathvar_LZTR1_type_15.csv", row.names=FALSE)
system('dx upload pathvar_LZTR1_type_15.csv --path Emily-folder/results_2/pathvar_LZTR1_type_15.csv')










# KDM6A Type ----------------------------------------------------------------- 




LR_oral_KDM6A <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+KDM6A, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_KDM6A)
# exp(cbind(OR = coef(LR_oral_KDM6A), confint(LR_oral_KDM6A)))
pathvar_KDM6A_type_1 <- tidy(LR_oral_KDM6A, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_KDM6A_type_1, "pathvar_KDM6A_type_1.csv", row.names=FALSE)
system('dx upload pathvar_KDM6A_type_1.csv --path Emily-folder/results_2/pathvar_KDM6A_type_1.csv')



LR_digestive_KDM6A <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+KDM6A, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_KDM6A)
# exp(cbind(OR = coef(LR_digestive_KDM6A), confint(LR_digestive_KDM6A)))
pathvar_KDM6A_type_2 <- tidy(LR_digestive_KDM6A, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_KDM6A_type_2, "pathvar_KDM6A_type_2.csv", row.names=FALSE)
system('dx upload pathvar_KDM6A_type_2.csv --path Emily-folder/results_2/pathvar_KDM6A_type_2.csv')



LR_respiratory_KDM6A <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+KDM6A, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_KDM6A)
# exp(cbind(OR = coef(LR_respiratory_KDM6A), confint(LR_respiratory_KDM6A)))
pathvar_KDM6A_type_3 <- tidy(LR_respiratory_KDM6A, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_KDM6A_type_3, "pathvar_KDM6A_type_3.csv", row.names=FALSE)
system('dx upload pathvar_KDM6A_type_3.csv --path Emily-folder/results_2/pathvar_KDM6A_type_3.csv')



LR_bone_KDM6A <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+KDM6A, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_KDM6A)
# exp(cbind(OR = coef(LR_bone_KDM6A), confint(LR_bone_KDM6A)))
pathvar_KDM6A_type_4 <- tidy(LR_bone_KDM6A, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_KDM6A_type_4, "pathvar_KDM6A_type_4.csv", row.names=FALSE)
system('dx upload pathvar_KDM6A_type_4.csv --path Emily-folder/results_2/pathvar_KDM6A_type_4.csv')



LR_skin_KDM6A <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+KDM6A, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_KDM6A)
# exp(cbind(OR = coef(LR_skin_KDM6A), confint(LR_skin_KDM6A)))
pathvar_KDM6A_type_5 <- tidy(LR_skin_KDM6A, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_KDM6A_type_5, "pathvar_KDM6A_type_5.csv", row.names=FALSE)
system('dx upload pathvar_KDM6A_type_5.csv --path Emily-folder/results_2/pathvar_KDM6A_type_5.csv')



LR_mesothelium_KDM6A <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+KDM6A, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_KDM6A)
# exp(cbind(OR = coef(LR_mesothelium_KDM6A), confint(LR_mesothelium_KDM6A)))
pathvar_KDM6A_type_6 <- tidy(LR_mesothelium_KDM6A, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_KDM6A_type_6, "pathvar_KDM6A_type_6.csv", row.names=FALSE)
system('dx upload pathvar_KDM6A_type_6.csv --path Emily-folder/results_2/pathvar_KDM6A_type_6.csv')



LR_breast_KDM6A <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+KDM6A, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_KDM6A)
# exp(cbind(OR = coef(LR_breast_KDM6A), confint(LR_breast_KDM6A)))
pathvar_KDM6A_type_7 <- tidy(LR_breast_KDM6A, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_KDM6A_type_7, "pathvar_KDM6A_type_7.csv", row.names=FALSE)
system('dx upload pathvar_KDM6A_type_7.csv --path Emily-folder/results_2/pathvar_KDM6A_type_7.csv')



LR_femalegenital_KDM6A <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+KDM6A, 
                              data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_KDM6A)
# exp(cbind(OR = coef(LR_femalegenital_KDM6A), confint(LR_femalegenital_KDM6A)))
pathvar_KDM6A_type_8 <- tidy(LR_femalegenital_KDM6A, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_KDM6A_type_8, "pathvar_KDM6A_type_8.csv", row.names=FALSE)
system('dx upload pathvar_KDM6A_type_8.csv --path Emily-folder/results_2/pathvar_KDM6A_type_8.csv')



LR_malegenital_KDM6A <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+KDM6A, 
                            data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_KDM6A)
# exp(cbind(OR = coef(LR_malegenital_KDM6A), confint(LR_malegenital_KDM6A)))
pathvar_KDM6A_type_9 <- tidy(LR_malegenital_KDM6A, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_KDM6A_type_9, "pathvar_KDM6A_type_9.csv", row.names=FALSE)
system('dx upload pathvar_KDM6A_type_9.csv --path Emily-folder/results_2/pathvar_KDM6A_type_9.csv')



LR_urinary_KDM6A <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+KDM6A, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_KDM6A)
# exp(cbind(OR = coef(LR_urinary_KDM6A), confint(LR_urinary_KDM6A)))
pathvar_KDM6A_type_10 <- tidy(LR_urinary_KDM6A, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_KDM6A_type_10, "pathvar_KDM6A_type_10.csv", row.names=FALSE)
system('dx upload pathvar_KDM6A_type_10.csv --path Emily-folder/results_2/pathvar_KDM6A_type_10.csv')



LR_cns_KDM6A <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+KDM6A, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_KDM6A)
# exp(cbind(OR = coef(LR_cns_KDM6A), confint(LR_cns_KDM6A)))
pathvar_KDM6A_type_11 <- tidy(LR_cns_KDM6A, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_KDM6A_type_11, "pathvar_KDM6A_type_11.csv", row.names=FALSE)
system('dx upload pathvar_KDM6A_type_11.csv --path Emily-folder/results_2/pathvar_KDM6A_type_11.csv')



LR_endocrine_KDM6A <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+KDM6A, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_KDM6A)
# exp(cbind(OR = coef(LR_endocrine_KDM6A), confint(LR_endocrine_KDM6A)))
pathvar_KDM6A_type_12 <- tidy(LR_endocrine_KDM6A, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_KDM6A_type_12, "pathvar_KDM6A_type_12.csv", row.names=FALSE)
system('dx upload pathvar_KDM6A_type_12.csv --path Emily-folder/results_2/pathvar_KDM6A_type_12.csv')



LR_lymphatic_KDM6A <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+KDM6A, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_KDM6A)
# exp(cbind(OR = coef(LR_lymphatic_KDM6A), confint(LR_lymphatic_KDM6A)))
pathvar_KDM6A_type_13 <- tidy(LR_lymphatic_KDM6A, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_KDM6A_type_13, "pathvar_KDM6A_type_13.csv", row.names=FALSE)
system('dx upload pathvar_KDM6A_type_13.csv --path Emily-folder/results_2/pathvar_KDM6A_type_13.csv')



LR_secondary_KDM6A <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+KDM6A, 
                          data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_KDM6A)
# exp(cbind(OR = coef(LR_secondary_KDM6A), confint(LR_secondary_KDM6A)))
pathvar_KDM6A_type_14 <- tidy(LR_secondary_KDM6A, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_KDM6A_type_14, "pathvar_KDM6A_type_14.csv", row.names=FALSE)
system('dx upload pathvar_KDM6A_type_14.csv --path Emily-folder/results_2/pathvar_KDM6A_type_14.csv')



LR_other_KDM6A <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+KDM6A, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_KDM6A)
# exp(cbind(OR = coef(LR_other_KDM6A), confint(LR_other_KDM6A)))
pathvar_KDM6A_type_15 <- tidy(LR_other_KDM6A, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_KDM6A_type_15, "pathvar_KDM6A_type_15.csv", row.names=FALSE)
system('dx upload pathvar_KDM6A_type_15.csv --path Emily-folder/results_2/pathvar_KDM6A_type_15.csv')











# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_CAgenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_CAgenes.R')
