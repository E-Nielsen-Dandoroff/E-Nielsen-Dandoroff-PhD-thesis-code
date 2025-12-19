# EntireCohort_Diagnoses

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results
install.packages('broom')
library(broom)


# EntireCohort - Diagnoses ------------------------------------------------



EntireCohort_Diagnoses <- fread('/mnt/project/EntireCohort_Diagnoses.csv', data.table = FALSE)
EntireCohort_Diagnoses_ICD10 <- EntireCohort_Diagnoses[,-3]

# Q02 Microcephaly
# 
# Q04 Other congenital malformations of brain
# 
# Q77 Osteochondrodysplasia with defects of growth of tubular bones and spine
# 
# Q78 Other osteochondrodysplasias
# 
# Q79.6 Ehlers-Danlos syndrome
# Q79.8 Other congenital malformations of musculoskeletal system
# Q79.9 Congenital malformation of musculoskeletal system, unspecified
# 
# Q87.1 Congenital malformation syndromes predominantly associated with short stature
# Q87.2 Congenital malformation syndromes predominantly involving limbs
# 
# Q96 Turner's syndrome

SS_M_diagnoses <- grepl("Q02|Q04|Q77|Q78|Q79.6|Q79.8|Q79.9|Q87.1|Q87.2|Q96", EntireCohort_Diagnoses_ICD10$p41270)

EntireCohort_Diagnoses_ICD10$logical <- SS_M_diagnoses

EntireCohort_Diagnoses_ICD10 <- EntireCohort_Diagnoses_ICD10[EntireCohort_Diagnoses_ICD10$logical == "TRUE",]
# 531 participants with SS/M type diagnoses



# MGSvar ----------------------------

# Using data.table to load in the EntireCohort_MGS dataset
EntireCohort_MGSvar <- fread('/mnt/project/Emily-folder/MGSvar/EntireCohort_MGSvar.csv', data.table = FALSE)

# The dataframe gets loaded in with an extra column called 'V1' which I don't want
EntireCohort_MGSvar <- EntireCohort_MGSvar[-1]
# removed :)







# CAgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

# go thr all the genes to check that they're present
grep("XRCC4", EntireCohort_pathvar$GENE)
# all seem to be present exccept for RAF1

# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$CAgene <- NA
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE1 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1','THRB','TNFRSF11B','TRIM37','XRCC4')] <- 1
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE2 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1','THRB','TNFRSF11B','TRIM37','XRCC4')] <- 1
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE3 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1','THRB','TNFRSF11B','TRIM37','XRCC4')] <- 1
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE4 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1','THRB','TNFRSF11B','TRIM37','XRCC4')] <- 1
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE5 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1','THRB','TNFRSF11B','TRIM37','XRCC4')] <- 1
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE6 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1','THRB','TNFRSF11B','TRIM37','XRCC4')] <- 1
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE7 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1','THRB','TNFRSF11B','TRIM37','XRCC4')] <- 1
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE8 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1','THRB','TNFRSF11B','TRIM37','XRCC4')] <- 1
EntireCohort_pathvar$CAgene[EntireCohort_pathvar$GENE9 %in% c('BRAF','BUB1B','CBL','CREBBP','EP300','ERCC3','FANCA','KDM6A','KRAS','LIG4','LZTR1','NRAS','POLE','PTPN11','RAF1','RECQL3','RECQL4','RIT1','RRAS2','SBDS','SMAD4','SMARCAL1','SOS1','TALDO1','THRB','TNFRSF11B','TRIM37','XRCC4')] <- 1


EntireCohort_pathvar$CAgene[is.na(EntireCohort_pathvar$CAgene)] <- 0
EntireCohort_pathvar$CAgene <- as.numeric(EntireCohort_pathvar$CAgene)
unique(EntireCohort_pathvar$CAgene)
# 0 1

table(EntireCohort_pathvar$CAgene)
#      0      1 
# 436335  33464




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




### 'THRB'  
EntireCohort_pathvar$THRB <- NA
EntireCohort_pathvar$THRB[EntireCohort_pathvar$GENE1 %in% c("THRB")] <- 1
EntireCohort_pathvar$THRB[EntireCohort_pathvar$GENE2 %in% c("THRB")] <- 1
EntireCohort_pathvar$THRB[EntireCohort_pathvar$GENE3 %in% c("THRB")] <- 1
EntireCohort_pathvar$THRB[EntireCohort_pathvar$GENE4 %in% c("THRB")] <- 1
EntireCohort_pathvar$THRB[EntireCohort_pathvar$GENE5 %in% c("THRB")] <- 1
EntireCohort_pathvar$THRB[EntireCohort_pathvar$GENE6 %in% c("THRB")] <- 1
EntireCohort_pathvar$THRB[EntireCohort_pathvar$GENE7 %in% c("THRB")] <- 1
EntireCohort_pathvar$THRB[EntireCohort_pathvar$GENE8 %in% c("THRB")] <- 1
EntireCohort_pathvar$THRB[EntireCohort_pathvar$GENE9 %in% c("THRB")] <- 1

EntireCohort_pathvar$THRB[is.na(EntireCohort_pathvar$THRB)] <- 0
EntireCohort_pathvar$THRB <- as.numeric(EntireCohort_pathvar$THRB)
unique(EntireCohort_pathvar$THRB)
# 0 1

table(EntireCohort_pathvar$THRB)
#      0      1 
# 469790      9




### 'TNFRSF11B'
EntireCohort_pathvar$TNFRSF11B <- NA
EntireCohort_pathvar$TNFRSF11B[EntireCohort_pathvar$GENE1 %in% c("TNFRSF11B")] <- 1
EntireCohort_pathvar$TNFRSF11B[EntireCohort_pathvar$GENE2 %in% c("TNFRSF11B")] <- 1
EntireCohort_pathvar$TNFRSF11B[EntireCohort_pathvar$GENE3 %in% c("TNFRSF11B")] <- 1
EntireCohort_pathvar$TNFRSF11B[EntireCohort_pathvar$GENE4 %in% c("TNFRSF11B")] <- 1
EntireCohort_pathvar$TNFRSF11B[EntireCohort_pathvar$GENE5 %in% c("TNFRSF11B")] <- 1
EntireCohort_pathvar$TNFRSF11B[EntireCohort_pathvar$GENE6 %in% c("TNFRSF11B")] <- 1
EntireCohort_pathvar$TNFRSF11B[EntireCohort_pathvar$GENE7 %in% c("TNFRSF11B")] <- 1
EntireCohort_pathvar$TNFRSF11B[EntireCohort_pathvar$GENE8 %in% c("TNFRSF11B")] <- 1
EntireCohort_pathvar$TNFRSF11B[EntireCohort_pathvar$GENE9 %in% c("TNFRSF11B")] <- 1

EntireCohort_pathvar$TNFRSF11B[is.na(EntireCohort_pathvar$TNFRSF11B)] <- 0
EntireCohort_pathvar$TNFRSF11B <- as.numeric(EntireCohort_pathvar$TNFRSF11B)
unique(EntireCohort_pathvar$TNFRSF11B)
# 0 1

table(EntireCohort_pathvar$TNFRSF11B)
#      0      1 
# 469521    278





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





# Biosynthgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

# go thr all the genes to check that they're present
grep("NEPRO", EntireCohort_pathvar$GENE)
# all seem to be present, except for:
# PTDSS1 = 0
# PEX5 = 0



# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$Biosynthgene <- NA
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE1 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE2 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE3 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE4 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE5 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE6 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE7 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE8 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE9 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1


EntireCohort_pathvar$Biosynthgene[is.na(EntireCohort_pathvar$Biosynthgene)] <- 0
EntireCohort_pathvar$Biosynthgene <- as.numeric(EntireCohort_pathvar$Biosynthgene)
unique(EntireCohort_pathvar$Biosynthgene)
# 0 1 

table(EntireCohort_pathvar$Biosynthgene)
# 0      1 
# 405291  64508 




# now to make a column for each of the Biosynthgenes
### 'INSR'
EntireCohort_pathvar$INSR <- NA
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE1 %in% c("INSR")] <- 1
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE2 %in% c("INSR")] <- 1
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE3 %in% c("INSR")] <- 1
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE4 %in% c("INSR")] <- 1
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE5 %in% c("INSR")] <- 1
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE6 %in% c("INSR")] <- 1
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE7 %in% c("INSR")] <- 1
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE8 %in% c("INSR")] <- 1
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE9 %in% c("INSR")] <- 1

EntireCohort_pathvar$INSR[is.na(EntireCohort_pathvar$INSR)] <- 0
EntireCohort_pathvar$INSR <- as.numeric(EntireCohort_pathvar$INSR)
unique(EntireCohort_pathvar$INSR)
# 0 1

table(EntireCohort_pathvar$INSR)
#  0      1 
# 468326   1473





### 'PIK3R1'
EntireCohort_pathvar$PIK3R1 <- NA
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE1 %in% c("PIK3R1")] <- 1
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE2 %in% c("PIK3R1")] <- 1
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE3 %in% c("PIK3R1")] <- 1
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE4 %in% c("PIK3R1")] <- 1
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE5 %in% c("PIK3R1")] <- 1
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE6 %in% c("PIK3R1")] <- 1
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE7 %in% c("PIK3R1")] <- 1
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE8 %in% c("PIK3R1")] <- 1
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE9 %in% c("PIK3R1")] <- 1

EntireCohort_pathvar$PIK3R1[is.na(EntireCohort_pathvar$PIK3R1)] <- 0
EntireCohort_pathvar$PIK3R1 <- as.numeric(EntireCohort_pathvar$PIK3R1)
unique(EntireCohort_pathvar$PIK3R1)
# 0 1

table(EntireCohort_pathvar$PIK3R1)
#  0      1 
# 469410    389




### 'PRMT7'
EntireCohort_pathvar$PRMT7 <- NA
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE1 %in% c("PRMT7")] <- 1
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE2 %in% c("PRMT7")] <- 1
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE3 %in% c("PRMT7")] <- 1
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE4 %in% c("PRMT7")] <- 1
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE5 %in% c("PRMT7")] <- 1
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE6 %in% c("PRMT7")] <- 1
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE7 %in% c("PRMT7")] <- 1
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE8 %in% c("PRMT7")] <- 1
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE9 %in% c("PRMT7")] <- 1

EntireCohort_pathvar$PRMT7[is.na(EntireCohort_pathvar$PRMT7)] <- 0
EntireCohort_pathvar$PRMT7 <- as.numeric(EntireCohort_pathvar$PRMT7)
unique(EntireCohort_pathvar$PRMT7)
# 0 1

table(EntireCohort_pathvar$PRMT7)
#  0      1 
# 468091   1708




### 'PEX7'
EntireCohort_pathvar$PEX7 <- NA
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE1 %in% c("PEX7")] <- 1
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE2 %in% c("PEX7")] <- 1
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE3 %in% c("PEX7")] <- 1
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE4 %in% c("PEX7")] <- 1
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE5 %in% c("PEX7")] <- 1
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE6 %in% c("PEX7")] <- 1
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE7 %in% c("PEX7")] <- 1
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE8 %in% c("PEX7")] <- 1
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE9 %in% c("PEX7")] <- 1

EntireCohort_pathvar$PEX7[is.na(EntireCohort_pathvar$PEX7)] <- 0
EntireCohort_pathvar$PEX7 <- as.numeric(EntireCohort_pathvar$PEX7)
unique(EntireCohort_pathvar$PEX7)
# 0 1

table(EntireCohort_pathvar$PEX7)
#  0      1 
# 468124   1675 




### 'SMPD4'
EntireCohort_pathvar$SMPD4 <- NA
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE1 %in% c("SMPD4")] <- 1
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE2 %in% c("SMPD4")] <- 1
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE3 %in% c("SMPD4")] <- 1
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE4 %in% c("SMPD4")] <- 1
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE5 %in% c("SMPD4")] <- 1
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE6 %in% c("SMPD4")] <- 1
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE7 %in% c("SMPD4")] <- 1
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE8 %in% c("SMPD4")] <- 1
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE9 %in% c("SMPD4")] <- 1

EntireCohort_pathvar$SMPD4[is.na(EntireCohort_pathvar$SMPD4)] <- 0
EntireCohort_pathvar$SMPD4 <- as.numeric(EntireCohort_pathvar$SMPD4)
unique(EntireCohort_pathvar$SMPD4)
# 0 1

table(EntireCohort_pathvar$SMPD4)
#  0      1 
# 456726  13073 




### 'PAPSS2'
EntireCohort_pathvar$PAPSS2 <- NA
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE1 %in% c("PAPSS2")] <- 1
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE2 %in% c("PAPSS2")] <- 1
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE3 %in% c("PAPSS2")] <- 1
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE4 %in% c("PAPSS2")] <- 1
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE5 %in% c("PAPSS2")] <- 1
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE6 %in% c("PAPSS2")] <- 1
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE7 %in% c("PAPSS2")] <- 1
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE8 %in% c("PAPSS2")] <- 1
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE9 %in% c("PAPSS2")] <- 1

EntireCohort_pathvar$PAPSS2[is.na(EntireCohort_pathvar$PAPSS2)] <- 0
EntireCohort_pathvar$PAPSS2 <- as.numeric(EntireCohort_pathvar$PAPSS2)
unique(EntireCohort_pathvar$PAPSS2)
# 0 1

table(EntireCohort_pathvar$PAPSS2)
#  0      1 
# 467166   2633 




### 'B3GLCT'    
EntireCohort_pathvar$B3GLCT <- NA
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE1 %in% c("B3GLCT")] <- 1
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE2 %in% c("B3GLCT")] <- 1
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE3 %in% c("B3GLCT")] <- 1
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE4 %in% c("B3GLCT")] <- 1
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE5 %in% c("B3GLCT")] <- 1
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE6 %in% c("B3GLCT")] <- 1
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE7 %in% c("B3GLCT")] <- 1
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE8 %in% c("B3GLCT")] <- 1
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE9 %in% c("B3GLCT")] <- 1

EntireCohort_pathvar$B3GLCT[is.na(EntireCohort_pathvar$B3GLCT)] <- 0
EntireCohort_pathvar$B3GLCT <- as.numeric(EntireCohort_pathvar$B3GLCT)
unique(EntireCohort_pathvar$B3GLCT)
# 0 1

table(EntireCohort_pathvar$B3GLCT)
#  0      1 
# 467332   2467 




### 'MINPP1'
EntireCohort_pathvar$MINPP1 <- NA
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE1 %in% c("MINPP1")] <- 1
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE2 %in% c("MINPP1")] <- 1
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE3 %in% c("MINPP1")] <- 1
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE4 %in% c("MINPP1")] <- 1
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE5 %in% c("MINPP1")] <- 1
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE6 %in% c("MINPP1")] <- 1
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE7 %in% c("MINPP1")] <- 1
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE8 %in% c("MINPP1")] <- 1
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE9 %in% c("MINPP1")] <- 1

EntireCohort_pathvar$MINPP1[is.na(EntireCohort_pathvar$MINPP1)] <- 0
EntireCohort_pathvar$MINPP1 <- as.numeric(EntireCohort_pathvar$MINPP1)
unique(EntireCohort_pathvar$MINPP1)
# 0 1

table(EntireCohort_pathvar$MINPP1)
#  0      1 
# 469166    633




### 'BPNT2'
EntireCohort_pathvar$BPNT2 <- NA
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE1 %in% c("BPNT2")] <- 1
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE2 %in% c("BPNT2")] <- 1
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE3 %in% c("BPNT2")] <- 1
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE4 %in% c("BPNT2")] <- 1
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE5 %in% c("BPNT2")] <- 1
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE6 %in% c("BPNT2")] <- 1
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE7 %in% c("BPNT2")] <- 1
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE8 %in% c("BPNT2")] <- 1
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE9 %in% c("BPNT2")] <- 1

EntireCohort_pathvar$BPNT2[is.na(EntireCohort_pathvar$BPNT2)] <- 0
EntireCohort_pathvar$BPNT2 <- as.numeric(EntireCohort_pathvar$BPNT2)
unique(EntireCohort_pathvar$BPNT2)
# 0 1

table(EntireCohort_pathvar$BPNT2)
#  0      1 
# 469298    501 





### 'TMEM165'
EntireCohort_pathvar$TMEM165 <- NA
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE1 %in% c("TMEM165")] <- 1
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE2 %in% c("TMEM165")] <- 1
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE3 %in% c("TMEM165")] <- 1
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE4 %in% c("TMEM165")] <- 1
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE5 %in% c("TMEM165")] <- 1
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE6 %in% c("TMEM165")] <- 1
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE7 %in% c("TMEM165")] <- 1
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE8 %in% c("TMEM165")] <- 1
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE9 %in% c("TMEM165")] <- 1

EntireCohort_pathvar$TMEM165[is.na(EntireCohort_pathvar$TMEM165)] <- 0
EntireCohort_pathvar$TMEM165 <- as.numeric(EntireCohort_pathvar$TMEM165)
unique(EntireCohort_pathvar$TMEM165)
# 0 1

table(EntireCohort_pathvar$TMEM165)
#  0      1 
# 469329    470




### 'ALDH18A1'
EntireCohort_pathvar$ALDH18A1 <- NA
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE1 %in% c("ALDH18A1")] <- 1
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE2 %in% c("ALDH18A1")] <- 1
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE3 %in% c("ALDH18A1")] <- 1
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE4 %in% c("ALDH18A1")] <- 1
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE5 %in% c("ALDH18A1")] <- 1
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE6 %in% c("ALDH18A1")] <- 1
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE7 %in% c("ALDH18A1")] <- 1
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE8 %in% c("ALDH18A1")] <- 1
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE9 %in% c("ALDH18A1")] <- 1

EntireCohort_pathvar$ALDH18A1[is.na(EntireCohort_pathvar$ALDH18A1)] <- 0
EntireCohort_pathvar$ALDH18A1 <- as.numeric(EntireCohort_pathvar$ALDH18A1)
unique(EntireCohort_pathvar$ALDH18A1)
# 0 1

table(EntireCohort_pathvar$ALDH18A1)
#  0      1 
# 468481   1318 




### 'COASY'
EntireCohort_pathvar$COASY <- NA
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE1 %in% c("COASY")] <- 1
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE2 %in% c("COASY")] <- 1
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE3 %in% c("COASY")] <- 1
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE4 %in% c("COASY")] <- 1
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE5 %in% c("COASY")] <- 1
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE6 %in% c("COASY")] <- 1
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE7 %in% c("COASY")] <- 1
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE8 %in% c("COASY")] <- 1
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE9 %in% c("COASY")] <- 1

EntireCohort_pathvar$COASY[is.na(EntireCohort_pathvar$COASY)] <- 0
EntireCohort_pathvar$COASY <- as.numeric(EntireCohort_pathvar$COASY)
unique(EntireCohort_pathvar$COASY)
# 0 1

table(EntireCohort_pathvar$COASY)
#  0      1 
# 468634   1165 




### 'GNPAT'
EntireCohort_pathvar$GNPAT <- NA
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE1 %in% c("GNPAT")] <- 1
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE2 %in% c("GNPAT")] <- 1
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE3 %in% c("GNPAT")] <- 1
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE4 %in% c("GNPAT")] <- 1
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE5 %in% c("GNPAT")] <- 1
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE6 %in% c("GNPAT")] <- 1
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE7 %in% c("GNPAT")] <- 1
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE8 %in% c("GNPAT")] <- 1
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE9 %in% c("GNPAT")] <- 1

EntireCohort_pathvar$GNPAT[is.na(EntireCohort_pathvar$GNPAT)] <- 0
EntireCohort_pathvar$GNPAT <- as.numeric(EntireCohort_pathvar$GNPAT)
unique(EntireCohort_pathvar$GNPAT)
# 0 1

table(EntireCohort_pathvar$GNPAT)
#  0      1 
# 468842    957 




### 'PTDSS1'
# no variants present




### 'MSMO1'
EntireCohort_pathvar$MSMO1 <- NA
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE1 %in% c("MSMO1")] <- 1
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE2 %in% c("MSMO1")] <- 1
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE3 %in% c("MSMO1")] <- 1
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE4 %in% c("MSMO1")] <- 1
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE5 %in% c("MSMO1")] <- 1
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE6 %in% c("MSMO1")] <- 1
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE7 %in% c("MSMO1")] <- 1
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE8 %in% c("MSMO1")] <- 1
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE9 %in% c("MSMO1")] <- 1

EntireCohort_pathvar$MSMO1[is.na(EntireCohort_pathvar$MSMO1)] <- 0
EntireCohort_pathvar$MSMO1 <- as.numeric(EntireCohort_pathvar$MSMO1)
unique(EntireCohort_pathvar$MSMO1)
# 0 1

table(EntireCohort_pathvar$MSMO1)
#  0      1 
# 469177    622 





### 'DHCR24'
EntireCohort_pathvar$DHCR24 <- NA
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE1 %in% c("DHCR24")] <- 1
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE2 %in% c("DHCR24")] <- 1
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE3 %in% c("DHCR24")] <- 1
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE4 %in% c("DHCR24")] <- 1
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE5 %in% c("DHCR24")] <- 1
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE6 %in% c("DHCR24")] <- 1
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE7 %in% c("DHCR24")] <- 1
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE8 %in% c("DHCR24")] <- 1
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE9 %in% c("DHCR24")] <- 1

EntireCohort_pathvar$DHCR24[is.na(EntireCohort_pathvar$DHCR24)] <- 0
EntireCohort_pathvar$DHCR24 <- as.numeric(EntireCohort_pathvar$DHCR24)
unique(EntireCohort_pathvar$DHCR24)
# 0 1

table(EntireCohort_pathvar$DHCR24)
#  0      1 
# 466322   3477 





### 'DHCR7'
EntireCohort_pathvar$DHCR7 <- NA
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE1 %in% c("DHCR7")] <- 1
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE2 %in% c("DHCR7")] <- 1
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE3 %in% c("DHCR7")] <- 1
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE4 %in% c("DHCR7")] <- 1
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE5 %in% c("DHCR7")] <- 1
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE6 %in% c("DHCR7")] <- 1
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE7 %in% c("DHCR7")] <- 1
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE8 %in% c("DHCR7")] <- 1
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE9 %in% c("DHCR7")] <- 1

EntireCohort_pathvar$DHCR7[is.na(EntireCohort_pathvar$DHCR7)] <- 0
EntireCohort_pathvar$DHCR7 <- as.numeric(EntireCohort_pathvar$DHCR7)
unique(EntireCohort_pathvar$DHCR7)
# 0 1

table(EntireCohort_pathvar$DHCR7)
#  0      1 
# 457683  12116 




### 'PCYT1A'
EntireCohort_pathvar$PCYT1A <- NA
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE1 %in% c("PCYT1A")] <- 1
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE2 %in% c("PCYT1A")] <- 1
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE3 %in% c("PCYT1A")] <- 1
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE4 %in% c("PCYT1A")] <- 1
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE5 %in% c("PCYT1A")] <- 1
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE6 %in% c("PCYT1A")] <- 1
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE7 %in% c("PCYT1A")] <- 1
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE8 %in% c("PCYT1A")] <- 1
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE9 %in% c("PCYT1A")] <- 1

EntireCohort_pathvar$PCYT1A[is.na(EntireCohort_pathvar$PCYT1A)] <- 0
EntireCohort_pathvar$PCYT1A <- as.numeric(EntireCohort_pathvar$PCYT1A)
unique(EntireCohort_pathvar$PCYT1A)
# 0 1

table(EntireCohort_pathvar$PCYT1A)
#  0      1 
# 468548   1251 




### 'MTHFS'
EntireCohort_pathvar$MTHFS <- NA
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE1 %in% c("MTHFS")] <- 1
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE2 %in% c("MTHFS")] <- 1
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE3 %in% c("MTHFS")] <- 1
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE4 %in% c("MTHFS")] <- 1
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE5 %in% c("MTHFS")] <- 1
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE6 %in% c("MTHFS")] <- 1
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE7 %in% c("MTHFS")] <- 1
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE8 %in% c("MTHFS")] <- 1
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE9 %in% c("MTHFS")] <- 1

EntireCohort_pathvar$MTHFS[is.na(EntireCohort_pathvar$MTHFS)] <- 0
EntireCohort_pathvar$MTHFS <- as.numeric(EntireCohort_pathvar$MTHFS)
unique(EntireCohort_pathvar$MTHFS)
# 0 1

table(EntireCohort_pathvar$MTHFS)
#  0      1 
# 469150    649 




### 'NANS'
EntireCohort_pathvar$NANS <- NA
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE1 %in% c("NANS")] <- 1
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE2 %in% c("NANS")] <- 1
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE3 %in% c("NANS")] <- 1
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE4 %in% c("NANS")] <- 1
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE5 %in% c("NANS")] <- 1
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE6 %in% c("NANS")] <- 1
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE7 %in% c("NANS")] <- 1
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE8 %in% c("NANS")] <- 1
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE9 %in% c("NANS")] <- 1

EntireCohort_pathvar$NANS[is.na(EntireCohort_pathvar$NANS)] <- 0
EntireCohort_pathvar$NANS <- as.numeric(EntireCohort_pathvar$NANS)
unique(EntireCohort_pathvar$NANS)
# 0 1

table(EntireCohort_pathvar$NANS)
#  0      1 
# 467588   2211




### 'EBP'
EntireCohort_pathvar$EBP <- NA
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE1 %in% c("EBP")] <- 1
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE2 %in% c("EBP")] <- 1
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE3 %in% c("EBP")] <- 1
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE4 %in% c("EBP")] <- 1
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE5 %in% c("EBP")] <- 1
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE6 %in% c("EBP")] <- 1
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE7 %in% c("EBP")] <- 1
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE8 %in% c("EBP")] <- 1
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE9 %in% c("EBP")] <- 1

EntireCohort_pathvar$EBP[is.na(EntireCohort_pathvar$EBP)] <- 0
EntireCohort_pathvar$EBP <- as.numeric(EntireCohort_pathvar$EBP)
unique(EntireCohort_pathvar$EBP)
# 0 1

table(EntireCohort_pathvar$EBP)
#  0      1 
# 469556    243




### 'AMPD2' 
EntireCohort_pathvar$AMPD2 <- NA
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE1 %in% c("AMPD2")] <- 1
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE2 %in% c("AMPD2")] <- 1
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE3 %in% c("AMPD2")] <- 1
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE4 %in% c("AMPD2")] <- 1
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE5 %in% c("AMPD2")] <- 1
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE6 %in% c("AMPD2")] <- 1
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE7 %in% c("AMPD2")] <- 1
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE8 %in% c("AMPD2")] <- 1
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE9 %in% c("AMPD2")] <- 1

EntireCohort_pathvar$AMPD2[is.na(EntireCohort_pathvar$AMPD2)] <- 0
EntireCohort_pathvar$AMPD2 <- as.numeric(EntireCohort_pathvar$AMPD2)
unique(EntireCohort_pathvar$AMPD2)
# 0 1

table(EntireCohort_pathvar$AMPD2)
#  0      1 
# 468772   1027 






### 'LBR' 
EntireCohort_pathvar$LBR <- NA
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE1 %in% c("LBR")] <- 1
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE2 %in% c("LBR")] <- 1
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE3 %in% c("LBR")] <- 1
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE4 %in% c("LBR")] <- 1
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE5 %in% c("LBR")] <- 1
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE6 %in% c("LBR")] <- 1
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE7 %in% c("LBR")] <- 1
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE8 %in% c("LBR")] <- 1
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE9 %in% c("LBR")] <- 1

EntireCohort_pathvar$LBR[is.na(EntireCohort_pathvar$LBR)] <- 0
EntireCohort_pathvar$LBR <- as.numeric(EntireCohort_pathvar$LBR)
unique(EntireCohort_pathvar$LBR)
# 0 1

table(EntireCohort_pathvar$LBR)
#  0      1 
# 467126   2673





### 'AGPS' 
EntireCohort_pathvar$AGPS <- NA
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE1 %in% c("AGPS")] <- 1
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE2 %in% c("AGPS")] <- 1
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE3 %in% c("AGPS")] <- 1
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE4 %in% c("AGPS")] <- 1
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE5 %in% c("AGPS")] <- 1
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE6 %in% c("AGPS")] <- 1
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE7 %in% c("AGPS")] <- 1
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE8 %in% c("AGPS")] <- 1
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE9 %in% c("AGPS")] <- 1

EntireCohort_pathvar$AGPS[is.na(EntireCohort_pathvar$AGPS)] <- 0
EntireCohort_pathvar$AGPS <- as.numeric(EntireCohort_pathvar$AGPS)
unique(EntireCohort_pathvar$AGPS)
# 0 1

table(EntireCohort_pathvar$AGPS)
#  0      1 
# 469338    461 




### 'SLC35C1' 
EntireCohort_pathvar$SLC35C1 <- NA
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE1 %in% c("SLC35C1")] <- 1
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE2 %in% c("SLC35C1")] <- 1
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE3 %in% c("SLC35C1")] <- 1
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE4 %in% c("SLC35C1")] <- 1
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE5 %in% c("SLC35C1")] <- 1
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE6 %in% c("SLC35C1")] <- 1
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE7 %in% c("SLC35C1")] <- 1
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE8 %in% c("SLC35C1")] <- 1
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE9 %in% c("SLC35C1")] <- 1

EntireCohort_pathvar$SLC35C1[is.na(EntireCohort_pathvar$SLC35C1)] <- 0
EntireCohort_pathvar$SLC35C1 <- as.numeric(EntireCohort_pathvar$SLC35C1)
unique(EntireCohort_pathvar$SLC35C1)
# 0 1

table(EntireCohort_pathvar$SLC35C1)
#  0      1 
# 469230    569 




### 'PEX5'         
# no variants present




### 'SLC39A8'
EntireCohort_pathvar$SLC39A8 <- NA
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE1 %in% c("SLC39A8")] <- 1
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE2 %in% c("SLC39A8")] <- 1
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE3 %in% c("SLC39A8")] <- 1
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE4 %in% c("SLC39A8")] <- 1
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE5 %in% c("SLC39A8")] <- 1
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE6 %in% c("SLC39A8")] <- 1
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE7 %in% c("SLC39A8")] <- 1
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE8 %in% c("SLC39A8")] <- 1
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE9 %in% c("SLC39A8")] <- 1

EntireCohort_pathvar$SLC39A8[is.na(EntireCohort_pathvar$SLC39A8)] <- 0
EntireCohort_pathvar$SLC39A8 <- as.numeric(EntireCohort_pathvar$SLC39A8)
unique(EntireCohort_pathvar$SLC39A8)
# 0 1

table(EntireCohort_pathvar$SLC39A8)
#  0      1 
# 469387    412




### 'ASNS'  
EntireCohort_pathvar$ASNS <- NA
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE1 %in% c("ASNS")] <- 1
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE2 %in% c("ASNS")] <- 1
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE3 %in% c("ASNS")] <- 1
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE4 %in% c("ASNS")] <- 1
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE5 %in% c("ASNS")] <- 1
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE6 %in% c("ASNS")] <- 1
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE7 %in% c("ASNS")] <- 1
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE8 %in% c("ASNS")] <- 1
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE9 %in% c("ASNS")] <- 1

EntireCohort_pathvar$ASNS[is.na(EntireCohort_pathvar$ASNS)] <- 0
EntireCohort_pathvar$ASNS <- as.numeric(EntireCohort_pathvar$ASNS)
unique(EntireCohort_pathvar$ASNS)
# 0 1

table(EntireCohort_pathvar$ASNS)
#  0      1 
# 469052    747




### 'PSAT1'  
EntireCohort_pathvar$PSAT1 <- NA
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE1 %in% c("PSAT1")] <- 1
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE2 %in% c("PSAT1")] <- 1
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE3 %in% c("PSAT1")] <- 1
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE4 %in% c("PSAT1")] <- 1
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE5 %in% c("PSAT1")] <- 1
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE6 %in% c("PSAT1")] <- 1
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE7 %in% c("PSAT1")] <- 1
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE8 %in% c("PSAT1")] <- 1
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE9 %in% c("PSAT1")] <- 1

EntireCohort_pathvar$PSAT1[is.na(EntireCohort_pathvar$PSAT1)] <- 0
EntireCohort_pathvar$PSAT1 <- as.numeric(EntireCohort_pathvar$PSAT1)
unique(EntireCohort_pathvar$PSAT1)
# 0 1

table(EntireCohort_pathvar$PSAT1)
#   0      1 
# 468749   1050





### 'DPH1'
EntireCohort_pathvar$DPH1 <- NA
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE1 %in% c("DPH1")] <- 1
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE2 %in% c("DPH1")] <- 1
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE3 %in% c("DPH1")] <- 1
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE4 %in% c("DPH1")] <- 1
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE5 %in% c("DPH1")] <- 1
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE6 %in% c("DPH1")] <- 1
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE7 %in% c("DPH1")] <- 1
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE8 %in% c("DPH1")] <- 1
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE9 %in% c("DPH1")] <- 1

EntireCohort_pathvar$DPH1[is.na(EntireCohort_pathvar$DPH1)] <- 0
EntireCohort_pathvar$DPH1 <- as.numeric(EntireCohort_pathvar$DPH1)
unique(EntireCohort_pathvar$DPH1)
# 0 1

table(EntireCohort_pathvar$DPH1)
#  0      1 
# 463052   6747





### 'EIF2AK3'
EntireCohort_pathvar$EIF2AK3 <- NA
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE1 %in% c("EIF2AK3")] <- 1
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE2 %in% c("EIF2AK3")] <- 1
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE3 %in% c("EIF2AK3")] <- 1
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE4 %in% c("EIF2AK3")] <- 1
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE5 %in% c("EIF2AK3")] <- 1
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE6 %in% c("EIF2AK3")] <- 1
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE7 %in% c("EIF2AK3")] <- 1
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE8 %in% c("EIF2AK3")] <- 1
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE9 %in% c("EIF2AK3")] <- 1

EntireCohort_pathvar$EIF2AK3[is.na(EntireCohort_pathvar$EIF2AK3)] <- 0
EntireCohort_pathvar$EIF2AK3 <- as.numeric(EntireCohort_pathvar$EIF2AK3)
unique(EntireCohort_pathvar$EIF2AK3)
# 0 1

table(EntireCohort_pathvar$EIF2AK3)
#  0      1 
# 468939    860




### 'EIF2S3'  
EntireCohort_pathvar$EIF2S3 <- NA
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE1 %in% c("EIF2S3")] <- 1
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE2 %in% c("EIF2S3")] <- 1
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE3 %in% c("EIF2S3")] <- 1
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE4 %in% c("EIF2S3")] <- 1
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE5 %in% c("EIF2S3")] <- 1
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE6 %in% c("EIF2S3")] <- 1
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE7 %in% c("EIF2S3")] <- 1
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE8 %in% c("EIF2S3")] <- 1
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE9 %in% c("EIF2S3")] <- 1

EntireCohort_pathvar$EIF2S3[is.na(EntireCohort_pathvar$EIF2S3)] <- 0
EntireCohort_pathvar$EIF2S3 <- as.numeric(EntireCohort_pathvar$EIF2S3)
unique(EntireCohort_pathvar$EIF2S3)
# 0 1

table(EntireCohort_pathvar$EIF2S3)
#  0      1 
# 469784     15 




### 'EIF5A'         
EntireCohort_pathvar$EIF5A <- NA
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE1 %in% c("EIF5A")] <- 1
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE2 %in% c("EIF5A")] <- 1
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE3 %in% c("EIF5A")] <- 1
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE4 %in% c("EIF5A")] <- 1
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE5 %in% c("EIF5A")] <- 1
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE6 %in% c("EIF5A")] <- 1
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE7 %in% c("EIF5A")] <- 1
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE8 %in% c("EIF5A")] <- 1
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE9 %in% c("EIF5A")] <- 1

EntireCohort_pathvar$EIF5A[is.na(EntireCohort_pathvar$EIF5A)] <- 0
EntireCohort_pathvar$EIF5A <- as.numeric(EntireCohort_pathvar$EIF5A)
unique(EntireCohort_pathvar$EIF5A)
# 0 1

table(EntireCohort_pathvar$EIF5A)
#  0      1 
# 469760     39 




### 'RPL13'
EntireCohort_pathvar$RPL13 <- NA
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE1 %in% c("RPL13")] <- 1
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE2 %in% c("RPL13")] <- 1
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE3 %in% c("RPL13")] <- 1
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE4 %in% c("RPL13")] <- 1
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE5 %in% c("RPL13")] <- 1
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE6 %in% c("RPL13")] <- 1
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE7 %in% c("RPL13")] <- 1
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE8 %in% c("RPL13")] <- 1
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE9 %in% c("RPL13")] <- 1

EntireCohort_pathvar$RPL13[is.na(EntireCohort_pathvar$RPL13)] <- 0
EntireCohort_pathvar$RPL13 <- as.numeric(EntireCohort_pathvar$RPL13)
unique(EntireCohort_pathvar$RPL13)
# 0 1

table(EntireCohort_pathvar$RPL13)
#   0      1 
# 468905    894




### 'PHGDH'
EntireCohort_pathvar$PHGDH <- NA
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE1 %in% c("PHGDH")] <- 1
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE2 %in% c("PHGDH")] <- 1
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE3 %in% c("PHGDH")] <- 1
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE4 %in% c("PHGDH")] <- 1
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE5 %in% c("PHGDH")] <- 1
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE6 %in% c("PHGDH")] <- 1
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE7 %in% c("PHGDH")] <- 1
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE8 %in% c("PHGDH")] <- 1
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE9 %in% c("PHGDH")] <- 1

EntireCohort_pathvar$PHGDH[is.na(EntireCohort_pathvar$PHGDH)] <- 0
EntireCohort_pathvar$PHGDH <- as.numeric(EntireCohort_pathvar$PHGDH)
unique(EntireCohort_pathvar$PHGDH)
# 0 1

table(EntireCohort_pathvar$PHGDH)
#  0      1 
# 466675   3124




### 'NEPRO'   
EntireCohort_pathvar$NEPRO <- NA
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE1 %in% c("NEPRO")] <- 1
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE2 %in% c("NEPRO")] <- 1
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE3 %in% c("NEPRO")] <- 1
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE4 %in% c("NEPRO")] <- 1
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE5 %in% c("NEPRO")] <- 1
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE6 %in% c("NEPRO")] <- 1
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE7 %in% c("NEPRO")] <- 1
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE8 %in% c("NEPRO")] <- 1
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE9 %in% c("NEPRO")] <- 1

EntireCohort_pathvar$NEPRO[is.na(EntireCohort_pathvar$NEPRO)] <- 0
EntireCohort_pathvar$NEPRO <- as.numeric(EntireCohort_pathvar$NEPRO)
unique(EntireCohort_pathvar$NEPRO)
# 0 1

table(EntireCohort_pathvar$NEPRO)
#  0      1 
# 468580   1219






# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)






# BMgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

# go thr all the genes to check that they're present
grep("SERPINF1", EntireCohort_pathvar$GENE)
# all seem to be present except for the following:
# IFITM5 = 0
# COL1A2 = 0
# COL9A3 = 0
# MBTPS2 = 0


# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$BMgene <- NA
EntireCohort_pathvar$BMgene[EntireCohort_pathvar$GENE1 %in% c('FN1', 'ACAN', 'ADAMTS2', 'B3GALT6', 'B3GAT3', 'B4GALT7', 'CHST3', 'COL10A1', 'COL11A1', 'COL11A2', 'COL1A1', 'COL27A1', 'COL2A1', 'CRTAP', 'FAM20C', 'FKBP10', 'KDELR2', 'MATN3', 'MMP13', 'P3H1', 'SERPINH1', 'SLC10A7', 'SLC26A2', 'SLC35D1', 'CCN6', 'SERPINF1')] <- 1
EntireCohort_pathvar$BMgene[EntireCohort_pathvar$GENE2 %in% c('FN1', 'ACAN', 'ADAMTS2', 'B3GALT6', 'B3GAT3', 'B4GALT7', 'CHST3', 'COL10A1', 'COL11A1', 'COL11A2', 'COL1A1', 'COL27A1', 'COL2A1', 'CRTAP', 'FAM20C', 'FKBP10', 'KDELR2', 'MATN3', 'MMP13', 'P3H1', 'SERPINH1', 'SLC10A7', 'SLC26A2', 'SLC35D1', 'CCN6', 'SERPINF1')] <- 1
EntireCohort_pathvar$BMgene[EntireCohort_pathvar$GENE3 %in% c('FN1', 'ACAN', 'ADAMTS2', 'B3GALT6', 'B3GAT3', 'B4GALT7', 'CHST3', 'COL10A1', 'COL11A1', 'COL11A2', 'COL1A1', 'COL27A1', 'COL2A1', 'CRTAP', 'FAM20C', 'FKBP10', 'KDELR2', 'MATN3', 'MMP13', 'P3H1', 'SERPINH1', 'SLC10A7', 'SLC26A2', 'SLC35D1', 'CCN6', 'SERPINF1')] <- 1
EntireCohort_pathvar$BMgene[EntireCohort_pathvar$GENE4 %in% c('FN1', 'ACAN', 'ADAMTS2', 'B3GALT6', 'B3GAT3', 'B4GALT7', 'CHST3', 'COL10A1', 'COL11A1', 'COL11A2', 'COL1A1', 'COL27A1', 'COL2A1', 'CRTAP', 'FAM20C', 'FKBP10', 'KDELR2', 'MATN3', 'MMP13', 'P3H1', 'SERPINH1', 'SLC10A7', 'SLC26A2', 'SLC35D1', 'CCN6', 'SERPINF1')] <- 1
EntireCohort_pathvar$BMgene[EntireCohort_pathvar$GENE5 %in% c('FN1', 'ACAN', 'ADAMTS2', 'B3GALT6', 'B3GAT3', 'B4GALT7', 'CHST3', 'COL10A1', 'COL11A1', 'COL11A2', 'COL1A1', 'COL27A1', 'COL2A1', 'CRTAP', 'FAM20C', 'FKBP10', 'KDELR2', 'MATN3', 'MMP13', 'P3H1', 'SERPINH1', 'SLC10A7', 'SLC26A2', 'SLC35D1', 'CCN6', 'SERPINF1')] <- 1
EntireCohort_pathvar$BMgene[EntireCohort_pathvar$GENE6 %in% c('FN1', 'ACAN', 'ADAMTS2', 'B3GALT6', 'B3GAT3', 'B4GALT7', 'CHST3', 'COL10A1', 'COL11A1', 'COL11A2', 'COL1A1', 'COL27A1', 'COL2A1', 'CRTAP', 'FAM20C', 'FKBP10', 'KDELR2', 'MATN3', 'MMP13', 'P3H1', 'SERPINH1', 'SLC10A7', 'SLC26A2', 'SLC35D1', 'CCN6', 'SERPINF1')] <- 1
EntireCohort_pathvar$BMgene[EntireCohort_pathvar$GENE7 %in% c('FN1', 'ACAN', 'ADAMTS2', 'B3GALT6', 'B3GAT3', 'B4GALT7', 'CHST3', 'COL10A1', 'COL11A1', 'COL11A2', 'COL1A1', 'COL27A1', 'COL2A1', 'CRTAP', 'FAM20C', 'FKBP10', 'KDELR2', 'MATN3', 'MMP13', 'P3H1', 'SERPINH1', 'SLC10A7', 'SLC26A2', 'SLC35D1', 'CCN6', 'SERPINF1')] <- 1
EntireCohort_pathvar$BMgene[EntireCohort_pathvar$GENE8 %in% c('FN1', 'ACAN', 'ADAMTS2', 'B3GALT6', 'B3GAT3', 'B4GALT7', 'CHST3', 'COL10A1', 'COL11A1', 'COL11A2', 'COL1A1', 'COL27A1', 'COL2A1', 'CRTAP', 'FAM20C', 'FKBP10', 'KDELR2', 'MATN3', 'MMP13', 'P3H1', 'SERPINH1', 'SLC10A7', 'SLC26A2', 'SLC35D1', 'CCN6', 'SERPINF1')] <- 1
EntireCohort_pathvar$BMgene[EntireCohort_pathvar$GENE9 %in% c('FN1', 'ACAN', 'ADAMTS2', 'B3GALT6', 'B3GAT3', 'B4GALT7', 'CHST3', 'COL10A1', 'COL11A1', 'COL11A2', 'COL1A1', 'COL27A1', 'COL2A1', 'CRTAP', 'FAM20C', 'FKBP10', 'KDELR2', 'MATN3', 'MMP13', 'P3H1', 'SERPINH1', 'SLC10A7', 'SLC26A2', 'SLC35D1', 'CCN6', 'SERPINF1')] <- 1


EntireCohort_pathvar$BMgene[is.na(EntireCohort_pathvar$BMgene)] <- 0
EntireCohort_pathvar$BMgene <- as.numeric(EntireCohort_pathvar$BMgene)
unique(EntireCohort_pathvar$BMgene)
# 0 1

table(EntireCohort_pathvar$BMgene)
# 0      1 
# 385854  83945




# now to make a column for each of the BMgenes
### 'IFITM5'
# no variants present





### 'FN1', 
EntireCohort_pathvar$FN1 <- NA
EntireCohort_pathvar$FN1[EntireCohort_pathvar$GENE1 %in% c("FN1")] <- 1
EntireCohort_pathvar$FN1[EntireCohort_pathvar$GENE2 %in% c("FN1")] <- 1
EntireCohort_pathvar$FN1[EntireCohort_pathvar$GENE3 %in% c("FN1")] <- 1
EntireCohort_pathvar$FN1[EntireCohort_pathvar$GENE4 %in% c("FN1")] <- 1
EntireCohort_pathvar$FN1[EntireCohort_pathvar$GENE5 %in% c("FN1")] <- 1
EntireCohort_pathvar$FN1[EntireCohort_pathvar$GENE6 %in% c("FN1")] <- 1
EntireCohort_pathvar$FN1[EntireCohort_pathvar$GENE7 %in% c("FN1")] <- 1
EntireCohort_pathvar$FN1[EntireCohort_pathvar$GENE8 %in% c("FN1")] <- 1
EntireCohort_pathvar$FN1[EntireCohort_pathvar$GENE9 %in% c("FN1")] <- 1

EntireCohort_pathvar$FN1[is.na(EntireCohort_pathvar$FN1)] <- 0
EntireCohort_pathvar$FN1 <- as.numeric(EntireCohort_pathvar$FN1)
unique(EntireCohort_pathvar$FN1)
# 0 1

table(EntireCohort_pathvar$FN1)
# 0      1 
# 452633  17166




### 'ACAN', 
EntireCohort_pathvar$ACAN <- NA
EntireCohort_pathvar$ACAN[EntireCohort_pathvar$GENE1 %in% c("ACAN")] <- 1
EntireCohort_pathvar$ACAN[EntireCohort_pathvar$GENE2 %in% c("ACAN")] <- 1
EntireCohort_pathvar$ACAN[EntireCohort_pathvar$GENE3 %in% c("ACAN")] <- 1
EntireCohort_pathvar$ACAN[EntireCohort_pathvar$GENE4 %in% c("ACAN")] <- 1
EntireCohort_pathvar$ACAN[EntireCohort_pathvar$GENE5 %in% c("ACAN")] <- 1
EntireCohort_pathvar$ACAN[EntireCohort_pathvar$GENE6 %in% c("ACAN")] <- 1
EntireCohort_pathvar$ACAN[EntireCohort_pathvar$GENE7 %in% c("ACAN")] <- 1
EntireCohort_pathvar$ACAN[EntireCohort_pathvar$GENE8 %in% c("ACAN")] <- 1
EntireCohort_pathvar$ACAN[EntireCohort_pathvar$GENE9 %in% c("ACAN")] <- 1

EntireCohort_pathvar$ACAN[is.na(EntireCohort_pathvar$ACAN)] <- 0
EntireCohort_pathvar$ACAN <- as.numeric(EntireCohort_pathvar$ACAN)
unique(EntireCohort_pathvar$ACAN)
# 0 1

table(EntireCohort_pathvar$ACAN)
#  0      1 
# 467950   1849




### 'ADAMTS2', 
EntireCohort_pathvar$ADAMTS2 <- NA
EntireCohort_pathvar$ADAMTS2[EntireCohort_pathvar$GENE1 %in% c("ADAMTS2")] <- 1
EntireCohort_pathvar$ADAMTS2[EntireCohort_pathvar$GENE2 %in% c("ADAMTS2")] <- 1
EntireCohort_pathvar$ADAMTS2[EntireCohort_pathvar$GENE3 %in% c("ADAMTS2")] <- 1
EntireCohort_pathvar$ADAMTS2[EntireCohort_pathvar$GENE4 %in% c("ADAMTS2")] <- 1
EntireCohort_pathvar$ADAMTS2[EntireCohort_pathvar$GENE5 %in% c("ADAMTS2")] <- 1
EntireCohort_pathvar$ADAMTS2[EntireCohort_pathvar$GENE6 %in% c("ADAMTS2")] <- 1
EntireCohort_pathvar$ADAMTS2[EntireCohort_pathvar$GENE7 %in% c("ADAMTS2")] <- 1
EntireCohort_pathvar$ADAMTS2[EntireCohort_pathvar$GENE8 %in% c("ADAMTS2")] <- 1
EntireCohort_pathvar$ADAMTS2[EntireCohort_pathvar$GENE9 %in% c("ADAMTS2")] <- 1

EntireCohort_pathvar$ADAMTS2[is.na(EntireCohort_pathvar$ADAMTS2)] <- 0
EntireCohort_pathvar$ADAMTS2 <- as.numeric(EntireCohort_pathvar$ADAMTS2)
unique(EntireCohort_pathvar$ADAMTS2)
# 0 1

table(EntireCohort_pathvar$ADAMTS2)
#  0      1 
# 467635   2164




### 'B3GALT6', 
EntireCohort_pathvar$B3GALT6 <- NA
EntireCohort_pathvar$B3GALT6[EntireCohort_pathvar$GENE1 %in% c("B3GALT6")] <- 1
EntireCohort_pathvar$B3GALT6[EntireCohort_pathvar$GENE2 %in% c("B3GALT6")] <- 1
EntireCohort_pathvar$B3GALT6[EntireCohort_pathvar$GENE3 %in% c("B3GALT6")] <- 1
EntireCohort_pathvar$B3GALT6[EntireCohort_pathvar$GENE4 %in% c("B3GALT6")] <- 1
EntireCohort_pathvar$B3GALT6[EntireCohort_pathvar$GENE5 %in% c("B3GALT6")] <- 1
EntireCohort_pathvar$B3GALT6[EntireCohort_pathvar$GENE6 %in% c("B3GALT6")] <- 1
EntireCohort_pathvar$B3GALT6[EntireCohort_pathvar$GENE7 %in% c("B3GALT6")] <- 1
EntireCohort_pathvar$B3GALT6[EntireCohort_pathvar$GENE8 %in% c("B3GALT6")] <- 1
EntireCohort_pathvar$B3GALT6[EntireCohort_pathvar$GENE9 %in% c("B3GALT6")] <- 1

EntireCohort_pathvar$B3GALT6[is.na(EntireCohort_pathvar$B3GALT6)] <- 0
EntireCohort_pathvar$B3GALT6 <- as.numeric(EntireCohort_pathvar$B3GALT6)
unique(EntireCohort_pathvar$B3GALT6)
# 0 1

table(EntireCohort_pathvar$B3GALT6)
#   0      1 
# 469120    679




### 'B3GAT3', 
EntireCohort_pathvar$B3GAT3 <- NA
EntireCohort_pathvar$B3GAT3[EntireCohort_pathvar$GENE1 %in% c("B3GAT3")] <- 1
EntireCohort_pathvar$B3GAT3[EntireCohort_pathvar$GENE2 %in% c("B3GAT3")] <- 1
EntireCohort_pathvar$B3GAT3[EntireCohort_pathvar$GENE3 %in% c("B3GAT3")] <- 1
EntireCohort_pathvar$B3GAT3[EntireCohort_pathvar$GENE4 %in% c("B3GAT3")] <- 1
EntireCohort_pathvar$B3GAT3[EntireCohort_pathvar$GENE5 %in% c("B3GAT3")] <- 1
EntireCohort_pathvar$B3GAT3[EntireCohort_pathvar$GENE6 %in% c("B3GAT3")] <- 1
EntireCohort_pathvar$B3GAT3[EntireCohort_pathvar$GENE7 %in% c("B3GAT3")] <- 1
EntireCohort_pathvar$B3GAT3[EntireCohort_pathvar$GENE8 %in% c("B3GAT3")] <- 1
EntireCohort_pathvar$B3GAT3[EntireCohort_pathvar$GENE9 %in% c("B3GAT3")] <- 1

EntireCohort_pathvar$B3GAT3[is.na(EntireCohort_pathvar$B3GAT3)] <- 0
EntireCohort_pathvar$B3GAT3 <- as.numeric(EntireCohort_pathvar$B3GAT3)
unique(EntireCohort_pathvar$B3GAT3)
# 0 1

table(EntireCohort_pathvar$B3GAT3)
#  0      1 
# 469215    584




### 'B4GALT7', 
EntireCohort_pathvar$B4GALT7 <- NA
EntireCohort_pathvar$B4GALT7[EntireCohort_pathvar$GENE1 %in% c("B4GALT7")] <- 1
EntireCohort_pathvar$B4GALT7[EntireCohort_pathvar$GENE2 %in% c("B4GALT7")] <- 1
EntireCohort_pathvar$B4GALT7[EntireCohort_pathvar$GENE3 %in% c("B4GALT7")] <- 1
EntireCohort_pathvar$B4GALT7[EntireCohort_pathvar$GENE4 %in% c("B4GALT7")] <- 1
EntireCohort_pathvar$B4GALT7[EntireCohort_pathvar$GENE5 %in% c("B4GALT7")] <- 1
EntireCohort_pathvar$B4GALT7[EntireCohort_pathvar$GENE6 %in% c("B4GALT7")] <- 1
EntireCohort_pathvar$B4GALT7[EntireCohort_pathvar$GENE7 %in% c("B4GALT7")] <- 1
EntireCohort_pathvar$B4GALT7[EntireCohort_pathvar$GENE8 %in% c("B4GALT7")] <- 1
EntireCohort_pathvar$B4GALT7[EntireCohort_pathvar$GENE9 %in% c("B4GALT7")] <- 1

EntireCohort_pathvar$B4GALT7[is.na(EntireCohort_pathvar$B4GALT7)] <- 0
EntireCohort_pathvar$B4GALT7 <- as.numeric(EntireCohort_pathvar$B4GALT7)
unique(EntireCohort_pathvar$B4GALT7)
# 0 1

table(EntireCohort_pathvar$B4GALT7)
#   0      1 
# 467535   2264




### 'CHST3', 
EntireCohort_pathvar$CHST3 <- NA
EntireCohort_pathvar$CHST3[EntireCohort_pathvar$GENE1 %in% c("CHST3")] <- 1
EntireCohort_pathvar$CHST3[EntireCohort_pathvar$GENE2 %in% c("CHST3")] <- 1
EntireCohort_pathvar$CHST3[EntireCohort_pathvar$GENE3 %in% c("CHST3")] <- 1
EntireCohort_pathvar$CHST3[EntireCohort_pathvar$GENE4 %in% c("CHST3")] <- 1
EntireCohort_pathvar$CHST3[EntireCohort_pathvar$GENE5 %in% c("CHST3")] <- 1
EntireCohort_pathvar$CHST3[EntireCohort_pathvar$GENE6 %in% c("CHST3")] <- 1
EntireCohort_pathvar$CHST3[EntireCohort_pathvar$GENE7 %in% c("CHST3")] <- 1
EntireCohort_pathvar$CHST3[EntireCohort_pathvar$GENE8 %in% c("CHST3")] <- 1
EntireCohort_pathvar$CHST3[EntireCohort_pathvar$GENE9 %in% c("CHST3")] <- 1

EntireCohort_pathvar$CHST3[is.na(EntireCohort_pathvar$CHST3)] <- 0
EntireCohort_pathvar$CHST3 <- as.numeric(EntireCohort_pathvar$CHST3)
unique(EntireCohort_pathvar$CHST3)
# 0 1

table(EntireCohort_pathvar$CHST3)
#  0      1 
# 469150    649




### 'COL10A1', 
EntireCohort_pathvar$COL10A1 <- NA
EntireCohort_pathvar$COL10A1[EntireCohort_pathvar$GENE1 %in% c("COL10A1")] <- 1
EntireCohort_pathvar$COL10A1[EntireCohort_pathvar$GENE2 %in% c("COL10A1")] <- 1
EntireCohort_pathvar$COL10A1[EntireCohort_pathvar$GENE3 %in% c("COL10A1")] <- 1
EntireCohort_pathvar$COL10A1[EntireCohort_pathvar$GENE4 %in% c("COL10A1")] <- 1
EntireCohort_pathvar$COL10A1[EntireCohort_pathvar$GENE5 %in% c("COL10A1")] <- 1
EntireCohort_pathvar$COL10A1[EntireCohort_pathvar$GENE6 %in% c("COL10A1")] <- 1
EntireCohort_pathvar$COL10A1[EntireCohort_pathvar$GENE7 %in% c("COL10A1")] <- 1
EntireCohort_pathvar$COL10A1[EntireCohort_pathvar$GENE8 %in% c("COL10A1")] <- 1
EntireCohort_pathvar$COL10A1[EntireCohort_pathvar$GENE9 %in% c("COL10A1")] <- 1

EntireCohort_pathvar$COL10A1[is.na(EntireCohort_pathvar$COL10A1)] <- 0
EntireCohort_pathvar$COL10A1 <- as.numeric(EntireCohort_pathvar$COL10A1)
unique(EntireCohort_pathvar$COL10A1)
# 0 1

table(EntireCohort_pathvar$COL10A1)
#    0      1 
# 467697   2102





### 'COL11A1', 
EntireCohort_pathvar$COL11A1 <- NA
EntireCohort_pathvar$COL11A1[EntireCohort_pathvar$GENE1 %in% c("COL11A1")] <- 1
EntireCohort_pathvar$COL11A1[EntireCohort_pathvar$GENE2 %in% c("COL11A1")] <- 1
EntireCohort_pathvar$COL11A1[EntireCohort_pathvar$GENE3 %in% c("COL11A1")] <- 1
EntireCohort_pathvar$COL11A1[EntireCohort_pathvar$GENE4 %in% c("COL11A1")] <- 1
EntireCohort_pathvar$COL11A1[EntireCohort_pathvar$GENE5 %in% c("COL11A1")] <- 1
EntireCohort_pathvar$COL11A1[EntireCohort_pathvar$GENE6 %in% c("COL11A1")] <- 1
EntireCohort_pathvar$COL11A1[EntireCohort_pathvar$GENE7 %in% c("COL11A1")] <- 1
EntireCohort_pathvar$COL11A1[EntireCohort_pathvar$GENE8 %in% c("COL11A1")] <- 1
EntireCohort_pathvar$COL11A1[EntireCohort_pathvar$GENE9 %in% c("COL11A1")] <- 1

EntireCohort_pathvar$COL11A1[is.na(EntireCohort_pathvar$COL11A1)] <- 0
EntireCohort_pathvar$COL11A1 <- as.numeric(EntireCohort_pathvar$COL11A1)
unique(EntireCohort_pathvar$COL11A1)
# 0 1

table(EntireCohort_pathvar$COL11A1)
#  0      1 
# 461913   7886 




### 'COL11A2', 
EntireCohort_pathvar$COL11A2 <- NA
EntireCohort_pathvar$COL11A2[EntireCohort_pathvar$GENE1 %in% c("COL11A2")] <- 1
EntireCohort_pathvar$COL11A2[EntireCohort_pathvar$GENE2 %in% c("COL11A2")] <- 1
EntireCohort_pathvar$COL11A2[EntireCohort_pathvar$GENE3 %in% c("COL11A2")] <- 1
EntireCohort_pathvar$COL11A2[EntireCohort_pathvar$GENE4 %in% c("COL11A2")] <- 1
EntireCohort_pathvar$COL11A2[EntireCohort_pathvar$GENE5 %in% c("COL11A2")] <- 1
EntireCohort_pathvar$COL11A2[EntireCohort_pathvar$GENE6 %in% c("COL11A2")] <- 1
EntireCohort_pathvar$COL11A2[EntireCohort_pathvar$GENE7 %in% c("COL11A2")] <- 1
EntireCohort_pathvar$COL11A2[EntireCohort_pathvar$GENE8 %in% c("COL11A2")] <- 1
EntireCohort_pathvar$COL11A2[EntireCohort_pathvar$GENE9 %in% c("COL11A2")] <- 1

EntireCohort_pathvar$COL11A2[is.na(EntireCohort_pathvar$COL11A2)] <- 0
EntireCohort_pathvar$COL11A2 <- as.numeric(EntireCohort_pathvar$COL11A2)
unique(EntireCohort_pathvar$COL11A2)
# 0 1

table(EntireCohort_pathvar$COL11A2)
#  0      1 
# 462483   7316




### 'COL1A1', 
EntireCohort_pathvar$COL1A1 <- NA
EntireCohort_pathvar$COL1A1[EntireCohort_pathvar$GENE1 %in% c("COL1A1")] <- 1
EntireCohort_pathvar$COL1A1[EntireCohort_pathvar$GENE2 %in% c("COL1A1")] <- 1
EntireCohort_pathvar$COL1A1[EntireCohort_pathvar$GENE3 %in% c("COL1A1")] <- 1
EntireCohort_pathvar$COL1A1[EntireCohort_pathvar$GENE4 %in% c("COL1A1")] <- 1
EntireCohort_pathvar$COL1A1[EntireCohort_pathvar$GENE5 %in% c("COL1A1")] <- 1
EntireCohort_pathvar$COL1A1[EntireCohort_pathvar$GENE6 %in% c("COL1A1")] <- 1
EntireCohort_pathvar$COL1A1[EntireCohort_pathvar$GENE7 %in% c("COL1A1")] <- 1
EntireCohort_pathvar$COL1A1[EntireCohort_pathvar$GENE8 %in% c("COL1A1")] <- 1
EntireCohort_pathvar$COL1A1[EntireCohort_pathvar$GENE9 %in% c("COL1A1")] <- 1

EntireCohort_pathvar$COL1A1[is.na(EntireCohort_pathvar$COL1A1)] <- 0
EntireCohort_pathvar$COL1A1 <- as.numeric(EntireCohort_pathvar$COL1A1)
unique(EntireCohort_pathvar$COL1A1)
# 0 1

table(EntireCohort_pathvar$COL1A1)
#  0      1 
# 466541   3258 




### 'COL1A2', 
# no variants present




### 'COL27A1', 
EntireCohort_pathvar$COL27A1 <- NA
EntireCohort_pathvar$COL27A1[EntireCohort_pathvar$GENE1 %in% c("COL27A1")] <- 1
EntireCohort_pathvar$COL27A1[EntireCohort_pathvar$GENE2 %in% c("COL27A1")] <- 1
EntireCohort_pathvar$COL27A1[EntireCohort_pathvar$GENE3 %in% c("COL27A1")] <- 1
EntireCohort_pathvar$COL27A1[EntireCohort_pathvar$GENE4 %in% c("COL27A1")] <- 1
EntireCohort_pathvar$COL27A1[EntireCohort_pathvar$GENE5 %in% c("COL27A1")] <- 1
EntireCohort_pathvar$COL27A1[EntireCohort_pathvar$GENE6 %in% c("COL27A1")] <- 1
EntireCohort_pathvar$COL27A1[EntireCohort_pathvar$GENE7 %in% c("COL27A1")] <- 1
EntireCohort_pathvar$COL27A1[EntireCohort_pathvar$GENE8 %in% c("COL27A1")] <- 1
EntireCohort_pathvar$COL27A1[EntireCohort_pathvar$GENE9 %in% c("COL27A1")] <- 1

EntireCohort_pathvar$COL27A1[is.na(EntireCohort_pathvar$COL27A1)] <- 0
EntireCohort_pathvar$COL27A1 <- as.numeric(EntireCohort_pathvar$COL27A1)
unique(EntireCohort_pathvar$COL27A1)
# 0 1

table(EntireCohort_pathvar$COL27A1)
#   0      1 
# 463895   5904




### 'COL2A1', 
EntireCohort_pathvar$COL2A1 <- NA
EntireCohort_pathvar$COL2A1[EntireCohort_pathvar$GENE1 %in% c("COL2A1")] <- 1
EntireCohort_pathvar$COL2A1[EntireCohort_pathvar$GENE2 %in% c("COL2A1")] <- 1
EntireCohort_pathvar$COL2A1[EntireCohort_pathvar$GENE3 %in% c("COL2A1")] <- 1
EntireCohort_pathvar$COL2A1[EntireCohort_pathvar$GENE4 %in% c("COL2A1")] <- 1
EntireCohort_pathvar$COL2A1[EntireCohort_pathvar$GENE5 %in% c("COL2A1")] <- 1
EntireCohort_pathvar$COL2A1[EntireCohort_pathvar$GENE6 %in% c("COL2A1")] <- 1
EntireCohort_pathvar$COL2A1[EntireCohort_pathvar$GENE7 %in% c("COL2A1")] <- 1
EntireCohort_pathvar$COL2A1[EntireCohort_pathvar$GENE8 %in% c("COL2A1")] <- 1
EntireCohort_pathvar$COL2A1[EntireCohort_pathvar$GENE9 %in% c("COL2A1")] <- 1

EntireCohort_pathvar$COL2A1[is.na(EntireCohort_pathvar$COL2A1)] <- 0
EntireCohort_pathvar$COL2A1 <- as.numeric(EntireCohort_pathvar$COL2A1)
unique(EntireCohort_pathvar$COL2A1)
# 0 1

table(EntireCohort_pathvar$COL2A1)
#  0      1 
# 467491   2308





### 'COL9A3', 
# no varians present





### 'CRTAP', 
EntireCohort_pathvar$CRTAP <- NA
EntireCohort_pathvar$CRTAP[EntireCohort_pathvar$GENE1 %in% c("CRTAP")] <- 1
EntireCohort_pathvar$CRTAP[EntireCohort_pathvar$GENE2 %in% c("CRTAP")] <- 1
EntireCohort_pathvar$CRTAP[EntireCohort_pathvar$GENE3 %in% c("CRTAP")] <- 1
EntireCohort_pathvar$CRTAP[EntireCohort_pathvar$GENE4 %in% c("CRTAP")] <- 1
EntireCohort_pathvar$CRTAP[EntireCohort_pathvar$GENE5 %in% c("CRTAP")] <- 1
EntireCohort_pathvar$CRTAP[EntireCohort_pathvar$GENE6 %in% c("CRTAP")] <- 1
EntireCohort_pathvar$CRTAP[EntireCohort_pathvar$GENE7 %in% c("CRTAP")] <- 1
EntireCohort_pathvar$CRTAP[EntireCohort_pathvar$GENE8 %in% c("CRTAP")] <- 1
EntireCohort_pathvar$CRTAP[EntireCohort_pathvar$GENE9 %in% c("CRTAP")] <- 1

EntireCohort_pathvar$CRTAP[is.na(EntireCohort_pathvar$CRTAP)] <- 0
EntireCohort_pathvar$CRTAP <- as.numeric(EntireCohort_pathvar$CRTAP)
unique(EntireCohort_pathvar$CRTAP)
# 0 1

table(EntireCohort_pathvar$CRTAP)
#  0      1 
# 467451   2348




### 'FAM20C', 
EntireCohort_pathvar$FAM20C <- NA
EntireCohort_pathvar$FAM20C[EntireCohort_pathvar$GENE1 %in% c("FAM20C")] <- 1
EntireCohort_pathvar$FAM20C[EntireCohort_pathvar$GENE2 %in% c("FAM20C")] <- 1
EntireCohort_pathvar$FAM20C[EntireCohort_pathvar$GENE3 %in% c("FAM20C")] <- 1
EntireCohort_pathvar$FAM20C[EntireCohort_pathvar$GENE4 %in% c("FAM20C")] <- 1
EntireCohort_pathvar$FAM20C[EntireCohort_pathvar$GENE5 %in% c("FAM20C")] <- 1
EntireCohort_pathvar$FAM20C[EntireCohort_pathvar$GENE6 %in% c("FAM20C")] <- 1
EntireCohort_pathvar$FAM20C[EntireCohort_pathvar$GENE7 %in% c("FAM20C")] <- 1
EntireCohort_pathvar$FAM20C[EntireCohort_pathvar$GENE8 %in% c("FAM20C")] <- 1
EntireCohort_pathvar$FAM20C[EntireCohort_pathvar$GENE9 %in% c("FAM20C")] <- 1

EntireCohort_pathvar$FAM20C[is.na(EntireCohort_pathvar$FAM20C)] <- 0
EntireCohort_pathvar$FAM20C <- as.numeric(EntireCohort_pathvar$FAM20C)
unique(EntireCohort_pathvar$FAM20C)
# 0 1

table(EntireCohort_pathvar$FAM20C)
# 0      1 
# 468728   1071 




### 'FKBP10', 
EntireCohort_pathvar$FKBP10 <- NA
EntireCohort_pathvar$FKBP10[EntireCohort_pathvar$GENE1 %in% c("FKBP10")] <- 1
EntireCohort_pathvar$FKBP10[EntireCohort_pathvar$GENE2 %in% c("FKBP10")] <- 1
EntireCohort_pathvar$FKBP10[EntireCohort_pathvar$GENE3 %in% c("FKBP10")] <- 1
EntireCohort_pathvar$FKBP10[EntireCohort_pathvar$GENE4 %in% c("FKBP10")] <- 1
EntireCohort_pathvar$FKBP10[EntireCohort_pathvar$GENE5 %in% c("FKBP10")] <- 1
EntireCohort_pathvar$FKBP10[EntireCohort_pathvar$GENE6 %in% c("FKBP10")] <- 1
EntireCohort_pathvar$FKBP10[EntireCohort_pathvar$GENE7 %in% c("FKBP10")] <- 1
EntireCohort_pathvar$FKBP10[EntireCohort_pathvar$GENE8 %in% c("FKBP10")] <- 1
EntireCohort_pathvar$FKBP10[EntireCohort_pathvar$GENE9 %in% c("FKBP10")] <- 1

EntireCohort_pathvar$FKBP10[is.na(EntireCohort_pathvar$FKBP10)] <- 0
EntireCohort_pathvar$FKBP10 <- as.numeric(EntireCohort_pathvar$FKBP10)
unique(EntireCohort_pathvar$FKBP10)
# 0 1

table(EntireCohort_pathvar$FKBP10)
#  0      1 
# 463020   6779




### 'KDELR2', 
EntireCohort_pathvar$KDELR2 <- NA
EntireCohort_pathvar$KDELR2[EntireCohort_pathvar$GENE1 %in% c("KDELR2")] <- 1
EntireCohort_pathvar$KDELR2[EntireCohort_pathvar$GENE2 %in% c("KDELR2")] <- 1
EntireCohort_pathvar$KDELR2[EntireCohort_pathvar$GENE3 %in% c("KDELR2")] <- 1
EntireCohort_pathvar$KDELR2[EntireCohort_pathvar$GENE4 %in% c("KDELR2")] <- 1
EntireCohort_pathvar$KDELR2[EntireCohort_pathvar$GENE5 %in% c("KDELR2")] <- 1
EntireCohort_pathvar$KDELR2[EntireCohort_pathvar$GENE6 %in% c("KDELR2")] <- 1
EntireCohort_pathvar$KDELR2[EntireCohort_pathvar$GENE7 %in% c("KDELR2")] <- 1
EntireCohort_pathvar$KDELR2[EntireCohort_pathvar$GENE8 %in% c("KDELR2")] <- 1
EntireCohort_pathvar$KDELR2[EntireCohort_pathvar$GENE9 %in% c("KDELR2")] <- 1

EntireCohort_pathvar$KDELR2[is.na(EntireCohort_pathvar$KDELR2)] <- 0
EntireCohort_pathvar$KDELR2 <- as.numeric(EntireCohort_pathvar$KDELR2)
unique(EntireCohort_pathvar$KDELR2)
# 0 1

table(EntireCohort_pathvar$KDELR2)
#  0      1 
# 469308    491 




### 'MATN3', 
EntireCohort_pathvar$MATN3 <- NA
EntireCohort_pathvar$MATN3[EntireCohort_pathvar$GENE1 %in% c("MATN3")] <- 1
EntireCohort_pathvar$MATN3[EntireCohort_pathvar$GENE2 %in% c("MATN3")] <- 1
EntireCohort_pathvar$MATN3[EntireCohort_pathvar$GENE3 %in% c("MATN3")] <- 1
EntireCohort_pathvar$MATN3[EntireCohort_pathvar$GENE4 %in% c("MATN3")] <- 1
EntireCohort_pathvar$MATN3[EntireCohort_pathvar$GENE5 %in% c("MATN3")] <- 1
EntireCohort_pathvar$MATN3[EntireCohort_pathvar$GENE6 %in% c("MATN3")] <- 1
EntireCohort_pathvar$MATN3[EntireCohort_pathvar$GENE7 %in% c("MATN3")] <- 1
EntireCohort_pathvar$MATN3[EntireCohort_pathvar$GENE8 %in% c("MATN3")] <- 1
EntireCohort_pathvar$MATN3[EntireCohort_pathvar$GENE9 %in% c("MATN3")] <- 1

EntireCohort_pathvar$MATN3[is.na(EntireCohort_pathvar$MATN3)] <- 0
EntireCohort_pathvar$MATN3 <- as.numeric(EntireCohort_pathvar$MATN3)
unique(EntireCohort_pathvar$MATN3)
# 0 1

table(EntireCohort_pathvar$MATN3)
#  0      1 
# 453775  16024




### 'MBTPS2', 
# no variants present






### 'MMP13', 
EntireCohort_pathvar$MMP13 <- NA
EntireCohort_pathvar$MMP13[EntireCohort_pathvar$GENE1 %in% c("MMP13")] <- 1
EntireCohort_pathvar$MMP13[EntireCohort_pathvar$GENE2 %in% c("MMP13")] <- 1
EntireCohort_pathvar$MMP13[EntireCohort_pathvar$GENE3 %in% c("MMP13")] <- 1
EntireCohort_pathvar$MMP13[EntireCohort_pathvar$GENE4 %in% c("MMP13")] <- 1
EntireCohort_pathvar$MMP13[EntireCohort_pathvar$GENE5 %in% c("MMP13")] <- 1
EntireCohort_pathvar$MMP13[EntireCohort_pathvar$GENE6 %in% c("MMP13")] <- 1
EntireCohort_pathvar$MMP13[EntireCohort_pathvar$GENE7 %in% c("MMP13")] <- 1
EntireCohort_pathvar$MMP13[EntireCohort_pathvar$GENE8 %in% c("MMP13")] <- 1
EntireCohort_pathvar$MMP13[EntireCohort_pathvar$GENE9 %in% c("MMP13")] <- 1

EntireCohort_pathvar$MMP13[is.na(EntireCohort_pathvar$MMP13)] <- 0
EntireCohort_pathvar$MMP13 <- as.numeric(EntireCohort_pathvar$MMP13)
unique(EntireCohort_pathvar$MMP13)
# 0 1

table(EntireCohort_pathvar$MMP13)
#  0      1 
# 468932    867





### 'P3H1', 
EntireCohort_pathvar$P3H1 <- NA
EntireCohort_pathvar$P3H1[EntireCohort_pathvar$GENE1 %in% c("P3H1")] <- 1
EntireCohort_pathvar$P3H1[EntireCohort_pathvar$GENE2 %in% c("P3H1")] <- 1
EntireCohort_pathvar$P3H1[EntireCohort_pathvar$GENE3 %in% c("P3H1")] <- 1
EntireCohort_pathvar$P3H1[EntireCohort_pathvar$GENE4 %in% c("P3H1")] <- 1
EntireCohort_pathvar$P3H1[EntireCohort_pathvar$GENE5 %in% c("P3H1")] <- 1
EntireCohort_pathvar$P3H1[EntireCohort_pathvar$GENE6 %in% c("P3H1")] <- 1
EntireCohort_pathvar$P3H1[EntireCohort_pathvar$GENE7 %in% c("P3H1")] <- 1
EntireCohort_pathvar$P3H1[EntireCohort_pathvar$GENE8 %in% c("P3H1")] <- 1
EntireCohort_pathvar$P3H1[EntireCohort_pathvar$GENE9 %in% c("P3H1")] <- 1

EntireCohort_pathvar$P3H1[is.na(EntireCohort_pathvar$P3H1)] <- 0
EntireCohort_pathvar$P3H1 <- as.numeric(EntireCohort_pathvar$P3H1)
unique(EntireCohort_pathvar$P3H1)
# 0 1

table(EntireCohort_pathvar$P3H1)
#  0      1 
# 466683   3116




### 'SERPINH1', 
EntireCohort_pathvar$SERPINH1 <- NA
EntireCohort_pathvar$SERPINH1[EntireCohort_pathvar$GENE1 %in% c("SERPINH1")] <- 1
EntireCohort_pathvar$SERPINH1[EntireCohort_pathvar$GENE2 %in% c("SERPINH1")] <- 1
EntireCohort_pathvar$SERPINH1[EntireCohort_pathvar$GENE3 %in% c("SERPINH1")] <- 1
EntireCohort_pathvar$SERPINH1[EntireCohort_pathvar$GENE4 %in% c("SERPINH1")] <- 1
EntireCohort_pathvar$SERPINH1[EntireCohort_pathvar$GENE5 %in% c("SERPINH1")] <- 1
EntireCohort_pathvar$SERPINH1[EntireCohort_pathvar$GENE6 %in% c("SERPINH1")] <- 1
EntireCohort_pathvar$SERPINH1[EntireCohort_pathvar$GENE7 %in% c("SERPINH1")] <- 1
EntireCohort_pathvar$SERPINH1[EntireCohort_pathvar$GENE8 %in% c("SERPINH1")] <- 1
EntireCohort_pathvar$SERPINH1[EntireCohort_pathvar$GENE9 %in% c("SERPINH1")] <- 1

EntireCohort_pathvar$SERPINH1[is.na(EntireCohort_pathvar$SERPINH1)] <- 0
EntireCohort_pathvar$SERPINH1 <- as.numeric(EntireCohort_pathvar$SERPINH1)
unique(EntireCohort_pathvar$SERPINH1)
# 0 1

table(EntireCohort_pathvar$SERPINH1)
#  0      1 
# 469092    707




### 'SLC10A7', 
EntireCohort_pathvar$SLC10A7 <- NA
EntireCohort_pathvar$SLC10A7[EntireCohort_pathvar$GENE1 %in% c("SLC10A7")] <- 1
EntireCohort_pathvar$SLC10A7[EntireCohort_pathvar$GENE2 %in% c("SLC10A7")] <- 1
EntireCohort_pathvar$SLC10A7[EntireCohort_pathvar$GENE3 %in% c("SLC10A7")] <- 1
EntireCohort_pathvar$SLC10A7[EntireCohort_pathvar$GENE4 %in% c("SLC10A7")] <- 1
EntireCohort_pathvar$SLC10A7[EntireCohort_pathvar$GENE5 %in% c("SLC10A7")] <- 1
EntireCohort_pathvar$SLC10A7[EntireCohort_pathvar$GENE6 %in% c("SLC10A7")] <- 1
EntireCohort_pathvar$SLC10A7[EntireCohort_pathvar$GENE7 %in% c("SLC10A7")] <- 1
EntireCohort_pathvar$SLC10A7[EntireCohort_pathvar$GENE8 %in% c("SLC10A7")] <- 1
EntireCohort_pathvar$SLC10A7[EntireCohort_pathvar$GENE9 %in% c("SLC10A7")] <- 1

EntireCohort_pathvar$SLC10A7[is.na(EntireCohort_pathvar$SLC10A7)] <- 0
EntireCohort_pathvar$SLC10A7 <- as.numeric(EntireCohort_pathvar$SLC10A7)
unique(EntireCohort_pathvar$SLC10A7)
# 0 1

table(EntireCohort_pathvar$SLC10A7)
#  0      1 
# 469330    469 




### 'SLC26A2', 
EntireCohort_pathvar$SLC26A2 <- NA
EntireCohort_pathvar$SLC26A2[EntireCohort_pathvar$GENE1 %in% c("SLC26A2")] <- 1
EntireCohort_pathvar$SLC26A2[EntireCohort_pathvar$GENE2 %in% c("SLC26A2")] <- 1
EntireCohort_pathvar$SLC26A2[EntireCohort_pathvar$GENE3 %in% c("SLC26A2")] <- 1
EntireCohort_pathvar$SLC26A2[EntireCohort_pathvar$GENE4 %in% c("SLC26A2")] <- 1
EntireCohort_pathvar$SLC26A2[EntireCohort_pathvar$GENE5 %in% c("SLC26A2")] <- 1
EntireCohort_pathvar$SLC26A2[EntireCohort_pathvar$GENE6 %in% c("SLC26A2")] <- 1
EntireCohort_pathvar$SLC26A2[EntireCohort_pathvar$GENE7 %in% c("SLC26A2")] <- 1
EntireCohort_pathvar$SLC26A2[EntireCohort_pathvar$GENE8 %in% c("SLC26A2")] <- 1
EntireCohort_pathvar$SLC26A2[EntireCohort_pathvar$GENE9 %in% c("SLC26A2")] <- 1

EntireCohort_pathvar$SLC26A2[is.na(EntireCohort_pathvar$SLC26A2)] <- 0
EntireCohort_pathvar$SLC26A2 <- as.numeric(EntireCohort_pathvar$SLC26A2)
unique(EntireCohort_pathvar$SLC26A2)
# 0 1

table(EntireCohort_pathvar$SLC26A2)
#  0      1 
# 466746   3053 




### 'SLC35D1', 
EntireCohort_pathvar$SLC35D1 <- NA
EntireCohort_pathvar$SLC35D1[EntireCohort_pathvar$GENE1 %in% c("SLC35D1")] <- 1
EntireCohort_pathvar$SLC35D1[EntireCohort_pathvar$GENE2 %in% c("SLC35D1")] <- 1
EntireCohort_pathvar$SLC35D1[EntireCohort_pathvar$GENE3 %in% c("SLC35D1")] <- 1
EntireCohort_pathvar$SLC35D1[EntireCohort_pathvar$GENE4 %in% c("SLC35D1")] <- 1
EntireCohort_pathvar$SLC35D1[EntireCohort_pathvar$GENE5 %in% c("SLC35D1")] <- 1
EntireCohort_pathvar$SLC35D1[EntireCohort_pathvar$GENE6 %in% c("SLC35D1")] <- 1
EntireCohort_pathvar$SLC35D1[EntireCohort_pathvar$GENE7 %in% c("SLC35D1")] <- 1
EntireCohort_pathvar$SLC35D1[EntireCohort_pathvar$GENE8 %in% c("SLC35D1")] <- 1
EntireCohort_pathvar$SLC35D1[EntireCohort_pathvar$GENE9 %in% c("SLC35D1")] <- 1

EntireCohort_pathvar$SLC35D1[is.na(EntireCohort_pathvar$SLC35D1)] <- 0
EntireCohort_pathvar$SLC35D1 <- as.numeric(EntireCohort_pathvar$SLC35D1)
unique(EntireCohort_pathvar$SLC35D1)
# 0 1

table(EntireCohort_pathvar$SLC35D1)
#  0      1 
# 469327    472




### 'CCN6', 
EntireCohort_pathvar$CCN6 <- NA
EntireCohort_pathvar$CCN6[EntireCohort_pathvar$GENE1 %in% c("CCN6")] <- 1
EntireCohort_pathvar$CCN6[EntireCohort_pathvar$GENE2 %in% c("CCN6")] <- 1
EntireCohort_pathvar$CCN6[EntireCohort_pathvar$GENE3 %in% c("CCN6")] <- 1
EntireCohort_pathvar$CCN6[EntireCohort_pathvar$GENE4 %in% c("CCN6")] <- 1
EntireCohort_pathvar$CCN6[EntireCohort_pathvar$GENE5 %in% c("CCN6")] <- 1
EntireCohort_pathvar$CCN6[EntireCohort_pathvar$GENE6 %in% c("CCN6")] <- 1
EntireCohort_pathvar$CCN6[EntireCohort_pathvar$GENE7 %in% c("CCN6")] <- 1
EntireCohort_pathvar$CCN6[EntireCohort_pathvar$GENE8 %in% c("CCN6")] <- 1
EntireCohort_pathvar$CCN6[EntireCohort_pathvar$GENE9 %in% c("CCN6")] <- 1

EntireCohort_pathvar$CCN6[is.na(EntireCohort_pathvar$CCN6)] <- 0
EntireCohort_pathvar$CCN6 <- as.numeric(EntireCohort_pathvar$CCN6)
unique(EntireCohort_pathvar$CCN6)
# 0 1

table(EntireCohort_pathvar$CCN6)
#  0      1 
# 468805    994




### 'SERPINF1'
EntireCohort_pathvar$SERPINF1 <- NA
EntireCohort_pathvar$SERPINF1[EntireCohort_pathvar$GENE1 %in% c("SERPINF1")] <- 1
EntireCohort_pathvar$SERPINF1[EntireCohort_pathvar$GENE2 %in% c("SERPINF1")] <- 1
EntireCohort_pathvar$SERPINF1[EntireCohort_pathvar$GENE3 %in% c("SERPINF1")] <- 1
EntireCohort_pathvar$SERPINF1[EntireCohort_pathvar$GENE4 %in% c("SERPINF1")] <- 1
EntireCohort_pathvar$SERPINF1[EntireCohort_pathvar$GENE5 %in% c("SERPINF1")] <- 1
EntireCohort_pathvar$SERPINF1[EntireCohort_pathvar$GENE6 %in% c("SERPINF1")] <- 1
EntireCohort_pathvar$SERPINF1[EntireCohort_pathvar$GENE7 %in% c("SERPINF1")] <- 1
EntireCohort_pathvar$SERPINF1[EntireCohort_pathvar$GENE8 %in% c("SERPINF1")] <- 1
EntireCohort_pathvar$SERPINF1[EntireCohort_pathvar$GENE9 %in% c("SERPINF1")] <- 1

EntireCohort_pathvar$SERPINF1[is.na(EntireCohort_pathvar$SERPINF1)] <- 0
EntireCohort_pathvar$SERPINF1 <- as.numeric(EntireCohort_pathvar$SERPINF1)
unique(EntireCohort_pathvar$SERPINF1)
# 0 1

table(EntireCohort_pathvar$SERPINF1)
#  0      1 
# 468670   1129 







# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)






# Bonegenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

# go thr all the genes to check that they're present
grep("PPP3CA", EntireCohort_pathvar$GENE)
# all seem to be present, except for:
# PIEZO2 = 0
# NOTCH2 = 0
# PPP3CA = 0





# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$Bonegene <- NA
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE1 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2')] <- 1
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE2 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2')] <- 1
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE3 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2')] <- 1
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE4 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2')] <- 1
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE5 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2')] <- 1
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE6 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2')] <- 1
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE7 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2')] <- 1
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE8 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2')] <- 1
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE9 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2')] <- 1


EntireCohort_pathvar$Bonegene[is.na(EntireCohort_pathvar$Bonegene)] <- 0
EntireCohort_pathvar$Bonegene <- as.numeric(EntireCohort_pathvar$Bonegene)
unique(EntireCohort_pathvar$Bonegene)
# 0 1

table(EntireCohort_pathvar$Bonegene)
#  0      1 
# 443865  25934 




# now to make a column for each of the Bonegenes
### 'DYNC2H1'
EntireCohort_pathvar$DYNC2H1 <- NA
EntireCohort_pathvar$DYNC2H1[EntireCohort_pathvar$GENE1 %in% c("DYNC2H1")] <- 1
EntireCohort_pathvar$DYNC2H1[EntireCohort_pathvar$GENE2 %in% c("DYNC2H1")] <- 1
EntireCohort_pathvar$DYNC2H1[EntireCohort_pathvar$GENE3 %in% c("DYNC2H1")] <- 1
EntireCohort_pathvar$DYNC2H1[EntireCohort_pathvar$GENE4 %in% c("DYNC2H1")] <- 1
EntireCohort_pathvar$DYNC2H1[EntireCohort_pathvar$GENE5 %in% c("DYNC2H1")] <- 1
EntireCohort_pathvar$DYNC2H1[EntireCohort_pathvar$GENE6 %in% c("DYNC2H1")] <- 1
EntireCohort_pathvar$DYNC2H1[EntireCohort_pathvar$GENE7 %in% c("DYNC2H1")] <- 1
EntireCohort_pathvar$DYNC2H1[EntireCohort_pathvar$GENE8 %in% c("DYNC2H1")] <- 1
EntireCohort_pathvar$DYNC2H1[EntireCohort_pathvar$GENE9 %in% c("DYNC2H1")] <- 1

EntireCohort_pathvar$DYNC2H1[is.na(EntireCohort_pathvar$DYNC2H1)] <- 0
EntireCohort_pathvar$DYNC2H1 <- as.numeric(EntireCohort_pathvar$DYNC2H1)
unique(EntireCohort_pathvar$DYNC2H1)
# 0 1

table(EntireCohort_pathvar$DYNC2H1)
#  0      1 
# 462465   7334





### 'SGMS2'
EntireCohort_pathvar$SGMS2 <- NA
EntireCohort_pathvar$SGMS2[EntireCohort_pathvar$GENE1 %in% c("SGMS2")] <- 1
EntireCohort_pathvar$SGMS2[EntireCohort_pathvar$GENE2 %in% c("SGMS2")] <- 1
EntireCohort_pathvar$SGMS2[EntireCohort_pathvar$GENE3 %in% c("SGMS2")] <- 1
EntireCohort_pathvar$SGMS2[EntireCohort_pathvar$GENE4 %in% c("SGMS2")] <- 1
EntireCohort_pathvar$SGMS2[EntireCohort_pathvar$GENE5 %in% c("SGMS2")] <- 1
EntireCohort_pathvar$SGMS2[EntireCohort_pathvar$GENE6 %in% c("SGMS2")] <- 1
EntireCohort_pathvar$SGMS2[EntireCohort_pathvar$GENE7 %in% c("SGMS2")] <- 1
EntireCohort_pathvar$SGMS2[EntireCohort_pathvar$GENE8 %in% c("SGMS2")] <- 1
EntireCohort_pathvar$SGMS2[EntireCohort_pathvar$GENE9 %in% c("SGMS2")] <- 1

EntireCohort_pathvar$SGMS2[is.na(EntireCohort_pathvar$SGMS2)] <- 0
EntireCohort_pathvar$SGMS2 <- as.numeric(EntireCohort_pathvar$SGMS2)
unique(EntireCohort_pathvar$SGMS2)
# 0 1

table(EntireCohort_pathvar$SGMS2)
#   0      1 
# 469797      2 
## only 2 participants - so won't include this gene in the GOI analyses, 
## but the 2 participants will still be included in Bonegene



### 'TMEM38B'
EntireCohort_pathvar$TMEM38B <- NA
EntireCohort_pathvar$TMEM38B[EntireCohort_pathvar$GENE1 %in% c("TMEM38B")] <- 1
EntireCohort_pathvar$TMEM38B[EntireCohort_pathvar$GENE2 %in% c("TMEM38B")] <- 1
EntireCohort_pathvar$TMEM38B[EntireCohort_pathvar$GENE3 %in% c("TMEM38B")] <- 1
EntireCohort_pathvar$TMEM38B[EntireCohort_pathvar$GENE4 %in% c("TMEM38B")] <- 1
EntireCohort_pathvar$TMEM38B[EntireCohort_pathvar$GENE5 %in% c("TMEM38B")] <- 1
EntireCohort_pathvar$TMEM38B[EntireCohort_pathvar$GENE6 %in% c("TMEM38B")] <- 1
EntireCohort_pathvar$TMEM38B[EntireCohort_pathvar$GENE7 %in% c("TMEM38B")] <- 1
EntireCohort_pathvar$TMEM38B[EntireCohort_pathvar$GENE8 %in% c("TMEM38B")] <- 1
EntireCohort_pathvar$TMEM38B[EntireCohort_pathvar$GENE9 %in% c("TMEM38B")] <- 1

EntireCohort_pathvar$TMEM38B[is.na(EntireCohort_pathvar$TMEM38B)] <- 0
EntireCohort_pathvar$TMEM38B <- as.numeric(EntireCohort_pathvar$TMEM38B)
unique(EntireCohort_pathvar$TMEM38B)
# 0 1

table(EntireCohort_pathvar$TMEM38B)
#  0      1 
# 469509    290 




### 'TRPV4'
EntireCohort_pathvar$TRPV4 <- NA
EntireCohort_pathvar$TRPV4[EntireCohort_pathvar$GENE1 %in% c("TRPV4")] <- 1
EntireCohort_pathvar$TRPV4[EntireCohort_pathvar$GENE2 %in% c("TRPV4")] <- 1
EntireCohort_pathvar$TRPV4[EntireCohort_pathvar$GENE3 %in% c("TRPV4")] <- 1
EntireCohort_pathvar$TRPV4[EntireCohort_pathvar$GENE4 %in% c("TRPV4")] <- 1
EntireCohort_pathvar$TRPV4[EntireCohort_pathvar$GENE5 %in% c("TRPV4")] <- 1
EntireCohort_pathvar$TRPV4[EntireCohort_pathvar$GENE6 %in% c("TRPV4")] <- 1
EntireCohort_pathvar$TRPV4[EntireCohort_pathvar$GENE7 %in% c("TRPV4")] <- 1
EntireCohort_pathvar$TRPV4[EntireCohort_pathvar$GENE8 %in% c("TRPV4")] <- 1
EntireCohort_pathvar$TRPV4[EntireCohort_pathvar$GENE9 %in% c("TRPV4")] <- 1

EntireCohort_pathvar$TRPV4[is.na(EntireCohort_pathvar$TRPV4)] <- 0
EntireCohort_pathvar$TRPV4 <- as.numeric(EntireCohort_pathvar$TRPV4)
unique(EntireCohort_pathvar$TRPV4)
# 0 1

table(EntireCohort_pathvar$TRPV4)
#  0      1 
# 469795      4




### 'ARSL'
EntireCohort_pathvar$ARSL <- NA
EntireCohort_pathvar$ARSL[EntireCohort_pathvar$GENE1 %in% c("ARSL")] <- 1
EntireCohort_pathvar$ARSL[EntireCohort_pathvar$GENE2 %in% c("ARSL")] <- 1
EntireCohort_pathvar$ARSL[EntireCohort_pathvar$GENE3 %in% c("ARSL")] <- 1
EntireCohort_pathvar$ARSL[EntireCohort_pathvar$GENE4 %in% c("ARSL")] <- 1
EntireCohort_pathvar$ARSL[EntireCohort_pathvar$GENE5 %in% c("ARSL")] <- 1
EntireCohort_pathvar$ARSL[EntireCohort_pathvar$GENE6 %in% c("ARSL")] <- 1
EntireCohort_pathvar$ARSL[EntireCohort_pathvar$GENE7 %in% c("ARSL")] <- 1
EntireCohort_pathvar$ARSL[EntireCohort_pathvar$GENE8 %in% c("ARSL")] <- 1
EntireCohort_pathvar$ARSL[EntireCohort_pathvar$GENE9 %in% c("ARSL")] <- 1

EntireCohort_pathvar$ARSL[is.na(EntireCohort_pathvar$ARSL)] <- 0
EntireCohort_pathvar$ARSL <- as.numeric(EntireCohort_pathvar$ARSL)
unique(EntireCohort_pathvar$ARSL)
# 0 1

table(EntireCohort_pathvar$ARSL)
#  0      1 
# 469786     13




### 'HSPG2'
EntireCohort_pathvar$HSPG2 <- NA
EntireCohort_pathvar$HSPG2[EntireCohort_pathvar$GENE1 %in% c("HSPG2")] <- 1
EntireCohort_pathvar$HSPG2[EntireCohort_pathvar$GENE2 %in% c("HSPG2")] <- 1
EntireCohort_pathvar$HSPG2[EntireCohort_pathvar$GENE3 %in% c("HSPG2")] <- 1
EntireCohort_pathvar$HSPG2[EntireCohort_pathvar$GENE4 %in% c("HSPG2")] <- 1
EntireCohort_pathvar$HSPG2[EntireCohort_pathvar$GENE5 %in% c("HSPG2")] <- 1
EntireCohort_pathvar$HSPG2[EntireCohort_pathvar$GENE6 %in% c("HSPG2")] <- 1
EntireCohort_pathvar$HSPG2[EntireCohort_pathvar$GENE7 %in% c("HSPG2")] <- 1
EntireCohort_pathvar$HSPG2[EntireCohort_pathvar$GENE8 %in% c("HSPG2")] <- 1
EntireCohort_pathvar$HSPG2[EntireCohort_pathvar$GENE9 %in% c("HSPG2")] <- 1

EntireCohort_pathvar$HSPG2[is.na(EntireCohort_pathvar$HSPG2)] <- 0
EntireCohort_pathvar$HSPG2 <- as.numeric(EntireCohort_pathvar$HSPG2)
unique(EntireCohort_pathvar$HSPG2)
# 0 1

table(EntireCohort_pathvar$HSPG2)
#  0      1 
# 460241   9558




### 'RSPRY1'
EntireCohort_pathvar$RSPRY1 <- NA
EntireCohort_pathvar$RSPRY1[EntireCohort_pathvar$GENE1 %in% c("RSPRY1")] <- 1
EntireCohort_pathvar$RSPRY1[EntireCohort_pathvar$GENE2 %in% c("RSPRY1")] <- 1
EntireCohort_pathvar$RSPRY1[EntireCohort_pathvar$GENE3 %in% c("RSPRY1")] <- 1
EntireCohort_pathvar$RSPRY1[EntireCohort_pathvar$GENE4 %in% c("RSPRY1")] <- 1
EntireCohort_pathvar$RSPRY1[EntireCohort_pathvar$GENE5 %in% c("RSPRY1")] <- 1
EntireCohort_pathvar$RSPRY1[EntireCohort_pathvar$GENE6 %in% c("RSPRY1")] <- 1
EntireCohort_pathvar$RSPRY1[EntireCohort_pathvar$GENE7 %in% c("RSPRY1")] <- 1
EntireCohort_pathvar$RSPRY1[EntireCohort_pathvar$GENE8 %in% c("RSPRY1")] <- 1
EntireCohort_pathvar$RSPRY1[EntireCohort_pathvar$GENE9 %in% c("RSPRY1")] <- 1

EntireCohort_pathvar$RSPRY1[is.na(EntireCohort_pathvar$RSPRY1)] <- 0
EntireCohort_pathvar$RSPRY1 <- as.numeric(EntireCohort_pathvar$RSPRY1)
unique(EntireCohort_pathvar$RSPRY1)
# 0 1

table(EntireCohort_pathvar$RSPRY1)
#  0      1 
# 469376    423




### 'PIEZO2'
# no variants present




### 'CANT1'
EntireCohort_pathvar$CANT1 <- NA
EntireCohort_pathvar$CANT1[EntireCohort_pathvar$GENE1 %in% c("CANT1")] <- 1
EntireCohort_pathvar$CANT1[EntireCohort_pathvar$GENE2 %in% c("CANT1")] <- 1
EntireCohort_pathvar$CANT1[EntireCohort_pathvar$GENE3 %in% c("CANT1")] <- 1
EntireCohort_pathvar$CANT1[EntireCohort_pathvar$GENE4 %in% c("CANT1")] <- 1
EntireCohort_pathvar$CANT1[EntireCohort_pathvar$GENE5 %in% c("CANT1")] <- 1
EntireCohort_pathvar$CANT1[EntireCohort_pathvar$GENE6 %in% c("CANT1")] <- 1
EntireCohort_pathvar$CANT1[EntireCohort_pathvar$GENE7 %in% c("CANT1")] <- 1
EntireCohort_pathvar$CANT1[EntireCohort_pathvar$GENE8 %in% c("CANT1")] <- 1
EntireCohort_pathvar$CANT1[EntireCohort_pathvar$GENE9 %in% c("CANT1")] <- 1

EntireCohort_pathvar$CANT1[is.na(EntireCohort_pathvar$CANT1)] <- 0
EntireCohort_pathvar$CANT1 <- as.numeric(EntireCohort_pathvar$CANT1)
unique(EntireCohort_pathvar$CANT1)
# 0 1

table(EntireCohort_pathvar$CANT1)
#  0      1 
# 469023    776





### 'NOTCH2'
# no variants present




### 'ACP5'
EntireCohort_pathvar$ACP5 <- NA
EntireCohort_pathvar$ACP5[EntireCohort_pathvar$GENE1 %in% c("ACP5")] <- 1
EntireCohort_pathvar$ACP5[EntireCohort_pathvar$GENE2 %in% c("ACP5")] <- 1
EntireCohort_pathvar$ACP5[EntireCohort_pathvar$GENE3 %in% c("ACP5")] <- 1
EntireCohort_pathvar$ACP5[EntireCohort_pathvar$GENE4 %in% c("ACP5")] <- 1
EntireCohort_pathvar$ACP5[EntireCohort_pathvar$GENE5 %in% c("ACP5")] <- 1
EntireCohort_pathvar$ACP5[EntireCohort_pathvar$GENE6 %in% c("ACP5")] <- 1
EntireCohort_pathvar$ACP5[EntireCohort_pathvar$GENE7 %in% c("ACP5")] <- 1
EntireCohort_pathvar$ACP5[EntireCohort_pathvar$GENE8 %in% c("ACP5")] <- 1
EntireCohort_pathvar$ACP5[EntireCohort_pathvar$GENE9 %in% c("ACP5")] <- 1

EntireCohort_pathvar$ACP5[is.na(EntireCohort_pathvar$ACP5)] <- 0
EntireCohort_pathvar$ACP5 <- as.numeric(EntireCohort_pathvar$ACP5)
unique(EntireCohort_pathvar$ACP5)
# 0 1

table(EntireCohort_pathvar$ACP5)
#  0      1 
# 469369    430 





### 'DDR2'
EntireCohort_pathvar$DDR2 <- NA
EntireCohort_pathvar$DDR2[EntireCohort_pathvar$GENE1 %in% c("DDR2")] <- 1
EntireCohort_pathvar$DDR2[EntireCohort_pathvar$GENE2 %in% c("DDR2")] <- 1
EntireCohort_pathvar$DDR2[EntireCohort_pathvar$GENE3 %in% c("DDR2")] <- 1
EntireCohort_pathvar$DDR2[EntireCohort_pathvar$GENE4 %in% c("DDR2")] <- 1
EntireCohort_pathvar$DDR2[EntireCohort_pathvar$GENE5 %in% c("DDR2")] <- 1
EntireCohort_pathvar$DDR2[EntireCohort_pathvar$GENE6 %in% c("DDR2")] <- 1
EntireCohort_pathvar$DDR2[EntireCohort_pathvar$GENE7 %in% c("DDR2")] <- 1
EntireCohort_pathvar$DDR2[EntireCohort_pathvar$GENE8 %in% c("DDR2")] <- 1
EntireCohort_pathvar$DDR2[EntireCohort_pathvar$GENE9 %in% c("DDR2")] <- 1

EntireCohort_pathvar$DDR2[is.na(EntireCohort_pathvar$DDR2)] <- 0
EntireCohort_pathvar$DDR2 <- as.numeric(EntireCohort_pathvar$DDR2)
unique(EntireCohort_pathvar$DDR2)
# 0 1

table(EntireCohort_pathvar$DDR2)
#   0      1 
# 468104   1695 





### 'INPPL1'
EntireCohort_pathvar$INPPL1 <- NA
EntireCohort_pathvar$INPPL1[EntireCohort_pathvar$GENE1 %in% c("INPPL1")] <- 1
EntireCohort_pathvar$INPPL1[EntireCohort_pathvar$GENE2 %in% c("INPPL1")] <- 1
EntireCohort_pathvar$INPPL1[EntireCohort_pathvar$GENE3 %in% c("INPPL1")] <- 1
EntireCohort_pathvar$INPPL1[EntireCohort_pathvar$GENE4 %in% c("INPPL1")] <- 1
EntireCohort_pathvar$INPPL1[EntireCohort_pathvar$GENE5 %in% c("INPPL1")] <- 1
EntireCohort_pathvar$INPPL1[EntireCohort_pathvar$GENE6 %in% c("INPPL1")] <- 1
EntireCohort_pathvar$INPPL1[EntireCohort_pathvar$GENE7 %in% c("INPPL1")] <- 1
EntireCohort_pathvar$INPPL1[EntireCohort_pathvar$GENE8 %in% c("INPPL1")] <- 1
EntireCohort_pathvar$INPPL1[EntireCohort_pathvar$GENE9 %in% c("INPPL1")] <- 1

EntireCohort_pathvar$INPPL1[is.na(EntireCohort_pathvar$INPPL1)] <- 0
EntireCohort_pathvar$INPPL1 <- as.numeric(EntireCohort_pathvar$INPPL1)
unique(EntireCohort_pathvar$INPPL1)
# 0 1

table(EntireCohort_pathvar$INPPL1)
#    0      1 
# 467440   2359




### 'LIFR'
EntireCohort_pathvar$LIFR <- NA
EntireCohort_pathvar$LIFR[EntireCohort_pathvar$GENE1 %in% c("LIFR")] <- 1
EntireCohort_pathvar$LIFR[EntireCohort_pathvar$GENE2 %in% c("LIFR")] <- 1
EntireCohort_pathvar$LIFR[EntireCohort_pathvar$GENE3 %in% c("LIFR")] <- 1
EntireCohort_pathvar$LIFR[EntireCohort_pathvar$GENE4 %in% c("LIFR")] <- 1
EntireCohort_pathvar$LIFR[EntireCohort_pathvar$GENE5 %in% c("LIFR")] <- 1
EntireCohort_pathvar$LIFR[EntireCohort_pathvar$GENE6 %in% c("LIFR")] <- 1
EntireCohort_pathvar$LIFR[EntireCohort_pathvar$GENE7 %in% c("LIFR")] <- 1
EntireCohort_pathvar$LIFR[EntireCohort_pathvar$GENE8 %in% c("LIFR")] <- 1
EntireCohort_pathvar$LIFR[EntireCohort_pathvar$GENE9 %in% c("LIFR")] <- 1

EntireCohort_pathvar$LIFR[is.na(EntireCohort_pathvar$LIFR)] <- 0
EntireCohort_pathvar$LIFR <- as.numeric(EntireCohort_pathvar$LIFR)
unique(EntireCohort_pathvar$LIFR)
# 0 1

table(EntireCohort_pathvar$LIFR)
#  0      1 
# 468357   1442 




### 'NPPC'
EntireCohort_pathvar$NPPC <- NA
EntireCohort_pathvar$NPPC[EntireCohort_pathvar$GENE1 %in% c("NPPC")] <- 1
EntireCohort_pathvar$NPPC[EntireCohort_pathvar$GENE2 %in% c("NPPC")] <- 1
EntireCohort_pathvar$NPPC[EntireCohort_pathvar$GENE3 %in% c("NPPC")] <- 1
EntireCohort_pathvar$NPPC[EntireCohort_pathvar$GENE4 %in% c("NPPC")] <- 1
EntireCohort_pathvar$NPPC[EntireCohort_pathvar$GENE5 %in% c("NPPC")] <- 1
EntireCohort_pathvar$NPPC[EntireCohort_pathvar$GENE6 %in% c("NPPC")] <- 1
EntireCohort_pathvar$NPPC[EntireCohort_pathvar$GENE7 %in% c("NPPC")] <- 1
EntireCohort_pathvar$NPPC[EntireCohort_pathvar$GENE8 %in% c("NPPC")] <- 1
EntireCohort_pathvar$NPPC[EntireCohort_pathvar$GENE9 %in% c("NPPC")] <- 1

EntireCohort_pathvar$NPPC[is.na(EntireCohort_pathvar$NPPC)] <- 0
EntireCohort_pathvar$NPPC <- as.numeric(EntireCohort_pathvar$NPPC)
unique(EntireCohort_pathvar$NPPC)
# 0 1

table(EntireCohort_pathvar$NPPC)
#  0      1 
# 469754     45




### 'NPR2'
EntireCohort_pathvar$NPR2 <- NA
EntireCohort_pathvar$NPR2[EntireCohort_pathvar$GENE1 %in% c("NPR2")] <- 1
EntireCohort_pathvar$NPR2[EntireCohort_pathvar$GENE2 %in% c("NPR2")] <- 1
EntireCohort_pathvar$NPR2[EntireCohort_pathvar$GENE3 %in% c("NPR2")] <- 1
EntireCohort_pathvar$NPR2[EntireCohort_pathvar$GENE4 %in% c("NPR2")] <- 1
EntireCohort_pathvar$NPR2[EntireCohort_pathvar$GENE5 %in% c("NPR2")] <- 1
EntireCohort_pathvar$NPR2[EntireCohort_pathvar$GENE6 %in% c("NPR2")] <- 1
EntireCohort_pathvar$NPR2[EntireCohort_pathvar$GENE7 %in% c("NPR2")] <- 1
EntireCohort_pathvar$NPR2[EntireCohort_pathvar$GENE8 %in% c("NPR2")] <- 1
EntireCohort_pathvar$NPR2[EntireCohort_pathvar$GENE9 %in% c("NPR2")] <- 1

EntireCohort_pathvar$NPR2[is.na(EntireCohort_pathvar$NPR2)] <- 0
EntireCohort_pathvar$NPR2 <- as.numeric(EntireCohort_pathvar$NPR2)
unique(EntireCohort_pathvar$NPR2)
# 0 1

table(EntireCohort_pathvar$NPR2)
#  0      1 
# 467995   1804




### 'PRKG2'
EntireCohort_pathvar$PRKG2 <- NA
EntireCohort_pathvar$PRKG2[EntireCohort_pathvar$GENE1 %in% c("PRKG2")] <- 1
EntireCohort_pathvar$PRKG2[EntireCohort_pathvar$GENE2 %in% c("PRKG2")] <- 1
EntireCohort_pathvar$PRKG2[EntireCohort_pathvar$GENE3 %in% c("PRKG2")] <- 1
EntireCohort_pathvar$PRKG2[EntireCohort_pathvar$GENE4 %in% c("PRKG2")] <- 1
EntireCohort_pathvar$PRKG2[EntireCohort_pathvar$GENE5 %in% c("PRKG2")] <- 1
EntireCohort_pathvar$PRKG2[EntireCohort_pathvar$GENE6 %in% c("PRKG2")] <- 1
EntireCohort_pathvar$PRKG2[EntireCohort_pathvar$GENE7 %in% c("PRKG2")] <- 1
EntireCohort_pathvar$PRKG2[EntireCohort_pathvar$GENE8 %in% c("PRKG2")] <- 1
EntireCohort_pathvar$PRKG2[EntireCohort_pathvar$GENE9 %in% c("PRKG2")] <- 1

EntireCohort_pathvar$PRKG2[is.na(EntireCohort_pathvar$PRKG2)] <- 0
EntireCohort_pathvar$PRKG2 <- as.numeric(EntireCohort_pathvar$PRKG2)
unique(EntireCohort_pathvar$PRKG2)
# 0 1

table(EntireCohort_pathvar$PRKG2)
#  0      1 
# 469460    339 




### 'PPP3CA'
# no variants present





# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)




# Braingenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

# go thr all the genes to check that they're present
grep("NUP37", EntireCohort_pathvar$GENE)
# all seem to be present :)


# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$Braingene <- NA
EntireCohort_pathvar$Braingene[EntireCohort_pathvar$GENE1 %in% c('NUP133', 'DYNC1I2', 'CASK', 'SLC1A4', 'MFSD2A', 'PPFIBP1', 'PCLO', 'PCDH12', 'OCLN', 'FRMD4A', 'NUP37')] <- 1
EntireCohort_pathvar$Braingene[EntireCohort_pathvar$GENE2 %in% c('NUP133', 'DYNC1I2', 'CASK', 'SLC1A4', 'MFSD2A', 'PPFIBP1', 'PCLO', 'PCDH12', 'OCLN', 'FRMD4A', 'NUP37')] <- 1
EntireCohort_pathvar$Braingene[EntireCohort_pathvar$GENE3 %in% c('NUP133', 'DYNC1I2', 'CASK', 'SLC1A4', 'MFSD2A', 'PPFIBP1', 'PCLO', 'PCDH12', 'OCLN', 'FRMD4A', 'NUP37')] <- 1
EntireCohort_pathvar$Braingene[EntireCohort_pathvar$GENE4 %in% c('NUP133', 'DYNC1I2', 'CASK', 'SLC1A4', 'MFSD2A', 'PPFIBP1', 'PCLO', 'PCDH12', 'OCLN', 'FRMD4A', 'NUP37')] <- 1
EntireCohort_pathvar$Braingene[EntireCohort_pathvar$GENE5 %in% c('NUP133', 'DYNC1I2', 'CASK', 'SLC1A4', 'MFSD2A', 'PPFIBP1', 'PCLO', 'PCDH12', 'OCLN', 'FRMD4A', 'NUP37')] <- 1
EntireCohort_pathvar$Braingene[EntireCohort_pathvar$GENE6 %in% c('NUP133', 'DYNC1I2', 'CASK', 'SLC1A4', 'MFSD2A', 'PPFIBP1', 'PCLO', 'PCDH12', 'OCLN', 'FRMD4A', 'NUP37')] <- 1
EntireCohort_pathvar$Braingene[EntireCohort_pathvar$GENE7 %in% c('NUP133', 'DYNC1I2', 'CASK', 'SLC1A4', 'MFSD2A', 'PPFIBP1', 'PCLO', 'PCDH12', 'OCLN', 'FRMD4A', 'NUP37')] <- 1
EntireCohort_pathvar$Braingene[EntireCohort_pathvar$GENE8 %in% c('NUP133', 'DYNC1I2', 'CASK', 'SLC1A4', 'MFSD2A', 'PPFIBP1', 'PCLO', 'PCDH12', 'OCLN', 'FRMD4A', 'NUP37')] <- 1
EntireCohort_pathvar$Braingene[EntireCohort_pathvar$GENE9 %in% c('NUP133', 'DYNC1I2', 'CASK', 'SLC1A4', 'MFSD2A', 'PPFIBP1', 'PCLO', 'PCDH12', 'OCLN', 'FRMD4A', 'NUP37')] <- 1


EntireCohort_pathvar$Braingene[is.na(EntireCohort_pathvar$Braingene)] <- 0
EntireCohort_pathvar$Braingene <- as.numeric(EntireCohort_pathvar$Braingene)
unique(EntireCohort_pathvar$Braingene)
# 0 1

table(EntireCohort_pathvar$Braingene)
# 0       1 
# 445031  24768




# now to make a column for each of the Braingenes
### 'NUP133'
EntireCohort_pathvar$NUP133 <- NA
EntireCohort_pathvar$NUP133[EntireCohort_pathvar$GENE1 %in% c("NUP133")] <- 1
EntireCohort_pathvar$NUP133[EntireCohort_pathvar$GENE2 %in% c("NUP133")] <- 1
EntireCohort_pathvar$NUP133[EntireCohort_pathvar$GENE3 %in% c("NUP133")] <- 1
EntireCohort_pathvar$NUP133[EntireCohort_pathvar$GENE4 %in% c("NUP133")] <- 1
EntireCohort_pathvar$NUP133[EntireCohort_pathvar$GENE5 %in% c("NUP133")] <- 1
EntireCohort_pathvar$NUP133[EntireCohort_pathvar$GENE6 %in% c("NUP133")] <- 1
EntireCohort_pathvar$NUP133[EntireCohort_pathvar$GENE7 %in% c("NUP133")] <- 1
EntireCohort_pathvar$NUP133[EntireCohort_pathvar$GENE8 %in% c("NUP133")] <- 1
EntireCohort_pathvar$NUP133[EntireCohort_pathvar$GENE9 %in% c("NUP133")] <- 1

EntireCohort_pathvar$NUP133[is.na(EntireCohort_pathvar$NUP133)] <- 0
EntireCohort_pathvar$NUP133 <- as.numeric(EntireCohort_pathvar$NUP133)
unique(EntireCohort_pathvar$NUP133)
# 0 1

table(EntireCohort_pathvar$NUP133)
#  0      1 
# 468750   1049





### 'DYNC1I2'
EntireCohort_pathvar$DYNC1I2 <- NA
EntireCohort_pathvar$DYNC1I2[EntireCohort_pathvar$GENE1 %in% c("DYNC1I2")] <- 1
EntireCohort_pathvar$DYNC1I2[EntireCohort_pathvar$GENE2 %in% c("DYNC1I2")] <- 1
EntireCohort_pathvar$DYNC1I2[EntireCohort_pathvar$GENE3 %in% c("DYNC1I2")] <- 1
EntireCohort_pathvar$DYNC1I2[EntireCohort_pathvar$GENE4 %in% c("DYNC1I2")] <- 1
EntireCohort_pathvar$DYNC1I2[EntireCohort_pathvar$GENE5 %in% c("DYNC1I2")] <- 1
EntireCohort_pathvar$DYNC1I2[EntireCohort_pathvar$GENE6 %in% c("DYNC1I2")] <- 1
EntireCohort_pathvar$DYNC1I2[EntireCohort_pathvar$GENE7 %in% c("DYNC1I2")] <- 1
EntireCohort_pathvar$DYNC1I2[EntireCohort_pathvar$GENE8 %in% c("DYNC1I2")] <- 1
EntireCohort_pathvar$DYNC1I2[EntireCohort_pathvar$GENE9 %in% c("DYNC1I2")] <- 1

EntireCohort_pathvar$DYNC1I2[is.na(EntireCohort_pathvar$DYNC1I2)] <- 0
EntireCohort_pathvar$DYNC1I2 <- as.numeric(EntireCohort_pathvar$DYNC1I2)
unique(EntireCohort_pathvar$DYNC1I2)
# 0 1

table(EntireCohort_pathvar$DYNC1I2)
#  0      1 
# 468616   1183




### 'CASK'
EntireCohort_pathvar$CASK <- NA
EntireCohort_pathvar$CASK[EntireCohort_pathvar$GENE1 %in% c("CASK")] <- 1
EntireCohort_pathvar$CASK[EntireCohort_pathvar$GENE2 %in% c("CASK")] <- 1
EntireCohort_pathvar$CASK[EntireCohort_pathvar$GENE3 %in% c("CASK")] <- 1
EntireCohort_pathvar$CASK[EntireCohort_pathvar$GENE4 %in% c("CASK")] <- 1
EntireCohort_pathvar$CASK[EntireCohort_pathvar$GENE5 %in% c("CASK")] <- 1
EntireCohort_pathvar$CASK[EntireCohort_pathvar$GENE6 %in% c("CASK")] <- 1
EntireCohort_pathvar$CASK[EntireCohort_pathvar$GENE7 %in% c("CASK")] <- 1
EntireCohort_pathvar$CASK[EntireCohort_pathvar$GENE8 %in% c("CASK")] <- 1
EntireCohort_pathvar$CASK[EntireCohort_pathvar$GENE9 %in% c("CASK")] <- 1

EntireCohort_pathvar$CASK[is.na(EntireCohort_pathvar$CASK)] <- 0
EntireCohort_pathvar$CASK <- as.numeric(EntireCohort_pathvar$CASK)
unique(EntireCohort_pathvar$CASK)
# 0 1

table(EntireCohort_pathvar$CASK)
#  0      1 
# 469445    354 




### 'SLC1A4'
EntireCohort_pathvar$SLC1A4 <- NA
EntireCohort_pathvar$SLC1A4[EntireCohort_pathvar$GENE1 %in% c("SLC1A4")] <- 1
EntireCohort_pathvar$SLC1A4[EntireCohort_pathvar$GENE2 %in% c("SLC1A4")] <- 1
EntireCohort_pathvar$SLC1A4[EntireCohort_pathvar$GENE3 %in% c("SLC1A4")] <- 1
EntireCohort_pathvar$SLC1A4[EntireCohort_pathvar$GENE4 %in% c("SLC1A4")] <- 1
EntireCohort_pathvar$SLC1A4[EntireCohort_pathvar$GENE5 %in% c("SLC1A4")] <- 1
EntireCohort_pathvar$SLC1A4[EntireCohort_pathvar$GENE6 %in% c("SLC1A4")] <- 1
EntireCohort_pathvar$SLC1A4[EntireCohort_pathvar$GENE7 %in% c("SLC1A4")] <- 1
EntireCohort_pathvar$SLC1A4[EntireCohort_pathvar$GENE8 %in% c("SLC1A4")] <- 1
EntireCohort_pathvar$SLC1A4[EntireCohort_pathvar$GENE9 %in% c("SLC1A4")] <- 1

EntireCohort_pathvar$SLC1A4[is.na(EntireCohort_pathvar$SLC1A4)] <- 0
EntireCohort_pathvar$SLC1A4 <- as.numeric(EntireCohort_pathvar$SLC1A4)
unique(EntireCohort_pathvar$SLC1A4)
# 0 1

table(EntireCohort_pathvar$SLC1A4)
#  0      1 
# 467574   2225




### 'MFSD2A'
EntireCohort_pathvar$MFSD2A <- NA
EntireCohort_pathvar$MFSD2A[EntireCohort_pathvar$GENE1 %in% c("MFSD2A")] <- 1
EntireCohort_pathvar$MFSD2A[EntireCohort_pathvar$GENE2 %in% c("MFSD2A")] <- 1
EntireCohort_pathvar$MFSD2A[EntireCohort_pathvar$GENE3 %in% c("MFSD2A")] <- 1
EntireCohort_pathvar$MFSD2A[EntireCohort_pathvar$GENE4 %in% c("MFSD2A")] <- 1
EntireCohort_pathvar$MFSD2A[EntireCohort_pathvar$GENE5 %in% c("MFSD2A")] <- 1
EntireCohort_pathvar$MFSD2A[EntireCohort_pathvar$GENE6 %in% c("MFSD2A")] <- 1
EntireCohort_pathvar$MFSD2A[EntireCohort_pathvar$GENE7 %in% c("MFSD2A")] <- 1
EntireCohort_pathvar$MFSD2A[EntireCohort_pathvar$GENE8 %in% c("MFSD2A")] <- 1
EntireCohort_pathvar$MFSD2A[EntireCohort_pathvar$GENE9 %in% c("MFSD2A")] <- 1

EntireCohort_pathvar$MFSD2A[is.na(EntireCohort_pathvar$MFSD2A)] <- 0
EntireCohort_pathvar$MFSD2A <- as.numeric(EntireCohort_pathvar$MFSD2A)
unique(EntireCohort_pathvar$MFSD2A)
# 0 1

table(EntireCohort_pathvar$MFSD2A)
#   0      1 
# 468784   1015




### 'PPFIBP1'
EntireCohort_pathvar$PPFIBP1 <- NA
EntireCohort_pathvar$PPFIBP1[EntireCohort_pathvar$GENE1 %in% c("PPFIBP1")] <- 1
EntireCohort_pathvar$PPFIBP1[EntireCohort_pathvar$GENE2 %in% c("PPFIBP1")] <- 1
EntireCohort_pathvar$PPFIBP1[EntireCohort_pathvar$GENE3 %in% c("PPFIBP1")] <- 1
EntireCohort_pathvar$PPFIBP1[EntireCohort_pathvar$GENE4 %in% c("PPFIBP1")] <- 1
EntireCohort_pathvar$PPFIBP1[EntireCohort_pathvar$GENE5 %in% c("PPFIBP1")] <- 1
EntireCohort_pathvar$PPFIBP1[EntireCohort_pathvar$GENE6 %in% c("PPFIBP1")] <- 1
EntireCohort_pathvar$PPFIBP1[EntireCohort_pathvar$GENE7 %in% c("PPFIBP1")] <- 1
EntireCohort_pathvar$PPFIBP1[EntireCohort_pathvar$GENE8 %in% c("PPFIBP1")] <- 1
EntireCohort_pathvar$PPFIBP1[EntireCohort_pathvar$GENE9 %in% c("PPFIBP1")] <- 1

EntireCohort_pathvar$PPFIBP1[is.na(EntireCohort_pathvar$PPFIBP1)] <- 0
EntireCohort_pathvar$PPFIBP1 <- as.numeric(EntireCohort_pathvar$PPFIBP1)
unique(EntireCohort_pathvar$PPFIBP1)
# 0 1

table(EntireCohort_pathvar$PPFIBP1)
#  0      1 
# 466140   3659 




### 'PCLO'
EntireCohort_pathvar$PCLO <- NA
EntireCohort_pathvar$PCLO[EntireCohort_pathvar$GENE1 %in% c("PCLO")] <- 1
EntireCohort_pathvar$PCLO[EntireCohort_pathvar$GENE2 %in% c("PCLO")] <- 1
EntireCohort_pathvar$PCLO[EntireCohort_pathvar$GENE3 %in% c("PCLO")] <- 1
EntireCohort_pathvar$PCLO[EntireCohort_pathvar$GENE4 %in% c("PCLO")] <- 1
EntireCohort_pathvar$PCLO[EntireCohort_pathvar$GENE5 %in% c("PCLO")] <- 1
EntireCohort_pathvar$PCLO[EntireCohort_pathvar$GENE6 %in% c("PCLO")] <- 1
EntireCohort_pathvar$PCLO[EntireCohort_pathvar$GENE7 %in% c("PCLO")] <- 1
EntireCohort_pathvar$PCLO[EntireCohort_pathvar$GENE8 %in% c("PCLO")] <- 1
EntireCohort_pathvar$PCLO[EntireCohort_pathvar$GENE9 %in% c("PCLO")] <- 1

EntireCohort_pathvar$PCLO[is.na(EntireCohort_pathvar$PCLO)] <- 0
EntireCohort_pathvar$PCLO <- as.numeric(EntireCohort_pathvar$PCLO)
unique(EntireCohort_pathvar$PCLO)
# 0 1

table(EntireCohort_pathvar$PCLO)
#  0      1 
# 465984   3815




### 'PCDH12'
EntireCohort_pathvar$PCDH12 <- NA
EntireCohort_pathvar$PCDH12[EntireCohort_pathvar$GENE1 %in% c("PCDH12")] <- 1
EntireCohort_pathvar$PCDH12[EntireCohort_pathvar$GENE2 %in% c("PCDH12")] <- 1
EntireCohort_pathvar$PCDH12[EntireCohort_pathvar$GENE3 %in% c("PCDH12")] <- 1
EntireCohort_pathvar$PCDH12[EntireCohort_pathvar$GENE4 %in% c("PCDH12")] <- 1
EntireCohort_pathvar$PCDH12[EntireCohort_pathvar$GENE5 %in% c("PCDH12")] <- 1
EntireCohort_pathvar$PCDH12[EntireCohort_pathvar$GENE6 %in% c("PCDH12")] <- 1
EntireCohort_pathvar$PCDH12[EntireCohort_pathvar$GENE7 %in% c("PCDH12")] <- 1
EntireCohort_pathvar$PCDH12[EntireCohort_pathvar$GENE8 %in% c("PCDH12")] <- 1
EntireCohort_pathvar$PCDH12[EntireCohort_pathvar$GENE9 %in% c("PCDH12")] <- 1

EntireCohort_pathvar$PCDH12[is.na(EntireCohort_pathvar$PCDH12)] <- 0
EntireCohort_pathvar$PCDH12 <- as.numeric(EntireCohort_pathvar$PCDH12)
unique(EntireCohort_pathvar$PCDH12)
# 0 1

table(EntireCohort_pathvar$PCDH12)
#   0      1 
# 462325   7474 




### 'OCLN'
EntireCohort_pathvar$OCLN <- NA
EntireCohort_pathvar$OCLN[EntireCohort_pathvar$GENE1 %in% c("OCLN")] <- 1
EntireCohort_pathvar$OCLN[EntireCohort_pathvar$GENE2 %in% c("OCLN")] <- 1
EntireCohort_pathvar$OCLN[EntireCohort_pathvar$GENE3 %in% c("OCLN")] <- 1
EntireCohort_pathvar$OCLN[EntireCohort_pathvar$GENE4 %in% c("OCLN")] <- 1
EntireCohort_pathvar$OCLN[EntireCohort_pathvar$GENE5 %in% c("OCLN")] <- 1
EntireCohort_pathvar$OCLN[EntireCohort_pathvar$GENE6 %in% c("OCLN")] <- 1
EntireCohort_pathvar$OCLN[EntireCohort_pathvar$GENE7 %in% c("OCLN")] <- 1
EntireCohort_pathvar$OCLN[EntireCohort_pathvar$GENE8 %in% c("OCLN")] <- 1
EntireCohort_pathvar$OCLN[EntireCohort_pathvar$GENE9 %in% c("OCLN")] <- 1

EntireCohort_pathvar$OCLN[is.na(EntireCohort_pathvar$OCLN)] <- 0
EntireCohort_pathvar$OCLN <- as.numeric(EntireCohort_pathvar$OCLN)
unique(EntireCohort_pathvar$OCLN)
# 0 1

table(EntireCohort_pathvar$OCLN)
#  0      1 
# 469169    630 





### 'FRMD4A'
EntireCohort_pathvar$FRMD4A <- NA
EntireCohort_pathvar$FRMD4A[EntireCohort_pathvar$GENE1 %in% c("FRMD4A")] <- 1
EntireCohort_pathvar$FRMD4A[EntireCohort_pathvar$GENE2 %in% c("FRMD4A")] <- 1
EntireCohort_pathvar$FRMD4A[EntireCohort_pathvar$GENE3 %in% c("FRMD4A")] <- 1
EntireCohort_pathvar$FRMD4A[EntireCohort_pathvar$GENE4 %in% c("FRMD4A")] <- 1
EntireCohort_pathvar$FRMD4A[EntireCohort_pathvar$GENE5 %in% c("FRMD4A")] <- 1
EntireCohort_pathvar$FRMD4A[EntireCohort_pathvar$GENE6 %in% c("FRMD4A")] <- 1
EntireCohort_pathvar$FRMD4A[EntireCohort_pathvar$GENE7 %in% c("FRMD4A")] <- 1
EntireCohort_pathvar$FRMD4A[EntireCohort_pathvar$GENE8 %in% c("FRMD4A")] <- 1
EntireCohort_pathvar$FRMD4A[EntireCohort_pathvar$GENE9 %in% c("FRMD4A")] <- 1

EntireCohort_pathvar$FRMD4A[is.na(EntireCohort_pathvar$FRMD4A)] <- 0
EntireCohort_pathvar$FRMD4A <- as.numeric(EntireCohort_pathvar$FRMD4A)
unique(EntireCohort_pathvar$FRMD4A)
# 0 1

table(EntireCohort_pathvar$FRMD4A)
#   0      1 
# 466197   3602




### 'NUP37' 
EntireCohort_pathvar$NUP37 <- NA
EntireCohort_pathvar$NUP37[EntireCohort_pathvar$GENE1 %in% c("NUP37")] <- 1
EntireCohort_pathvar$NUP37[EntireCohort_pathvar$GENE2 %in% c("NUP37")] <- 1
EntireCohort_pathvar$NUP37[EntireCohort_pathvar$GENE3 %in% c("NUP37")] <- 1
EntireCohort_pathvar$NUP37[EntireCohort_pathvar$GENE4 %in% c("NUP37")] <- 1
EntireCohort_pathvar$NUP37[EntireCohort_pathvar$GENE5 %in% c("NUP37")] <- 1
EntireCohort_pathvar$NUP37[EntireCohort_pathvar$GENE6 %in% c("NUP37")] <- 1
EntireCohort_pathvar$NUP37[EntireCohort_pathvar$GENE7 %in% c("NUP37")] <- 1
EntireCohort_pathvar$NUP37[EntireCohort_pathvar$GENE8 %in% c("NUP37")] <- 1
EntireCohort_pathvar$NUP37[EntireCohort_pathvar$GENE9 %in% c("NUP37")] <- 1

EntireCohort_pathvar$NUP37[is.na(EntireCohort_pathvar$NUP37)] <- 0
EntireCohort_pathvar$NUP37 <- as.numeric(EntireCohort_pathvar$NUP37)
unique(EntireCohort_pathvar$NUP37)
# 0 1

table(EntireCohort_pathvar$NUP37)
#   0      1 
# 469451    348 






# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)









# CCgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

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



# CentCytgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

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




# DDRgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

# go thr all the genes to check that they're present
grep("ERCC6", EntireCohort_pathvar$GENE)
# all seem to be present


# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$DDRgene <- NA
EntireCohort_pathvar$DDRgene[EntireCohort_pathvar$GENE1 %in% c('PNKP', 'RAD50', 'TONSL', 'TOP3A', 'RBBP8', 'TRAIP', 'NHEJ1', 'NSMCE2', 'ERCC8', 'ERCC6')] <- 1
EntireCohort_pathvar$DDRgene[EntireCohort_pathvar$GENE2 %in% c('PNKP', 'RAD50', 'TONSL', 'TOP3A', 'RBBP8', 'TRAIP', 'NHEJ1', 'NSMCE2', 'ERCC8', 'ERCC6')] <- 1
EntireCohort_pathvar$DDRgene[EntireCohort_pathvar$GENE3 %in% c('PNKP', 'RAD50', 'TONSL', 'TOP3A', 'RBBP8', 'TRAIP', 'NHEJ1', 'NSMCE2', 'ERCC8', 'ERCC6')] <- 1
EntireCohort_pathvar$DDRgene[EntireCohort_pathvar$GENE4 %in% c('PNKP', 'RAD50', 'TONSL', 'TOP3A', 'RBBP8', 'TRAIP', 'NHEJ1', 'NSMCE2', 'ERCC8', 'ERCC6')] <- 1
EntireCohort_pathvar$DDRgene[EntireCohort_pathvar$GENE5 %in% c('PNKP', 'RAD50', 'TONSL', 'TOP3A', 'RBBP8', 'TRAIP', 'NHEJ1', 'NSMCE2', 'ERCC8', 'ERCC6')] <- 1
EntireCohort_pathvar$DDRgene[EntireCohort_pathvar$GENE6 %in% c('PNKP', 'RAD50', 'TONSL', 'TOP3A', 'RBBP8', 'TRAIP', 'NHEJ1', 'NSMCE2', 'ERCC8', 'ERCC6')] <- 1
EntireCohort_pathvar$DDRgene[EntireCohort_pathvar$GENE7 %in% c('PNKP', 'RAD50', 'TONSL', 'TOP3A', 'RBBP8', 'TRAIP', 'NHEJ1', 'NSMCE2', 'ERCC8', 'ERCC6')] <- 1
EntireCohort_pathvar$DDRgene[EntireCohort_pathvar$GENE8 %in% c('PNKP', 'RAD50', 'TONSL', 'TOP3A', 'RBBP8', 'TRAIP', 'NHEJ1', 'NSMCE2', 'ERCC8', 'ERCC6')] <- 1
EntireCohort_pathvar$DDRgene[EntireCohort_pathvar$GENE9 %in% c('PNKP', 'RAD50', 'TONSL', 'TOP3A', 'RBBP8', 'TRAIP', 'NHEJ1', 'NSMCE2', 'ERCC8', 'ERCC6')] <- 1


EntireCohort_pathvar$DDRgene[is.na(EntireCohort_pathvar$DDRgene)] <- 0
EntireCohort_pathvar$DDRgene <- as.numeric(EntireCohort_pathvar$DDRgene)
unique(EntireCohort_pathvar$DDRgene)
# 0 1

table(EntireCohort_pathvar$DDRgene)
# 0      1 
# 438867  30932 




# now to make a column for each of the DDRgenes
### 'PNKP'
EntireCohort_pathvar$PNKP <- NA
EntireCohort_pathvar$PNKP[EntireCohort_pathvar$GENE1 %in% c("PNKP")] <- 1
EntireCohort_pathvar$PNKP[EntireCohort_pathvar$GENE2 %in% c("PNKP")] <- 1
EntireCohort_pathvar$PNKP[EntireCohort_pathvar$GENE3 %in% c("PNKP")] <- 1
EntireCohort_pathvar$PNKP[EntireCohort_pathvar$GENE4 %in% c("PNKP")] <- 1
EntireCohort_pathvar$PNKP[EntireCohort_pathvar$GENE5 %in% c("PNKP")] <- 1
EntireCohort_pathvar$PNKP[EntireCohort_pathvar$GENE6 %in% c("PNKP")] <- 1
EntireCohort_pathvar$PNKP[EntireCohort_pathvar$GENE7 %in% c("PNKP")] <- 1
EntireCohort_pathvar$PNKP[EntireCohort_pathvar$GENE8 %in% c("PNKP")] <- 1
EntireCohort_pathvar$PNKP[EntireCohort_pathvar$GENE9 %in% c("PNKP")] <- 1

EntireCohort_pathvar$PNKP[is.na(EntireCohort_pathvar$PNKP)] <- 0
EntireCohort_pathvar$PNKP <- as.numeric(EntireCohort_pathvar$PNKP)
unique(EntireCohort_pathvar$PNKP)
# 0 1

table(EntireCohort_pathvar$PNKP)
#  0      1 
# 467338   2461





### 'RAD50'
EntireCohort_pathvar$RAD50 <- NA
EntireCohort_pathvar$RAD50[EntireCohort_pathvar$GENE1 %in% c("RAD50")] <- 1
EntireCohort_pathvar$RAD50[EntireCohort_pathvar$GENE2 %in% c("RAD50")] <- 1
EntireCohort_pathvar$RAD50[EntireCohort_pathvar$GENE3 %in% c("RAD50")] <- 1
EntireCohort_pathvar$RAD50[EntireCohort_pathvar$GENE4 %in% c("RAD50")] <- 1
EntireCohort_pathvar$RAD50[EntireCohort_pathvar$GENE5 %in% c("RAD50")] <- 1
EntireCohort_pathvar$RAD50[EntireCohort_pathvar$GENE6 %in% c("RAD50")] <- 1
EntireCohort_pathvar$RAD50[EntireCohort_pathvar$GENE7 %in% c("RAD50")] <- 1
EntireCohort_pathvar$RAD50[EntireCohort_pathvar$GENE8 %in% c("RAD50")] <- 1
EntireCohort_pathvar$RAD50[EntireCohort_pathvar$GENE9 %in% c("RAD50")] <- 1

EntireCohort_pathvar$RAD50[is.na(EntireCohort_pathvar$RAD50)] <- 0
EntireCohort_pathvar$RAD50 <- as.numeric(EntireCohort_pathvar$RAD50)
unique(EntireCohort_pathvar$RAD50)
# 0 1

table(EntireCohort_pathvar$RAD50)
#  0      1 
# 467354   2445




### 'TONSL'
EntireCohort_pathvar$TONSL <- NA
EntireCohort_pathvar$TONSL[EntireCohort_pathvar$GENE1 %in% c("TONSL")] <- 1
EntireCohort_pathvar$TONSL[EntireCohort_pathvar$GENE2 %in% c("TONSL")] <- 1
EntireCohort_pathvar$TONSL[EntireCohort_pathvar$GENE3 %in% c("TONSL")] <- 1
EntireCohort_pathvar$TONSL[EntireCohort_pathvar$GENE4 %in% c("TONSL")] <- 1
EntireCohort_pathvar$TONSL[EntireCohort_pathvar$GENE5 %in% c("TONSL")] <- 1
EntireCohort_pathvar$TONSL[EntireCohort_pathvar$GENE6 %in% c("TONSL")] <- 1
EntireCohort_pathvar$TONSL[EntireCohort_pathvar$GENE7 %in% c("TONSL")] <- 1
EntireCohort_pathvar$TONSL[EntireCohort_pathvar$GENE8 %in% c("TONSL")] <- 1
EntireCohort_pathvar$TONSL[EntireCohort_pathvar$GENE9 %in% c("TONSL")] <- 1

EntireCohort_pathvar$TONSL[is.na(EntireCohort_pathvar$TONSL)] <- 0
EntireCohort_pathvar$TONSL <- as.numeric(EntireCohort_pathvar$TONSL)
unique(EntireCohort_pathvar$TONSL)
# 0 1

table(EntireCohort_pathvar$TONSL)
#   0      1 
# 468164   1635




### 'TOP3A'
EntireCohort_pathvar$TOP3A <- NA
EntireCohort_pathvar$TOP3A[EntireCohort_pathvar$GENE1 %in% c("TOP3A")] <- 1
EntireCohort_pathvar$TOP3A[EntireCohort_pathvar$GENE2 %in% c("TOP3A")] <- 1
EntireCohort_pathvar$TOP3A[EntireCohort_pathvar$GENE3 %in% c("TOP3A")] <- 1
EntireCohort_pathvar$TOP3A[EntireCohort_pathvar$GENE4 %in% c("TOP3A")] <- 1
EntireCohort_pathvar$TOP3A[EntireCohort_pathvar$GENE5 %in% c("TOP3A")] <- 1
EntireCohort_pathvar$TOP3A[EntireCohort_pathvar$GENE6 %in% c("TOP3A")] <- 1
EntireCohort_pathvar$TOP3A[EntireCohort_pathvar$GENE7 %in% c("TOP3A")] <- 1
EntireCohort_pathvar$TOP3A[EntireCohort_pathvar$GENE8 %in% c("TOP3A")] <- 1
EntireCohort_pathvar$TOP3A[EntireCohort_pathvar$GENE9 %in% c("TOP3A")] <- 1

EntireCohort_pathvar$TOP3A[is.na(EntireCohort_pathvar$TOP3A)] <- 0
EntireCohort_pathvar$TOP3A <- as.numeric(EntireCohort_pathvar$TOP3A)
unique(EntireCohort_pathvar$TOP3A)
# 0 1

table(EntireCohort_pathvar$TOP3A)
# 0      1 
# 461930   7869




### 'RBBP8'
EntireCohort_pathvar$RBBP8 <- NA
EntireCohort_pathvar$RBBP8[EntireCohort_pathvar$GENE1 %in% c("RBBP8")] <- 1
EntireCohort_pathvar$RBBP8[EntireCohort_pathvar$GENE2 %in% c("RBBP8")] <- 1
EntireCohort_pathvar$RBBP8[EntireCohort_pathvar$GENE3 %in% c("RBBP8")] <- 1
EntireCohort_pathvar$RBBP8[EntireCohort_pathvar$GENE4 %in% c("RBBP8")] <- 1
EntireCohort_pathvar$RBBP8[EntireCohort_pathvar$GENE5 %in% c("RBBP8")] <- 1
EntireCohort_pathvar$RBBP8[EntireCohort_pathvar$GENE6 %in% c("RBBP8")] <- 1
EntireCohort_pathvar$RBBP8[EntireCohort_pathvar$GENE7 %in% c("RBBP8")] <- 1
EntireCohort_pathvar$RBBP8[EntireCohort_pathvar$GENE8 %in% c("RBBP8")] <- 1
EntireCohort_pathvar$RBBP8[EntireCohort_pathvar$GENE9 %in% c("RBBP8")] <- 1

EntireCohort_pathvar$RBBP8[is.na(EntireCohort_pathvar$RBBP8)] <- 0
EntireCohort_pathvar$RBBP8 <- as.numeric(EntireCohort_pathvar$RBBP8)
unique(EntireCohort_pathvar$RBBP8)
# 0 1

table(EntireCohort_pathvar$RBBP8)
#  0      1 
# 468755   1044




### 'TRAIP'
EntireCohort_pathvar$TRAIP <- NA
EntireCohort_pathvar$TRAIP[EntireCohort_pathvar$GENE1 %in% c("TRAIP")] <- 1
EntireCohort_pathvar$TRAIP[EntireCohort_pathvar$GENE2 %in% c("TRAIP")] <- 1
EntireCohort_pathvar$TRAIP[EntireCohort_pathvar$GENE3 %in% c("TRAIP")] <- 1
EntireCohort_pathvar$TRAIP[EntireCohort_pathvar$GENE4 %in% c("TRAIP")] <- 1
EntireCohort_pathvar$TRAIP[EntireCohort_pathvar$GENE5 %in% c("TRAIP")] <- 1
EntireCohort_pathvar$TRAIP[EntireCohort_pathvar$GENE6 %in% c("TRAIP")] <- 1
EntireCohort_pathvar$TRAIP[EntireCohort_pathvar$GENE7 %in% c("TRAIP")] <- 1
EntireCohort_pathvar$TRAIP[EntireCohort_pathvar$GENE8 %in% c("TRAIP")] <- 1
EntireCohort_pathvar$TRAIP[EntireCohort_pathvar$GENE9 %in% c("TRAIP")] <- 1

EntireCohort_pathvar$TRAIP[is.na(EntireCohort_pathvar$TRAIP)] <- 0
EntireCohort_pathvar$TRAIP <- as.numeric(EntireCohort_pathvar$TRAIP)
unique(EntireCohort_pathvar$TRAIP)
# 0 1

table(EntireCohort_pathvar$TRAIP)
#  0      1 
# 469122    677 




### 'NHEJ1'
EntireCohort_pathvar$NHEJ1 <- NA
EntireCohort_pathvar$NHEJ1[EntireCohort_pathvar$GENE1 %in% c("NHEJ1")] <- 1
EntireCohort_pathvar$NHEJ1[EntireCohort_pathvar$GENE2 %in% c("NHEJ1")] <- 1
EntireCohort_pathvar$NHEJ1[EntireCohort_pathvar$GENE3 %in% c("NHEJ1")] <- 1
EntireCohort_pathvar$NHEJ1[EntireCohort_pathvar$GENE4 %in% c("NHEJ1")] <- 1
EntireCohort_pathvar$NHEJ1[EntireCohort_pathvar$GENE5 %in% c("NHEJ1")] <- 1
EntireCohort_pathvar$NHEJ1[EntireCohort_pathvar$GENE6 %in% c("NHEJ1")] <- 1
EntireCohort_pathvar$NHEJ1[EntireCohort_pathvar$GENE7 %in% c("NHEJ1")] <- 1
EntireCohort_pathvar$NHEJ1[EntireCohort_pathvar$GENE8 %in% c("NHEJ1")] <- 1
EntireCohort_pathvar$NHEJ1[EntireCohort_pathvar$GENE9 %in% c("NHEJ1")] <- 1

EntireCohort_pathvar$NHEJ1[is.na(EntireCohort_pathvar$NHEJ1)] <- 0
EntireCohort_pathvar$NHEJ1 <- as.numeric(EntireCohort_pathvar$NHEJ1)
unique(EntireCohort_pathvar$NHEJ1)
# 0 1

table(EntireCohort_pathvar$NHEJ1)
#  0      1 
# 468601   1198




### 'NSMCE2'
EntireCohort_pathvar$NSMCE2 <- NA
EntireCohort_pathvar$NSMCE2[EntireCohort_pathvar$GENE1 %in% c("NSMCE2")] <- 1
EntireCohort_pathvar$NSMCE2[EntireCohort_pathvar$GENE2 %in% c("NSMCE2")] <- 1
EntireCohort_pathvar$NSMCE2[EntireCohort_pathvar$GENE3 %in% c("NSMCE2")] <- 1
EntireCohort_pathvar$NSMCE2[EntireCohort_pathvar$GENE4 %in% c("NSMCE2")] <- 1
EntireCohort_pathvar$NSMCE2[EntireCohort_pathvar$GENE5 %in% c("NSMCE2")] <- 1
EntireCohort_pathvar$NSMCE2[EntireCohort_pathvar$GENE6 %in% c("NSMCE2")] <- 1
EntireCohort_pathvar$NSMCE2[EntireCohort_pathvar$GENE7 %in% c("NSMCE2")] <- 1
EntireCohort_pathvar$NSMCE2[EntireCohort_pathvar$GENE8 %in% c("NSMCE2")] <- 1
EntireCohort_pathvar$NSMCE2[EntireCohort_pathvar$GENE9 %in% c("NSMCE2")] <- 1

EntireCohort_pathvar$NSMCE2[is.na(EntireCohort_pathvar$NSMCE2)] <- 0
EntireCohort_pathvar$NSMCE2 <- as.numeric(EntireCohort_pathvar$NSMCE2)
unique(EntireCohort_pathvar$NSMCE2)
# 0 1

table(EntireCohort_pathvar$NSMCE2)
#  0      1 
# 469671    128




### 'ERCC8' 
EntireCohort_pathvar$ERCC8 <- NA
EntireCohort_pathvar$ERCC8[EntireCohort_pathvar$GENE1 %in% c("ERCC8")] <- 1
EntireCohort_pathvar$ERCC8[EntireCohort_pathvar$GENE2 %in% c("ERCC8")] <- 1
EntireCohort_pathvar$ERCC8[EntireCohort_pathvar$GENE3 %in% c("ERCC8")] <- 1
EntireCohort_pathvar$ERCC8[EntireCohort_pathvar$GENE4 %in% c("ERCC8")] <- 1
EntireCohort_pathvar$ERCC8[EntireCohort_pathvar$GENE5 %in% c("ERCC8")] <- 1
EntireCohort_pathvar$ERCC8[EntireCohort_pathvar$GENE6 %in% c("ERCC8")] <- 1
EntireCohort_pathvar$ERCC8[EntireCohort_pathvar$GENE7 %in% c("ERCC8")] <- 1
EntireCohort_pathvar$ERCC8[EntireCohort_pathvar$GENE8 %in% c("ERCC8")] <- 1
EntireCohort_pathvar$ERCC8[EntireCohort_pathvar$GENE9 %in% c("ERCC8")] <- 1

EntireCohort_pathvar$ERCC8[is.na(EntireCohort_pathvar$ERCC8)] <- 0
EntireCohort_pathvar$ERCC8 <- as.numeric(EntireCohort_pathvar$ERCC8)
unique(EntireCohort_pathvar$ERCC8)
# 0 1

table(EntireCohort_pathvar$ERCC8)
#  0      1 
# 467164   2635





### 'ERCC6'     
EntireCohort_pathvar$ERCC6 <- NA
EntireCohort_pathvar$ERCC6[EntireCohort_pathvar$GENE1 %in% c("ERCC6")] <- 1
EntireCohort_pathvar$ERCC6[EntireCohort_pathvar$GENE2 %in% c("ERCC6")] <- 1
EntireCohort_pathvar$ERCC6[EntireCohort_pathvar$GENE3 %in% c("ERCC6")] <- 1
EntireCohort_pathvar$ERCC6[EntireCohort_pathvar$GENE4 %in% c("ERCC6")] <- 1
EntireCohort_pathvar$ERCC6[EntireCohort_pathvar$GENE5 %in% c("ERCC6")] <- 1
EntireCohort_pathvar$ERCC6[EntireCohort_pathvar$GENE6 %in% c("ERCC6")] <- 1
EntireCohort_pathvar$ERCC6[EntireCohort_pathvar$GENE7 %in% c("ERCC6")] <- 1
EntireCohort_pathvar$ERCC6[EntireCohort_pathvar$GENE8 %in% c("ERCC6")] <- 1
EntireCohort_pathvar$ERCC6[EntireCohort_pathvar$GENE9 %in% c("ERCC6")] <- 1

EntireCohort_pathvar$ERCC6[is.na(EntireCohort_pathvar$ERCC6)] <- 0
EntireCohort_pathvar$ERCC6 <- as.numeric(EntireCohort_pathvar$ERCC6)
unique(EntireCohort_pathvar$ERCC6)
# 0 1

table(EntireCohort_pathvar$ERCC6)
#  0      1 
# 458125  11674







# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)






# GHgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

# go thr all the genes to check that they're present
grep("STAT5B", EntireCohort_pathvar$GENE)
# all seem to be present

# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$GHgene <- NA
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE1 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE2 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE3 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE4 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE5 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE6 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE7 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE8 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE9 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1


EntireCohort_pathvar$GHgene[is.na(EntireCohort_pathvar$GHgene)] <- 0
EntireCohort_pathvar$GHgene <- as.numeric(EntireCohort_pathvar$GHgene)
unique(EntireCohort_pathvar$GHgene)
# 0 1

table(EntireCohort_pathvar$GHgene)
#      0      1 
# 456368  13431



# now to make a column for each of the GHgenes
### 'BTK'
EntireCohort_pathvar$BTK <- NA
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE1 %in% c("BTK")] <- 1
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE2 %in% c("BTK")] <- 1
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE3 %in% c("BTK")] <- 1
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE4 %in% c("BTK")] <- 1
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE5 %in% c("BTK")] <- 1
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE6 %in% c("BTK")] <- 1
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE7 %in% c("BTK")] <- 1
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE8 %in% c("BTK")] <- 1
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE9 %in% c("BTK")] <- 1

EntireCohort_pathvar$BTK[is.na(EntireCohort_pathvar$BTK)] <- 0
EntireCohort_pathvar$BTK <- as.numeric(EntireCohort_pathvar$BTK)
unique(EntireCohort_pathvar$BTK)
# 0 1

table(EntireCohort_pathvar$BTK)
# 0      1 
# 469710     89





### 'GH1'
EntireCohort_pathvar$GH1 <- NA
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE1 %in% c("GH1")] <- 1
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE2 %in% c("GH1")] <- 1
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE3 %in% c("GH1")] <- 1
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE4 %in% c("GH1")] <- 1
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE5 %in% c("GH1")] <- 1
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE6 %in% c("GH1")] <- 1
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE7 %in% c("GH1")] <- 1
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE8 %in% c("GH1")] <- 1
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE9 %in% c("GH1")] <- 1

EntireCohort_pathvar$GH1[is.na(EntireCohort_pathvar$GH1)] <- 0
EntireCohort_pathvar$GH1 <- as.numeric(EntireCohort_pathvar$GH1)
unique(EntireCohort_pathvar$GH1)
# 0 1

table(EntireCohort_pathvar$GH1)
# 0      1 
# 469552    247




### 'GHR'
EntireCohort_pathvar$GHR <- NA
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE1 %in% c("GHR")] <- 1
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE2 %in% c("GHR")] <- 1
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE3 %in% c("GHR")] <- 1
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE4 %in% c("GHR")] <- 1
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE5 %in% c("GHR")] <- 1
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE6 %in% c("GHR")] <- 1
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE7 %in% c("GHR")] <- 1
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE8 %in% c("GHR")] <- 1
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE9 %in% c("GHR")] <- 1

EntireCohort_pathvar$GHR[is.na(EntireCohort_pathvar$GHR)] <- 0
EntireCohort_pathvar$GHR <- as.numeric(EntireCohort_pathvar$GHR)
unique(EntireCohort_pathvar$GHR)
# 0 1

table(EntireCohort_pathvar$GHR)
# 0      1 
# 467645   2154




### 'GHRH'
EntireCohort_pathvar$GHRH <- NA
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE1 %in% c("GHRH")] <- 1
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE2 %in% c("GHRH")] <- 1
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE3 %in% c("GHRH")] <- 1
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE4 %in% c("GHRH")] <- 1
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE5 %in% c("GHRH")] <- 1
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE6 %in% c("GHRH")] <- 1
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE7 %in% c("GHRH")] <- 1
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE8 %in% c("GHRH")] <- 1
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE9 %in% c("GHRH")] <- 1

EntireCohort_pathvar$GHRH[is.na(EntireCohort_pathvar$GHRH)] <- 0
EntireCohort_pathvar$GHRH <- as.numeric(EntireCohort_pathvar$GHRH)
unique(EntireCohort_pathvar$GHRH)
# 0 1

table(EntireCohort_pathvar$GHRH)
# 0      1 
# 469712     87




### 'GHRHR'
EntireCohort_pathvar$GHRHR <- NA
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE1 %in% c("GHRHR")] <- 1
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE2 %in% c("GHRHR")] <- 1
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE3 %in% c("GHRHR")] <- 1
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE4 %in% c("GHRHR")] <- 1
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE5 %in% c("GHRHR")] <- 1
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE6 %in% c("GHRHR")] <- 1
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE7 %in% c("GHRHR")] <- 1
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE8 %in% c("GHRHR")] <- 1
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE9 %in% c("GHRHR")] <- 1

EntireCohort_pathvar$GHRHR[is.na(EntireCohort_pathvar$GHRHR)] <- 0
EntireCohort_pathvar$GHRHR <- as.numeric(EntireCohort_pathvar$GHRHR)
unique(EntireCohort_pathvar$GHRHR)
# 0 1

table(EntireCohort_pathvar$GHRHR)
# 0      1 
# 469138    661




### 'GHSR'
EntireCohort_pathvar$GHSR <- NA
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE1 %in% c("GHSR")] <- 1
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE2 %in% c("GHSR")] <- 1
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE3 %in% c("GHSR")] <- 1
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE4 %in% c("GHSR")] <- 1
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE5 %in% c("GHSR")] <- 1
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE6 %in% c("GHSR")] <- 1
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE7 %in% c("GHSR")] <- 1
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE8 %in% c("GHSR")] <- 1
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE9 %in% c("GHSR")] <- 1

EntireCohort_pathvar$GHSR[is.na(EntireCohort_pathvar$GHSR)] <- 0
EntireCohort_pathvar$GHSR <- as.numeric(EntireCohort_pathvar$GHSR)
unique(EntireCohort_pathvar$GHSR)
# 0 1

table(EntireCohort_pathvar$GHSR)
# 0      1 
# 468753   1046




### 'HESX1'
EntireCohort_pathvar$HESX1 <- NA
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE1 %in% c("HESX1")] <- 1
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE2 %in% c("HESX1")] <- 1
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE3 %in% c("HESX1")] <- 1
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE4 %in% c("HESX1")] <- 1
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE5 %in% c("HESX1")] <- 1
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE6 %in% c("HESX1")] <- 1
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE7 %in% c("HESX1")] <- 1
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE8 %in% c("HESX1")] <- 1
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE9 %in% c("HESX1")] <- 1

EntireCohort_pathvar$HESX1[is.na(EntireCohort_pathvar$HESX1)] <- 0
EntireCohort_pathvar$HESX1 <- as.numeric(EntireCohort_pathvar$HESX1)
unique(EntireCohort_pathvar$HESX1)
# 0 1

table(EntireCohort_pathvar$HESX1)
# 0      1 
# 469289    510




### 'IGF1'
EntireCohort_pathvar$IGF1 <- NA
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE1 %in% c("IGF1")] <- 1
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE2 %in% c("IGF1")] <- 1
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE3 %in% c("IGF1")] <- 1
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE4 %in% c("IGF1")] <- 1
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE5 %in% c("IGF1")] <- 1
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE6 %in% c("IGF1")] <- 1
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE7 %in% c("IGF1")] <- 1
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE8 %in% c("IGF1")] <- 1
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE9 %in% c("IGF1")] <- 1

EntireCohort_pathvar$IGF1[is.na(EntireCohort_pathvar$IGF1)] <- 0
EntireCohort_pathvar$IGF1 <- as.numeric(EntireCohort_pathvar$IGF1)
unique(EntireCohort_pathvar$IGF1)
# 0 1

table(EntireCohort_pathvar$IGF1)
# 0      1 
# 469390    409




### 'IGF1R'
EntireCohort_pathvar$IGF1R <- NA
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE1 %in% c("IGF1R")] <- 1
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE2 %in% c("IGF1R")] <- 1
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE3 %in% c("IGF1R")] <- 1
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE4 %in% c("IGF1R")] <- 1
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE5 %in% c("IGF1R")] <- 1
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE6 %in% c("IGF1R")] <- 1
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE7 %in% c("IGF1R")] <- 1
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE8 %in% c("IGF1R")] <- 1
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE9 %in% c("IGF1R")] <- 1

EntireCohort_pathvar$IGF1R[is.na(EntireCohort_pathvar$IGF1R)] <- 0
EntireCohort_pathvar$IGF1R <- as.numeric(EntireCohort_pathvar$IGF1R)
unique(EntireCohort_pathvar$IGF1R)
# 0 1

table(EntireCohort_pathvar$IGF1R)
# 0      1 
# 467800   1999




### 'IGFALS'
EntireCohort_pathvar$IGFALS <- NA
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE1 %in% c("IGFALS")] <- 1
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE2 %in% c("IGFALS")] <- 1
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE3 %in% c("IGFALS")] <- 1
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE4 %in% c("IGFALS")] <- 1
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE5 %in% c("IGFALS")] <- 1
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE6 %in% c("IGFALS")] <- 1
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE7 %in% c("IGFALS")] <- 1
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE8 %in% c("IGFALS")] <- 1
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE9 %in% c("IGFALS")] <- 1

EntireCohort_pathvar$IGFALS[is.na(EntireCohort_pathvar$IGFALS)] <- 0
EntireCohort_pathvar$IGFALS <- as.numeric(EntireCohort_pathvar$IGFALS)
unique(EntireCohort_pathvar$IGFALS)
# 0 1

table(EntireCohort_pathvar$IGFALS)
# 0      1 
# 468516   1283




### 'LHX3'
EntireCohort_pathvar$LHX3 <- NA
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE1 %in% c("LHX3")] <- 1
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE2 %in% c("LHX3")] <- 1
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE3 %in% c("LHX3")] <- 1
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE4 %in% c("LHX3")] <- 1
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE5 %in% c("LHX3")] <- 1
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE6 %in% c("LHX3")] <- 1
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE7 %in% c("LHX3")] <- 1
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE8 %in% c("LHX3")] <- 1
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE9 %in% c("LHX3")] <- 1

EntireCohort_pathvar$LHX3[is.na(EntireCohort_pathvar$LHX3)] <- 0
EntireCohort_pathvar$LHX3 <- as.numeric(EntireCohort_pathvar$LHX3)
unique(EntireCohort_pathvar$LHX3)
# 0 1

table(EntireCohort_pathvar$LHX3)
# 0      1 
# 468886    913




### 'LHX4'
EntireCohort_pathvar$LHX4 <- NA
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE1 %in% c("LHX4")] <- 1
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE2 %in% c("LHX4")] <- 1
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE3 %in% c("LHX4")] <- 1
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE4 %in% c("LHX4")] <- 1
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE5 %in% c("LHX4")] <- 1
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE6 %in% c("LHX4")] <- 1
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE7 %in% c("LHX4")] <- 1
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE8 %in% c("LHX4")] <- 1
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE9 %in% c("LHX4")] <- 1

EntireCohort_pathvar$LHX4[is.na(EntireCohort_pathvar$LHX4)] <- 0
EntireCohort_pathvar$LHX4 <- as.numeric(EntireCohort_pathvar$LHX4)
unique(EntireCohort_pathvar$LHX4)
# 0 1

table(EntireCohort_pathvar$LHX4)
# 0      1 
# 469232    567




### 'OTX2'
EntireCohort_pathvar$OTX2 <- NA
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE1 %in% c("OTX2")] <- 1
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE2 %in% c("OTX2")] <- 1
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE3 %in% c("OTX2")] <- 1
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE4 %in% c("OTX2")] <- 1
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE5 %in% c("OTX2")] <- 1
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE6 %in% c("OTX2")] <- 1
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE7 %in% c("OTX2")] <- 1
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE8 %in% c("OTX2")] <- 1
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE9 %in% c("OTX2")] <- 1

EntireCohort_pathvar$OTX2[is.na(EntireCohort_pathvar$OTX2)] <- 0
EntireCohort_pathvar$OTX2 <- as.numeric(EntireCohort_pathvar$OTX2)
unique(EntireCohort_pathvar$OTX2)
# 0 1

table(EntireCohort_pathvar$OTX2)
# 0      1 
# 469074    725




### 'PNPLA6'
EntireCohort_pathvar$PNPLA6 <- NA
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE1 %in% c("PNPLA6")] <- 1
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE2 %in% c("PNPLA6")] <- 1
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE3 %in% c("PNPLA6")] <- 1
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE4 %in% c("PNPLA6")] <- 1
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE5 %in% c("PNPLA6")] <- 1
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE6 %in% c("PNPLA6")] <- 1
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE7 %in% c("PNPLA6")] <- 1
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE8 %in% c("PNPLA6")] <- 1
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE9 %in% c("PNPLA6")] <- 1

EntireCohort_pathvar$PNPLA6[is.na(EntireCohort_pathvar$PNPLA6)] <- 0
EntireCohort_pathvar$PNPLA6 <- as.numeric(EntireCohort_pathvar$PNPLA6)
unique(EntireCohort_pathvar$PNPLA6)
# 0 1

table(EntireCohort_pathvar$PNPLA6)
# 0      1 
# 468124   1675




### 'POU1F1'
EntireCohort_pathvar$POU1F1 <- NA
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE1 %in% c("POU1F1")] <- 1
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE2 %in% c("POU1F1")] <- 1
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE3 %in% c("POU1F1")] <- 1
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE4 %in% c("POU1F1")] <- 1
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE5 %in% c("POU1F1")] <- 1
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE6 %in% c("POU1F1")] <- 1
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE7 %in% c("POU1F1")] <- 1
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE8 %in% c("POU1F1")] <- 1
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE9 %in% c("POU1F1")] <- 1

EntireCohort_pathvar$POU1F1[is.na(EntireCohort_pathvar$POU1F1)] <- 0
EntireCohort_pathvar$POU1F1 <- as.numeric(EntireCohort_pathvar$POU1F1)
unique(EntireCohort_pathvar$POU1F1)
# 0 1

table(EntireCohort_pathvar$POU1F1)
# 0      1 
# 469210    589




### 'PROP1'
EntireCohort_pathvar$PROP1 <- NA
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE1 %in% c("PROP1")] <- 1
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE2 %in% c("PROP1")] <- 1
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE3 %in% c("PROP1")] <- 1
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE4 %in% c("PROP1")] <- 1
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE5 %in% c("PROP1")] <- 1
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE6 %in% c("PROP1")] <- 1
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE7 %in% c("PROP1")] <- 1
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE8 %in% c("PROP1")] <- 1
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE9 %in% c("PROP1")] <- 1

EntireCohort_pathvar$PROP1[is.na(EntireCohort_pathvar$PROP1)] <- 0
EntireCohort_pathvar$PROP1 <- as.numeric(EntireCohort_pathvar$PROP1)
unique(EntireCohort_pathvar$PROP1)
# 0 1

table(EntireCohort_pathvar$PROP1)
# 0      1 
# 469703     96




### 'STAT5B'
EntireCohort_pathvar$STAT5B <- NA
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE1 %in% c("STAT5B")] <- 1
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE2 %in% c("STAT5B")] <- 1
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE3 %in% c("STAT5B")] <- 1
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE4 %in% c("STAT5B")] <- 1
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE5 %in% c("STAT5B")] <- 1
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE6 %in% c("STAT5B")] <- 1
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE7 %in% c("STAT5B")] <- 1
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE8 %in% c("STAT5B")] <- 1
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE9 %in% c("STAT5B")] <- 1

EntireCohort_pathvar$STAT5B[is.na(EntireCohort_pathvar$STAT5B)] <- 0
EntireCohort_pathvar$STAT5B <- as.numeric(EntireCohort_pathvar$STAT5B)
unique(EntireCohort_pathvar$STAT5B)
# 0 1

table(EntireCohort_pathvar$STAT5B)
# 0      1 
# 469249    550





# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)














# LTUgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

# go thr all the genes to check that they're present
grep("STAMBP", EntireCohort_pathvar$GENE)
# all seem to be present
# SPOP = 0
# RNF13 = 0



# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$LTUgene <- NA
EntireCohort_pathvar$LTUgene[EntireCohort_pathvar$GENE1 %in% c('GORAB', 'IER3IP1', 'COMP', 'UFSP2', 'OBSL1', 'DYM', 'SBF1', 'EXOC6B', 'RAB33B', 'TRAPPC6B', 'ARCN1', 'COPB2', 'CTSA', 'GALNS', 'GLB1', 'GNPTAB', 'GUSB', 'ARSB', 'CUL7', 'CCDC8', 'TBC1D23', 'UBR1', 'STAMBP')] <- 1
EntireCohort_pathvar$LTUgene[EntireCohort_pathvar$GENE2 %in% c('GORAB', 'IER3IP1', 'COMP', 'UFSP2', 'OBSL1', 'DYM', 'SBF1', 'EXOC6B', 'RAB33B', 'TRAPPC6B', 'ARCN1', 'COPB2', 'CTSA', 'GALNS', 'GLB1', 'GNPTAB', 'GUSB', 'ARSB', 'CUL7', 'CCDC8', 'TBC1D23', 'UBR1', 'STAMBP')] <- 1
EntireCohort_pathvar$LTUgene[EntireCohort_pathvar$GENE3 %in% c('GORAB', 'IER3IP1', 'COMP', 'UFSP2', 'OBSL1', 'DYM', 'SBF1', 'EXOC6B', 'RAB33B', 'TRAPPC6B', 'ARCN1', 'COPB2', 'CTSA', 'GALNS', 'GLB1', 'GNPTAB', 'GUSB', 'ARSB', 'CUL7', 'CCDC8', 'TBC1D23', 'UBR1', 'STAMBP')] <- 1
EntireCohort_pathvar$LTUgene[EntireCohort_pathvar$GENE4 %in% c('GORAB', 'IER3IP1', 'COMP', 'UFSP2', 'OBSL1', 'DYM', 'SBF1', 'EXOC6B', 'RAB33B', 'TRAPPC6B', 'ARCN1', 'COPB2', 'CTSA', 'GALNS', 'GLB1', 'GNPTAB', 'GUSB', 'ARSB', 'CUL7', 'CCDC8', 'TBC1D23', 'UBR1', 'STAMBP')] <- 1
EntireCohort_pathvar$LTUgene[EntireCohort_pathvar$GENE5 %in% c('GORAB', 'IER3IP1', 'COMP', 'UFSP2', 'OBSL1', 'DYM', 'SBF1', 'EXOC6B', 'RAB33B', 'TRAPPC6B', 'ARCN1', 'COPB2', 'CTSA', 'GALNS', 'GLB1', 'GNPTAB', 'GUSB', 'ARSB', 'CUL7', 'CCDC8', 'TBC1D23', 'UBR1', 'STAMBP')] <- 1
EntireCohort_pathvar$LTUgene[EntireCohort_pathvar$GENE6 %in% c('GORAB', 'IER3IP1', 'COMP', 'UFSP2', 'OBSL1', 'DYM', 'SBF1', 'EXOC6B', 'RAB33B', 'TRAPPC6B', 'ARCN1', 'COPB2', 'CTSA', 'GALNS', 'GLB1', 'GNPTAB', 'GUSB', 'ARSB', 'CUL7', 'CCDC8', 'TBC1D23', 'UBR1', 'STAMBP')] <- 1
EntireCohort_pathvar$LTUgene[EntireCohort_pathvar$GENE7 %in% c('GORAB', 'IER3IP1', 'COMP', 'UFSP2', 'OBSL1', 'DYM', 'SBF1', 'EXOC6B', 'RAB33B', 'TRAPPC6B', 'ARCN1', 'COPB2', 'CTSA', 'GALNS', 'GLB1', 'GNPTAB', 'GUSB', 'ARSB', 'CUL7', 'CCDC8', 'TBC1D23', 'UBR1', 'STAMBP')] <- 1
EntireCohort_pathvar$LTUgene[EntireCohort_pathvar$GENE8 %in% c('GORAB', 'IER3IP1', 'COMP', 'UFSP2', 'OBSL1', 'DYM', 'SBF1', 'EXOC6B', 'RAB33B', 'TRAPPC6B', 'ARCN1', 'COPB2', 'CTSA', 'GALNS', 'GLB1', 'GNPTAB', 'GUSB', 'ARSB', 'CUL7', 'CCDC8', 'TBC1D23', 'UBR1', 'STAMBP')] <- 1
EntireCohort_pathvar$LTUgene[EntireCohort_pathvar$GENE9 %in% c('GORAB', 'IER3IP1', 'COMP', 'UFSP2', 'OBSL1', 'DYM', 'SBF1', 'EXOC6B', 'RAB33B', 'TRAPPC6B', 'ARCN1', 'COPB2', 'CTSA', 'GALNS', 'GLB1', 'GNPTAB', 'GUSB', 'ARSB', 'CUL7', 'CCDC8', 'TBC1D23', 'UBR1', 'STAMBP')] <- 1


EntireCohort_pathvar$LTUgene[is.na(EntireCohort_pathvar$LTUgene)] <- 0
EntireCohort_pathvar$LTUgene <- as.numeric(EntireCohort_pathvar$LTUgene)
unique(EntireCohort_pathvar$LTUgene)
# 0 1

table(EntireCohort_pathvar$LTUgene)
# 0      1 
# 429354  40445




# now to make a column for each of the LTUgenes
### 'GORAB'
EntireCohort_pathvar$GORAB <- NA
EntireCohort_pathvar$GORAB[EntireCohort_pathvar$GENE1 %in% c("GORAB")] <- 1
EntireCohort_pathvar$GORAB[EntireCohort_pathvar$GENE2 %in% c("GORAB")] <- 1
EntireCohort_pathvar$GORAB[EntireCohort_pathvar$GENE3 %in% c("GORAB")] <- 1
EntireCohort_pathvar$GORAB[EntireCohort_pathvar$GENE4 %in% c("GORAB")] <- 1
EntireCohort_pathvar$GORAB[EntireCohort_pathvar$GENE5 %in% c("GORAB")] <- 1
EntireCohort_pathvar$GORAB[EntireCohort_pathvar$GENE6 %in% c("GORAB")] <- 1
EntireCohort_pathvar$GORAB[EntireCohort_pathvar$GENE7 %in% c("GORAB")] <- 1
EntireCohort_pathvar$GORAB[EntireCohort_pathvar$GENE8 %in% c("GORAB")] <- 1
EntireCohort_pathvar$GORAB[EntireCohort_pathvar$GENE9 %in% c("GORAB")] <- 1

EntireCohort_pathvar$GORAB[is.na(EntireCohort_pathvar$GORAB)] <- 0
EntireCohort_pathvar$GORAB <- as.numeric(EntireCohort_pathvar$GORAB)
unique(EntireCohort_pathvar$GORAB)
# 0 1

table(EntireCohort_pathvar$GORAB)
# 0      1 
# 469310    489





### 'IER3IP1'
EntireCohort_pathvar$IER3IP1 <- NA
EntireCohort_pathvar$IER3IP1[EntireCohort_pathvar$GENE1 %in% c("IER3IP1")] <- 1
EntireCohort_pathvar$IER3IP1[EntireCohort_pathvar$GENE2 %in% c("IER3IP1")] <- 1
EntireCohort_pathvar$IER3IP1[EntireCohort_pathvar$GENE3 %in% c("IER3IP1")] <- 1
EntireCohort_pathvar$IER3IP1[EntireCohort_pathvar$GENE4 %in% c("IER3IP1")] <- 1
EntireCohort_pathvar$IER3IP1[EntireCohort_pathvar$GENE5 %in% c("IER3IP1")] <- 1
EntireCohort_pathvar$IER3IP1[EntireCohort_pathvar$GENE6 %in% c("IER3IP1")] <- 1
EntireCohort_pathvar$IER3IP1[EntireCohort_pathvar$GENE7 %in% c("IER3IP1")] <- 1
EntireCohort_pathvar$IER3IP1[EntireCohort_pathvar$GENE8 %in% c("IER3IP1")] <- 1
EntireCohort_pathvar$IER3IP1[EntireCohort_pathvar$GENE9 %in% c("IER3IP1")] <- 1

EntireCohort_pathvar$IER3IP1[is.na(EntireCohort_pathvar$IER3IP1)] <- 0
EntireCohort_pathvar$IER3IP1 <- as.numeric(EntireCohort_pathvar$IER3IP1)
unique(EntireCohort_pathvar$IER3IP1)
# 0 1

table(EntireCohort_pathvar$IER3IP1)
#   0      1 
# 469611    188




### 'COMP'
EntireCohort_pathvar$COMP <- NA
EntireCohort_pathvar$COMP[EntireCohort_pathvar$GENE1 %in% c("COMP")] <- 1
EntireCohort_pathvar$COMP[EntireCohort_pathvar$GENE2 %in% c("COMP")] <- 1
EntireCohort_pathvar$COMP[EntireCohort_pathvar$GENE3 %in% c("COMP")] <- 1
EntireCohort_pathvar$COMP[EntireCohort_pathvar$GENE4 %in% c("COMP")] <- 1
EntireCohort_pathvar$COMP[EntireCohort_pathvar$GENE5 %in% c("COMP")] <- 1
EntireCohort_pathvar$COMP[EntireCohort_pathvar$GENE6 %in% c("COMP")] <- 1
EntireCohort_pathvar$COMP[EntireCohort_pathvar$GENE7 %in% c("COMP")] <- 1
EntireCohort_pathvar$COMP[EntireCohort_pathvar$GENE8 %in% c("COMP")] <- 1
EntireCohort_pathvar$COMP[EntireCohort_pathvar$GENE9 %in% c("COMP")] <- 1

EntireCohort_pathvar$COMP[is.na(EntireCohort_pathvar$COMP)] <- 0
EntireCohort_pathvar$COMP <- as.numeric(EntireCohort_pathvar$COMP)
unique(EntireCohort_pathvar$COMP)
# 0 1

table(EntireCohort_pathvar$COMP)
#  0      1 
# 468478   1321 




### 'UFSP2'
EntireCohort_pathvar$UFSP2 <- NA
EntireCohort_pathvar$UFSP2[EntireCohort_pathvar$GENE1 %in% c("UFSP2")] <- 1
EntireCohort_pathvar$UFSP2[EntireCohort_pathvar$GENE2 %in% c("UFSP2")] <- 1
EntireCohort_pathvar$UFSP2[EntireCohort_pathvar$GENE3 %in% c("UFSP2")] <- 1
EntireCohort_pathvar$UFSP2[EntireCohort_pathvar$GENE4 %in% c("UFSP2")] <- 1
EntireCohort_pathvar$UFSP2[EntireCohort_pathvar$GENE5 %in% c("UFSP2")] <- 1
EntireCohort_pathvar$UFSP2[EntireCohort_pathvar$GENE6 %in% c("UFSP2")] <- 1
EntireCohort_pathvar$UFSP2[EntireCohort_pathvar$GENE7 %in% c("UFSP2")] <- 1
EntireCohort_pathvar$UFSP2[EntireCohort_pathvar$GENE8 %in% c("UFSP2")] <- 1
EntireCohort_pathvar$UFSP2[EntireCohort_pathvar$GENE9 %in% c("UFSP2")] <- 1

EntireCohort_pathvar$UFSP2[is.na(EntireCohort_pathvar$UFSP2)] <- 0
EntireCohort_pathvar$UFSP2 <- as.numeric(EntireCohort_pathvar$UFSP2)
unique(EntireCohort_pathvar$UFSP2)
# 0 1

table(EntireCohort_pathvar$UFSP2)
#  0      1 
# 468992    807




### 'SPOP'
# no variants present




### 'OBSL1'  
EntireCohort_pathvar$OBSL1 <- NA
EntireCohort_pathvar$OBSL1[EntireCohort_pathvar$GENE1 %in% c("OBSL1")] <- 1
EntireCohort_pathvar$OBSL1[EntireCohort_pathvar$GENE2 %in% c("OBSL1")] <- 1
EntireCohort_pathvar$OBSL1[EntireCohort_pathvar$GENE3 %in% c("OBSL1")] <- 1
EntireCohort_pathvar$OBSL1[EntireCohort_pathvar$GENE4 %in% c("OBSL1")] <- 1
EntireCohort_pathvar$OBSL1[EntireCohort_pathvar$GENE5 %in% c("OBSL1")] <- 1
EntireCohort_pathvar$OBSL1[EntireCohort_pathvar$GENE6 %in% c("OBSL1")] <- 1
EntireCohort_pathvar$OBSL1[EntireCohort_pathvar$GENE7 %in% c("OBSL1")] <- 1
EntireCohort_pathvar$OBSL1[EntireCohort_pathvar$GENE8 %in% c("OBSL1")] <- 1
EntireCohort_pathvar$OBSL1[EntireCohort_pathvar$GENE9 %in% c("OBSL1")] <- 1

EntireCohort_pathvar$OBSL1[is.na(EntireCohort_pathvar$OBSL1)] <- 0
EntireCohort_pathvar$OBSL1 <- as.numeric(EntireCohort_pathvar$OBSL1)
unique(EntireCohort_pathvar$OBSL1)
# 0 1

table(EntireCohort_pathvar$OBSL1)
#  0      1 
# 464267   5532 




### 'DYM'
EntireCohort_pathvar$DYM <- NA
EntireCohort_pathvar$DYM[EntireCohort_pathvar$GENE1 %in% c("DYM")] <- 1
EntireCohort_pathvar$DYM[EntireCohort_pathvar$GENE2 %in% c("DYM")] <- 1
EntireCohort_pathvar$DYM[EntireCohort_pathvar$GENE3 %in% c("DYM")] <- 1
EntireCohort_pathvar$DYM[EntireCohort_pathvar$GENE4 %in% c("DYM")] <- 1
EntireCohort_pathvar$DYM[EntireCohort_pathvar$GENE5 %in% c("DYM")] <- 1
EntireCohort_pathvar$DYM[EntireCohort_pathvar$GENE6 %in% c("DYM")] <- 1
EntireCohort_pathvar$DYM[EntireCohort_pathvar$GENE7 %in% c("DYM")] <- 1
EntireCohort_pathvar$DYM[EntireCohort_pathvar$GENE8 %in% c("DYM")] <- 1
EntireCohort_pathvar$DYM[EntireCohort_pathvar$GENE9 %in% c("DYM")] <- 1

EntireCohort_pathvar$DYM[is.na(EntireCohort_pathvar$DYM)] <- 0
EntireCohort_pathvar$DYM <- as.numeric(EntireCohort_pathvar$DYM)
unique(EntireCohort_pathvar$DYM)
# 0 1

table(EntireCohort_pathvar$DYM)
#  0      1 
# 467951   1848 




### 'SBF1'
EntireCohort_pathvar$SBF1 <- NA
EntireCohort_pathvar$SBF1[EntireCohort_pathvar$GENE1 %in% c("SBF1")] <- 1
EntireCohort_pathvar$SBF1[EntireCohort_pathvar$GENE2 %in% c("SBF1")] <- 1
EntireCohort_pathvar$SBF1[EntireCohort_pathvar$GENE3 %in% c("SBF1")] <- 1
EntireCohort_pathvar$SBF1[EntireCohort_pathvar$GENE4 %in% c("SBF1")] <- 1
EntireCohort_pathvar$SBF1[EntireCohort_pathvar$GENE5 %in% c("SBF1")] <- 1
EntireCohort_pathvar$SBF1[EntireCohort_pathvar$GENE6 %in% c("SBF1")] <- 1
EntireCohort_pathvar$SBF1[EntireCohort_pathvar$GENE7 %in% c("SBF1")] <- 1
EntireCohort_pathvar$SBF1[EntireCohort_pathvar$GENE8 %in% c("SBF1")] <- 1
EntireCohort_pathvar$SBF1[EntireCohort_pathvar$GENE9 %in% c("SBF1")] <- 1

EntireCohort_pathvar$SBF1[is.na(EntireCohort_pathvar$SBF1)] <- 0
EntireCohort_pathvar$SBF1 <- as.numeric(EntireCohort_pathvar$SBF1)
unique(EntireCohort_pathvar$SBF1)
# 0 1

table(EntireCohort_pathvar$SBF1)
#  0      1 
# 464712   5087 




### 'RNF13'
# no variants present





### 'EXOC6B'   
EntireCohort_pathvar$EXOC6B <- NA
EntireCohort_pathvar$EXOC6B[EntireCohort_pathvar$GENE1 %in% c("EXOC6B")] <- 1
EntireCohort_pathvar$EXOC6B[EntireCohort_pathvar$GENE2 %in% c("EXOC6B")] <- 1
EntireCohort_pathvar$EXOC6B[EntireCohort_pathvar$GENE3 %in% c("EXOC6B")] <- 1
EntireCohort_pathvar$EXOC6B[EntireCohort_pathvar$GENE4 %in% c("EXOC6B")] <- 1
EntireCohort_pathvar$EXOC6B[EntireCohort_pathvar$GENE5 %in% c("EXOC6B")] <- 1
EntireCohort_pathvar$EXOC6B[EntireCohort_pathvar$GENE6 %in% c("EXOC6B")] <- 1
EntireCohort_pathvar$EXOC6B[EntireCohort_pathvar$GENE7 %in% c("EXOC6B")] <- 1
EntireCohort_pathvar$EXOC6B[EntireCohort_pathvar$GENE8 %in% c("EXOC6B")] <- 1
EntireCohort_pathvar$EXOC6B[EntireCohort_pathvar$GENE9 %in% c("EXOC6B")] <- 1

EntireCohort_pathvar$EXOC6B[is.na(EntireCohort_pathvar$EXOC6B)] <- 0
EntireCohort_pathvar$EXOC6B <- as.numeric(EntireCohort_pathvar$EXOC6B)
unique(EntireCohort_pathvar$EXOC6B)
# 0 1

table(EntireCohort_pathvar$EXOC6B)
#  0      1 
# 466016   3783




### 'RAB33B' 
EntireCohort_pathvar$RAB33B <- NA
EntireCohort_pathvar$RAB33B[EntireCohort_pathvar$GENE1 %in% c("RAB33B")] <- 1
EntireCohort_pathvar$RAB33B[EntireCohort_pathvar$GENE2 %in% c("RAB33B")] <- 1
EntireCohort_pathvar$RAB33B[EntireCohort_pathvar$GENE3 %in% c("RAB33B")] <- 1
EntireCohort_pathvar$RAB33B[EntireCohort_pathvar$GENE4 %in% c("RAB33B")] <- 1
EntireCohort_pathvar$RAB33B[EntireCohort_pathvar$GENE5 %in% c("RAB33B")] <- 1
EntireCohort_pathvar$RAB33B[EntireCohort_pathvar$GENE6 %in% c("RAB33B")] <- 1
EntireCohort_pathvar$RAB33B[EntireCohort_pathvar$GENE7 %in% c("RAB33B")] <- 1
EntireCohort_pathvar$RAB33B[EntireCohort_pathvar$GENE8 %in% c("RAB33B")] <- 1
EntireCohort_pathvar$RAB33B[EntireCohort_pathvar$GENE9 %in% c("RAB33B")] <- 1

EntireCohort_pathvar$RAB33B[is.na(EntireCohort_pathvar$RAB33B)] <- 0
EntireCohort_pathvar$RAB33B <- as.numeric(EntireCohort_pathvar$RAB33B)
unique(EntireCohort_pathvar$RAB33B)
# 0 1 

table(EntireCohort_pathvar$RAB33B)
#  0      1 
# 465932   3867




### 'TRAPPC6B'
EntireCohort_pathvar$TRAPPC6B <- NA
EntireCohort_pathvar$TRAPPC6B[EntireCohort_pathvar$GENE1 %in% c("TRAPPC6B")] <- 1
EntireCohort_pathvar$TRAPPC6B[EntireCohort_pathvar$GENE2 %in% c("TRAPPC6B")] <- 1
EntireCohort_pathvar$TRAPPC6B[EntireCohort_pathvar$GENE3 %in% c("TRAPPC6B")] <- 1
EntireCohort_pathvar$TRAPPC6B[EntireCohort_pathvar$GENE4 %in% c("TRAPPC6B")] <- 1
EntireCohort_pathvar$TRAPPC6B[EntireCohort_pathvar$GENE5 %in% c("TRAPPC6B")] <- 1
EntireCohort_pathvar$TRAPPC6B[EntireCohort_pathvar$GENE6 %in% c("TRAPPC6B")] <- 1
EntireCohort_pathvar$TRAPPC6B[EntireCohort_pathvar$GENE7 %in% c("TRAPPC6B")] <- 1
EntireCohort_pathvar$TRAPPC6B[EntireCohort_pathvar$GENE8 %in% c("TRAPPC6B")] <- 1
EntireCohort_pathvar$TRAPPC6B[EntireCohort_pathvar$GENE9 %in% c("TRAPPC6B")] <- 1

EntireCohort_pathvar$TRAPPC6B[is.na(EntireCohort_pathvar$TRAPPC6B)] <- 0
EntireCohort_pathvar$TRAPPC6B <- as.numeric(EntireCohort_pathvar$TRAPPC6B)
unique(EntireCohort_pathvar$TRAPPC6B)
# 0 1

table(EntireCohort_pathvar$TRAPPC6B)
# 0      1 
# 469287    512




### 'ARCN1'  
EntireCohort_pathvar$ARCN1 <- NA
EntireCohort_pathvar$ARCN1[EntireCohort_pathvar$GENE1 %in% c("ARCN1")] <- 1
EntireCohort_pathvar$ARCN1[EntireCohort_pathvar$GENE2 %in% c("ARCN1")] <- 1
EntireCohort_pathvar$ARCN1[EntireCohort_pathvar$GENE3 %in% c("ARCN1")] <- 1
EntireCohort_pathvar$ARCN1[EntireCohort_pathvar$GENE4 %in% c("ARCN1")] <- 1
EntireCohort_pathvar$ARCN1[EntireCohort_pathvar$GENE5 %in% c("ARCN1")] <- 1
EntireCohort_pathvar$ARCN1[EntireCohort_pathvar$GENE6 %in% c("ARCN1")] <- 1
EntireCohort_pathvar$ARCN1[EntireCohort_pathvar$GENE7 %in% c("ARCN1")] <- 1
EntireCohort_pathvar$ARCN1[EntireCohort_pathvar$GENE8 %in% c("ARCN1")] <- 1
EntireCohort_pathvar$ARCN1[EntireCohort_pathvar$GENE9 %in% c("ARCN1")] <- 1

EntireCohort_pathvar$ARCN1[is.na(EntireCohort_pathvar$ARCN1)] <- 0
EntireCohort_pathvar$ARCN1 <- as.numeric(EntireCohort_pathvar$ARCN1)
unique(EntireCohort_pathvar$ARCN1)
# 0 1

table(EntireCohort_pathvar$ARCN1)
#  0      1 
# 469344    455 




### 'COPB2' 
EntireCohort_pathvar$COPB2 <- NA
EntireCohort_pathvar$COPB2[EntireCohort_pathvar$GENE1 %in% c("COPB2")] <- 1
EntireCohort_pathvar$COPB2[EntireCohort_pathvar$GENE2 %in% c("COPB2")] <- 1
EntireCohort_pathvar$COPB2[EntireCohort_pathvar$GENE3 %in% c("COPB2")] <- 1
EntireCohort_pathvar$COPB2[EntireCohort_pathvar$GENE4 %in% c("COPB2")] <- 1
EntireCohort_pathvar$COPB2[EntireCohort_pathvar$GENE5 %in% c("COPB2")] <- 1
EntireCohort_pathvar$COPB2[EntireCohort_pathvar$GENE6 %in% c("COPB2")] <- 1
EntireCohort_pathvar$COPB2[EntireCohort_pathvar$GENE7 %in% c("COPB2")] <- 1
EntireCohort_pathvar$COPB2[EntireCohort_pathvar$GENE8 %in% c("COPB2")] <- 1
EntireCohort_pathvar$COPB2[EntireCohort_pathvar$GENE9 %in% c("COPB2")] <- 1

EntireCohort_pathvar$COPB2[is.na(EntireCohort_pathvar$COPB2)] <- 0
EntireCohort_pathvar$COPB2 <- as.numeric(EntireCohort_pathvar$COPB2)
unique(EntireCohort_pathvar$COPB2)
# 0 1

table(EntireCohort_pathvar$COPB2)
#  0      1 
# 468434   1365 




### 'CTSA'
EntireCohort_pathvar$CTSA <- NA
EntireCohort_pathvar$CTSA[EntireCohort_pathvar$GENE1 %in% c("CTSA")] <- 1
EntireCohort_pathvar$CTSA[EntireCohort_pathvar$GENE2 %in% c("CTSA")] <- 1
EntireCohort_pathvar$CTSA[EntireCohort_pathvar$GENE3 %in% c("CTSA")] <- 1
EntireCohort_pathvar$CTSA[EntireCohort_pathvar$GENE4 %in% c("CTSA")] <- 1
EntireCohort_pathvar$CTSA[EntireCohort_pathvar$GENE5 %in% c("CTSA")] <- 1
EntireCohort_pathvar$CTSA[EntireCohort_pathvar$GENE6 %in% c("CTSA")] <- 1
EntireCohort_pathvar$CTSA[EntireCohort_pathvar$GENE7 %in% c("CTSA")] <- 1
EntireCohort_pathvar$CTSA[EntireCohort_pathvar$GENE8 %in% c("CTSA")] <- 1
EntireCohort_pathvar$CTSA[EntireCohort_pathvar$GENE9 %in% c("CTSA")] <- 1

EntireCohort_pathvar$CTSA[is.na(EntireCohort_pathvar$CTSA)] <- 0
EntireCohort_pathvar$CTSA <- as.numeric(EntireCohort_pathvar$CTSA)
unique(EntireCohort_pathvar$CTSA)
# 0 1

table(EntireCohort_pathvar$CTSA)
#  0      1 
# 468185   1614 





### 'GALNS' 
EntireCohort_pathvar$GALNS <- NA
EntireCohort_pathvar$GALNS[EntireCohort_pathvar$GENE1 %in% c("GALNS")] <- 1
EntireCohort_pathvar$GALNS[EntireCohort_pathvar$GENE2 %in% c("GALNS")] <- 1
EntireCohort_pathvar$GALNS[EntireCohort_pathvar$GENE3 %in% c("GALNS")] <- 1
EntireCohort_pathvar$GALNS[EntireCohort_pathvar$GENE4 %in% c("GALNS")] <- 1
EntireCohort_pathvar$GALNS[EntireCohort_pathvar$GENE5 %in% c("GALNS")] <- 1
EntireCohort_pathvar$GALNS[EntireCohort_pathvar$GENE6 %in% c("GALNS")] <- 1
EntireCohort_pathvar$GALNS[EntireCohort_pathvar$GENE7 %in% c("GALNS")] <- 1
EntireCohort_pathvar$GALNS[EntireCohort_pathvar$GENE8 %in% c("GALNS")] <- 1
EntireCohort_pathvar$GALNS[EntireCohort_pathvar$GENE9 %in% c("GALNS")] <- 1

EntireCohort_pathvar$GALNS[is.na(EntireCohort_pathvar$GALNS)] <- 0
EntireCohort_pathvar$GALNS <- as.numeric(EntireCohort_pathvar$GALNS)
unique(EntireCohort_pathvar$GALNS)
# 0 1

table(EntireCohort_pathvar$GALNS)
# 0      1 
# 468826    973





### 'GLB1'
EntireCohort_pathvar$GLB1 <- NA
EntireCohort_pathvar$GLB1[EntireCohort_pathvar$GENE1 %in% c("GLB1")] <- 1
EntireCohort_pathvar$GLB1[EntireCohort_pathvar$GENE2 %in% c("GLB1")] <- 1
EntireCohort_pathvar$GLB1[EntireCohort_pathvar$GENE3 %in% c("GLB1")] <- 1
EntireCohort_pathvar$GLB1[EntireCohort_pathvar$GENE4 %in% c("GLB1")] <- 1
EntireCohort_pathvar$GLB1[EntireCohort_pathvar$GENE5 %in% c("GLB1")] <- 1
EntireCohort_pathvar$GLB1[EntireCohort_pathvar$GENE6 %in% c("GLB1")] <- 1
EntireCohort_pathvar$GLB1[EntireCohort_pathvar$GENE7 %in% c("GLB1")] <- 1
EntireCohort_pathvar$GLB1[EntireCohort_pathvar$GENE8 %in% c("GLB1")] <- 1
EntireCohort_pathvar$GLB1[EntireCohort_pathvar$GENE9 %in% c("GLB1")] <- 1

EntireCohort_pathvar$GLB1[is.na(EntireCohort_pathvar$GLB1)] <- 0
EntireCohort_pathvar$GLB1 <- as.numeric(EntireCohort_pathvar$GLB1)
unique(EntireCohort_pathvar$GLB1)
# 0 1

table(EntireCohort_pathvar$GLB1)
#  0      1 
# 468373   1426




### 'GNPTAB'
EntireCohort_pathvar$GNPTAB <- NA
EntireCohort_pathvar$GNPTAB[EntireCohort_pathvar$GENE1 %in% c("GNPTAB")] <- 1
EntireCohort_pathvar$GNPTAB[EntireCohort_pathvar$GENE2 %in% c("GNPTAB")] <- 1
EntireCohort_pathvar$GNPTAB[EntireCohort_pathvar$GENE3 %in% c("GNPTAB")] <- 1
EntireCohort_pathvar$GNPTAB[EntireCohort_pathvar$GENE4 %in% c("GNPTAB")] <- 1
EntireCohort_pathvar$GNPTAB[EntireCohort_pathvar$GENE5 %in% c("GNPTAB")] <- 1
EntireCohort_pathvar$GNPTAB[EntireCohort_pathvar$GENE6 %in% c("GNPTAB")] <- 1
EntireCohort_pathvar$GNPTAB[EntireCohort_pathvar$GENE7 %in% c("GNPTAB")] <- 1
EntireCohort_pathvar$GNPTAB[EntireCohort_pathvar$GENE8 %in% c("GNPTAB")] <- 1
EntireCohort_pathvar$GNPTAB[EntireCohort_pathvar$GENE9 %in% c("GNPTAB")] <- 1

EntireCohort_pathvar$GNPTAB[is.na(EntireCohort_pathvar$GNPTAB)] <- 0
EntireCohort_pathvar$GNPTAB <- as.numeric(EntireCohort_pathvar$GNPTAB)
unique(EntireCohort_pathvar$GNPTAB)
# 0 1

table(EntireCohort_pathvar$GNPTAB)
#  0      1 
# 467660   2139 




### 'GUSB'
EntireCohort_pathvar$GUSB <- NA
EntireCohort_pathvar$GUSB[EntireCohort_pathvar$GENE1 %in% c("GUSB")] <- 1
EntireCohort_pathvar$GUSB[EntireCohort_pathvar$GENE2 %in% c("GUSB")] <- 1
EntireCohort_pathvar$GUSB[EntireCohort_pathvar$GENE3 %in% c("GUSB")] <- 1
EntireCohort_pathvar$GUSB[EntireCohort_pathvar$GENE4 %in% c("GUSB")] <- 1
EntireCohort_pathvar$GUSB[EntireCohort_pathvar$GENE5 %in% c("GUSB")] <- 1
EntireCohort_pathvar$GUSB[EntireCohort_pathvar$GENE6 %in% c("GUSB")] <- 1
EntireCohort_pathvar$GUSB[EntireCohort_pathvar$GENE7 %in% c("GUSB")] <- 1
EntireCohort_pathvar$GUSB[EntireCohort_pathvar$GENE8 %in% c("GUSB")] <- 1
EntireCohort_pathvar$GUSB[EntireCohort_pathvar$GENE9 %in% c("GUSB")] <- 1

EntireCohort_pathvar$GUSB[is.na(EntireCohort_pathvar$GUSB)] <- 0
EntireCohort_pathvar$GUSB <- as.numeric(EntireCohort_pathvar$GUSB)
unique(EntireCohort_pathvar$GUSB)
# 0 1

table(EntireCohort_pathvar$GUSB)
#  0      1 
# 469264    535




### 'ARSB' 
EntireCohort_pathvar$ARSB <- NA
EntireCohort_pathvar$ARSB[EntireCohort_pathvar$GENE1 %in% c("ARSB")] <- 1
EntireCohort_pathvar$ARSB[EntireCohort_pathvar$GENE2 %in% c("ARSB")] <- 1
EntireCohort_pathvar$ARSB[EntireCohort_pathvar$GENE3 %in% c("ARSB")] <- 1
EntireCohort_pathvar$ARSB[EntireCohort_pathvar$GENE4 %in% c("ARSB")] <- 1
EntireCohort_pathvar$ARSB[EntireCohort_pathvar$GENE5 %in% c("ARSB")] <- 1
EntireCohort_pathvar$ARSB[EntireCohort_pathvar$GENE6 %in% c("ARSB")] <- 1
EntireCohort_pathvar$ARSB[EntireCohort_pathvar$GENE7 %in% c("ARSB")] <- 1
EntireCohort_pathvar$ARSB[EntireCohort_pathvar$GENE8 %in% c("ARSB")] <- 1
EntireCohort_pathvar$ARSB[EntireCohort_pathvar$GENE9 %in% c("ARSB")] <- 1

EntireCohort_pathvar$ARSB[is.na(EntireCohort_pathvar$ARSB)] <- 0
EntireCohort_pathvar$ARSB <- as.numeric(EntireCohort_pathvar$ARSB)
unique(EntireCohort_pathvar$ARSB)
# 0 1

table(EntireCohort_pathvar$ARSB)
#  0      1 
# 468174   1625




### 'CUL7'
EntireCohort_pathvar$CUL7 <- NA
EntireCohort_pathvar$CUL7[EntireCohort_pathvar$GENE1 %in% c("CUL7")] <- 1
EntireCohort_pathvar$CUL7[EntireCohort_pathvar$GENE2 %in% c("CUL7")] <- 1
EntireCohort_pathvar$CUL7[EntireCohort_pathvar$GENE3 %in% c("CUL7")] <- 1
EntireCohort_pathvar$CUL7[EntireCohort_pathvar$GENE4 %in% c("CUL7")] <- 1
EntireCohort_pathvar$CUL7[EntireCohort_pathvar$GENE5 %in% c("CUL7")] <- 1
EntireCohort_pathvar$CUL7[EntireCohort_pathvar$GENE6 %in% c("CUL7")] <- 1
EntireCohort_pathvar$CUL7[EntireCohort_pathvar$GENE7 %in% c("CUL7")] <- 1
EntireCohort_pathvar$CUL7[EntireCohort_pathvar$GENE8 %in% c("CUL7")] <- 1
EntireCohort_pathvar$CUL7[EntireCohort_pathvar$GENE9 %in% c("CUL7")] <- 1

EntireCohort_pathvar$CUL7[is.na(EntireCohort_pathvar$CUL7)] <- 0
EntireCohort_pathvar$CUL7 <- as.numeric(EntireCohort_pathvar$CUL7)
unique(EntireCohort_pathvar$CUL7)
# 0 1

table(EntireCohort_pathvar$CUL7)
#  0      1 
# 465659   4140 




### 'CCDC8' 
EntireCohort_pathvar$CCDC8 <- NA
EntireCohort_pathvar$CCDC8[EntireCohort_pathvar$GENE1 %in% c("CCDC8")] <- 1
EntireCohort_pathvar$CCDC8[EntireCohort_pathvar$GENE2 %in% c("CCDC8")] <- 1
EntireCohort_pathvar$CCDC8[EntireCohort_pathvar$GENE3 %in% c("CCDC8")] <- 1
EntireCohort_pathvar$CCDC8[EntireCohort_pathvar$GENE4 %in% c("CCDC8")] <- 1
EntireCohort_pathvar$CCDC8[EntireCohort_pathvar$GENE5 %in% c("CCDC8")] <- 1
EntireCohort_pathvar$CCDC8[EntireCohort_pathvar$GENE6 %in% c("CCDC8")] <- 1
EntireCohort_pathvar$CCDC8[EntireCohort_pathvar$GENE7 %in% c("CCDC8")] <- 1
EntireCohort_pathvar$CCDC8[EntireCohort_pathvar$GENE8 %in% c("CCDC8")] <- 1
EntireCohort_pathvar$CCDC8[EntireCohort_pathvar$GENE9 %in% c("CCDC8")] <- 1

EntireCohort_pathvar$CCDC8[is.na(EntireCohort_pathvar$CCDC8)] <- 0
EntireCohort_pathvar$CCDC8 <- as.numeric(EntireCohort_pathvar$CCDC8)
unique(EntireCohort_pathvar$CCDC8)
# 0 1

table(EntireCohort_pathvar$CCDC8)
#  0      1 
# 469359    440






### 'TBC1D23'
EntireCohort_pathvar$TBC1D23 <- NA
EntireCohort_pathvar$TBC1D23[EntireCohort_pathvar$GENE1 %in% c("TBC1D23")] <- 1
EntireCohort_pathvar$TBC1D23[EntireCohort_pathvar$GENE2 %in% c("TBC1D23")] <- 1
EntireCohort_pathvar$TBC1D23[EntireCohort_pathvar$GENE3 %in% c("TBC1D23")] <- 1
EntireCohort_pathvar$TBC1D23[EntireCohort_pathvar$GENE4 %in% c("TBC1D23")] <- 1
EntireCohort_pathvar$TBC1D23[EntireCohort_pathvar$GENE5 %in% c("TBC1D23")] <- 1
EntireCohort_pathvar$TBC1D23[EntireCohort_pathvar$GENE6 %in% c("TBC1D23")] <- 1
EntireCohort_pathvar$TBC1D23[EntireCohort_pathvar$GENE7 %in% c("TBC1D23")] <- 1
EntireCohort_pathvar$TBC1D23[EntireCohort_pathvar$GENE8 %in% c("TBC1D23")] <- 1
EntireCohort_pathvar$TBC1D23[EntireCohort_pathvar$GENE9 %in% c("TBC1D23")] <- 1

EntireCohort_pathvar$TBC1D23[is.na(EntireCohort_pathvar$TBC1D23)] <- 0
EntireCohort_pathvar$TBC1D23 <- as.numeric(EntireCohort_pathvar$TBC1D23)
unique(EntireCohort_pathvar$TBC1D23)
# 0 1

table(EntireCohort_pathvar$TBC1D23)
# 0      1 
# 468466   1333





### 'UBR1'
EntireCohort_pathvar$UBR1 <- NA
EntireCohort_pathvar$UBR1[EntireCohort_pathvar$GENE1 %in% c("UBR1")] <- 1
EntireCohort_pathvar$UBR1[EntireCohort_pathvar$GENE2 %in% c("UBR1")] <- 1
EntireCohort_pathvar$UBR1[EntireCohort_pathvar$GENE3 %in% c("UBR1")] <- 1
EntireCohort_pathvar$UBR1[EntireCohort_pathvar$GENE4 %in% c("UBR1")] <- 1
EntireCohort_pathvar$UBR1[EntireCohort_pathvar$GENE5 %in% c("UBR1")] <- 1
EntireCohort_pathvar$UBR1[EntireCohort_pathvar$GENE6 %in% c("UBR1")] <- 1
EntireCohort_pathvar$UBR1[EntireCohort_pathvar$GENE7 %in% c("UBR1")] <- 1
EntireCohort_pathvar$UBR1[EntireCohort_pathvar$GENE8 %in% c("UBR1")] <- 1
EntireCohort_pathvar$UBR1[EntireCohort_pathvar$GENE9 %in% c("UBR1")] <- 1

EntireCohort_pathvar$UBR1[is.na(EntireCohort_pathvar$UBR1)] <- 0
EntireCohort_pathvar$UBR1 <- as.numeric(EntireCohort_pathvar$UBR1)
unique(EntireCohort_pathvar$UBR1)
# 0 1

table(EntireCohort_pathvar$UBR1)
#  0      1 
# 467734   2065




### 'STAMBP'   
EntireCohort_pathvar$STAMBP <- NA
EntireCohort_pathvar$STAMBP[EntireCohort_pathvar$GENE1 %in% c("STAMBP")] <- 1
EntireCohort_pathvar$STAMBP[EntireCohort_pathvar$GENE2 %in% c("STAMBP")] <- 1
EntireCohort_pathvar$STAMBP[EntireCohort_pathvar$GENE3 %in% c("STAMBP")] <- 1
EntireCohort_pathvar$STAMBP[EntireCohort_pathvar$GENE4 %in% c("STAMBP")] <- 1
EntireCohort_pathvar$STAMBP[EntireCohort_pathvar$GENE5 %in% c("STAMBP")] <- 1
EntireCohort_pathvar$STAMBP[EntireCohort_pathvar$GENE6 %in% c("STAMBP")] <- 1
EntireCohort_pathvar$STAMBP[EntireCohort_pathvar$GENE7 %in% c("STAMBP")] <- 1
EntireCohort_pathvar$STAMBP[EntireCohort_pathvar$GENE8 %in% c("STAMBP")] <- 1
EntireCohort_pathvar$STAMBP[EntireCohort_pathvar$GENE9 %in% c("STAMBP")] <- 1

EntireCohort_pathvar$STAMBP[is.na(EntireCohort_pathvar$STAMBP)] <- 0
EntireCohort_pathvar$STAMBP <- as.numeric(EntireCohort_pathvar$STAMBP)
unique(EntireCohort_pathvar$STAMBP)
# 0 1

table(EntireCohort_pathvar$STAMBP)
#  0      1 
# 469157    642







# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)







# MGSgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

# tag by MGS gene
# could separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")

EntireCohort_pathvar$MGSgene <- NA
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE1 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE2 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE3 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE4 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE5 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE6 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE7 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE8 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE9 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1

EntireCohort_pathvar$MGSgene[is.na(EntireCohort_pathvar$MGSgene)] <- 0
EntireCohort_pathvar$MGSgene <- as.numeric(EntireCohort_pathvar$MGSgene)
unique(EntireCohort_pathvar$MGSgene)
# 0 1

table(EntireCohort_pathvar$MGSgene)
# 0      1 
# 448986  20813



# now to make a column for each of the MGS genes
### ORC1
EntireCohort_pathvar$ORC1 <- NA
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE1 %in% c("ORC1")] <- 1
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE2 %in% c("ORC1")] <- 1
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE3 %in% c("ORC1")] <- 1
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE4 %in% c("ORC1")] <- 1
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE5 %in% c("ORC1")] <- 1
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE6 %in% c("ORC1")] <- 1
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE7 %in% c("ORC1")] <- 1
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE8 %in% c("ORC1")] <- 1
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE9 %in% c("ORC1")] <- 1

EntireCohort_pathvar$ORC1[is.na(EntireCohort_pathvar$ORC1)] <- 0
EntireCohort_pathvar$ORC1 <- as.numeric(EntireCohort_pathvar$ORC1)
unique(EntireCohort_pathvar$ORC1)
# 0 1

table(EntireCohort_pathvar$ORC1)
# 0      1 
# 468182   1617

### ORC4
EntireCohort_pathvar$ORC4 <- NA
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE1 %in% c("ORC4")] <- 1
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE2 %in% c("ORC4")] <- 1
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE3 %in% c("ORC4")] <- 1
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE4 %in% c("ORC4")] <- 1
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE5 %in% c("ORC4")] <- 1
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE6 %in% c("ORC4")] <- 1
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE7 %in% c("ORC4")] <- 1
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE8 %in% c("ORC4")] <- 1
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE9 %in% c("ORC4")] <- 1

EntireCohort_pathvar$ORC4[is.na(EntireCohort_pathvar$ORC4)] <- 0
EntireCohort_pathvar$ORC4 <- as.numeric(EntireCohort_pathvar$ORC4)
unique(EntireCohort_pathvar$ORC4)
# 0 1

table(EntireCohort_pathvar$ORC4)
# 0      1 
# 466067   3732 

### ORC6
EntireCohort_pathvar$ORC6 <- NA
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE1 %in% c("ORC6")] <- 1
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE2 %in% c("ORC6")] <- 1
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE3 %in% c("ORC6")] <- 1
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE4 %in% c("ORC6")] <- 1
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE5 %in% c("ORC6")] <- 1
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE6 %in% c("ORC6")] <- 1
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE7 %in% c("ORC6")] <- 1
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE8 %in% c("ORC6")] <- 1
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE9 %in% c("ORC6")] <- 1

EntireCohort_pathvar$ORC6[is.na(EntireCohort_pathvar$ORC6)] <- 0
EntireCohort_pathvar$ORC6 <- as.numeric(EntireCohort_pathvar$ORC6)
unique(EntireCohort_pathvar$ORC6)
# 0 1

table(EntireCohort_pathvar$ORC6)
# 0      1 
# 468884    915

### CDT1
EntireCohort_pathvar$CDT1 <- NA
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE1 %in% c("CDT1")] <- 1
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE2 %in% c("CDT1")] <- 1
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE3 %in% c("CDT1")] <- 1
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE4 %in% c("CDT1")] <- 1
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE5 %in% c("CDT1")] <- 1
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE6 %in% c("CDT1")] <- 1
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE7 %in% c("CDT1")] <- 1
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE8 %in% c("CDT1")] <- 1
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE9 %in% c("CDT1")] <- 1

EntireCohort_pathvar$CDT1[is.na(EntireCohort_pathvar$CDT1)] <- 0
EntireCohort_pathvar$CDT1 <- as.numeric(EntireCohort_pathvar$CDT1)
unique(EntireCohort_pathvar$CDT1)
# 0 1

table(EntireCohort_pathvar$CDT1)
# 0      1 
# 468174   1625

### CDC6
EntireCohort_pathvar$CDC6 <- NA
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE1 %in% c("CDC6")] <- 1
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE2 %in% c("CDC6")] <- 1
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE3 %in% c("CDC6")] <- 1
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE4 %in% c("CDC6")] <- 1
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE5 %in% c("CDC6")] <- 1
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE6 %in% c("CDC6")] <- 1
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE7 %in% c("CDC6")] <- 1
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE8 %in% c("CDC6")] <- 1
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE9 %in% c("CDC6")] <- 1

EntireCohort_pathvar$CDC6[is.na(EntireCohort_pathvar$CDC6)] <- 0
EntireCohort_pathvar$CDC6 <- as.numeric(EntireCohort_pathvar$CDC6)
unique(EntireCohort_pathvar$CDC6)
# 0 1

table(EntireCohort_pathvar$CDC6)
# 0      1 
# 469347    452

### CDC45
EntireCohort_pathvar$CDC45 <- NA
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE1 %in% c("CDC45")] <- 1
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE2 %in% c("CDC45")] <- 1
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE3 %in% c("CDC45")] <- 1
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE4 %in% c("CDC45")] <- 1
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE5 %in% c("CDC45")] <- 1
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE6 %in% c("CDC45")] <- 1
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE7 %in% c("CDC45")] <- 1
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE8 %in% c("CDC45")] <- 1
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE9 %in% c("CDC45")] <- 1

EntireCohort_pathvar$CDC45[is.na(EntireCohort_pathvar$CDC45)] <- 0
EntireCohort_pathvar$CDC45 <- as.numeric(EntireCohort_pathvar$CDC45)
unique(EntireCohort_pathvar$CDC45)
# 0 1

table(EntireCohort_pathvar$CDC45)
# 0      1 
# 467756   2043 

### DONSON
EntireCohort_pathvar$DONSON <- NA
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE1 %in% c("DONSON")] <- 1
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE2 %in% c("DONSON")] <- 1
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE3 %in% c("DONSON")] <- 1
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE4 %in% c("DONSON")] <- 1
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE5 %in% c("DONSON")] <- 1
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE6 %in% c("DONSON")] <- 1
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE7 %in% c("DONSON")] <- 1
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE8 %in% c("DONSON")] <- 1
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE9 %in% c("DONSON")] <- 1

EntireCohort_pathvar$DONSON[is.na(EntireCohort_pathvar$DONSON)] <- 0
EntireCohort_pathvar$DONSON <- as.numeric(EntireCohort_pathvar$DONSON)
unique(EntireCohort_pathvar$DONSON)
# 0 1

table(EntireCohort_pathvar$DONSON)
# 0      1 
# 466255   3544 

### GINS2
EntireCohort_pathvar$GINS2 <- NA
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE1 %in% c("GINS2")] <- 1
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE2 %in% c("GINS2")] <- 1
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE3 %in% c("GINS2")] <- 1
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE4 %in% c("GINS2")] <- 1
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE5 %in% c("GINS2")] <- 1
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE6 %in% c("GINS2")] <- 1
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE7 %in% c("GINS2")] <- 1
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE8 %in% c("GINS2")] <- 1
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE9 %in% c("GINS2")] <- 1

EntireCohort_pathvar$GINS2[is.na(EntireCohort_pathvar$GINS2)] <- 0
EntireCohort_pathvar$GINS2 <- as.numeric(EntireCohort_pathvar$GINS2)
unique(EntireCohort_pathvar$GINS2)
# 0 1

table(EntireCohort_pathvar$GINS2)
# 0      1 
# 469245    554 

### GINS3
EntireCohort_pathvar$GINS3 <- NA
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE1 %in% c("GINS3")] <- 1
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE2 %in% c("GINS3")] <- 1
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE3 %in% c("GINS3")] <- 1
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE4 %in% c("GINS3")] <- 1
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE5 %in% c("GINS3")] <- 1
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE6 %in% c("GINS3")] <- 1
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE7 %in% c("GINS3")] <- 1
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE8 %in% c("GINS3")] <- 1
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE9 %in% c("GINS3")] <- 1

EntireCohort_pathvar$GINS3[is.na(EntireCohort_pathvar$GINS3)] <- 0
EntireCohort_pathvar$GINS3 <- as.numeric(EntireCohort_pathvar$GINS3)
unique(EntireCohort_pathvar$GINS3)
# 0 1

table(EntireCohort_pathvar$GINS3)
# 0      1 
# 469413    386 

### MCM3
EntireCohort_pathvar$MCM3 <- NA
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE1 %in% c("MCM3")] <- 1
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE2 %in% c("MCM3")] <- 1
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE3 %in% c("MCM3")] <- 1
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE4 %in% c("MCM3")] <- 1
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE5 %in% c("MCM3")] <- 1
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE6 %in% c("MCM3")] <- 1
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE7 %in% c("MCM3")] <- 1
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE8 %in% c("MCM3")] <- 1
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE9 %in% c("MCM3")] <- 1

EntireCohort_pathvar$MCM3[is.na(EntireCohort_pathvar$MCM3)] <- 0
EntireCohort_pathvar$MCM3 <- as.numeric(EntireCohort_pathvar$MCM3)
unique(EntireCohort_pathvar$MCM3)
# 0 1

table(EntireCohort_pathvar$MCM3)
# 0      1 
# 467331   2468

### MCM5
EntireCohort_pathvar$MCM5 <- NA
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE1 %in% c("MCM5")] <- 1
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE2 %in% c("MCM5")] <- 1
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE3 %in% c("MCM5")] <- 1
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE4 %in% c("MCM5")] <- 1
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE5 %in% c("MCM5")] <- 1
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE6 %in% c("MCM5")] <- 1
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE7 %in% c("MCM5")] <- 1
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE8 %in% c("MCM5")] <- 1
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE9 %in% c("MCM5")] <- 1

EntireCohort_pathvar$MCM5[is.na(EntireCohort_pathvar$MCM5)] <- 0
EntireCohort_pathvar$MCM5 <- as.numeric(EntireCohort_pathvar$MCM5)
unique(EntireCohort_pathvar$MCM5)
# 0 1

table(EntireCohort_pathvar$MCM5)
# 0      1 
# 467771   2028

### MCM7
EntireCohort_pathvar$MCM7 <- NA
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE1 %in% c("MCM7")] <- 1
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE2 %in% c("MCM7")] <- 1
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE3 %in% c("MCM7")] <- 1
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE4 %in% c("MCM7")] <- 1
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE5 %in% c("MCM7")] <- 1
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE6 %in% c("MCM7")] <- 1
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE7 %in% c("MCM7")] <- 1
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE8 %in% c("MCM7")] <- 1
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE9 %in% c("MCM7")] <- 1

EntireCohort_pathvar$MCM7[is.na(EntireCohort_pathvar$MCM7)] <- 0
EntireCohort_pathvar$MCM7 <- as.numeric(EntireCohort_pathvar$MCM7)
unique(EntireCohort_pathvar$MCM7)
# 0 1

table(EntireCohort_pathvar$MCM7)
# 0      1 
# 467936   1863


# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)




# Mitogenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

# go thr all the genes to check that they're present
grep("AIFM1", EntireCohort_pathvar$GENE)
# all seem to be present, except for:
# AIFM1


# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$Mitogene <- NA
EntireCohort_pathvar$Mitogene[EntireCohort_pathvar$GENE1 %in% c('PYCR1', 'PISD', 'SLC25A19', 'BCS1L', 'ATAD3A', 'PDHA1', 'PAM16', 'GPX4')] <- 1
EntireCohort_pathvar$Mitogene[EntireCohort_pathvar$GENE2 %in% c('PYCR1', 'PISD', 'SLC25A19', 'BCS1L', 'ATAD3A', 'PDHA1', 'PAM16', 'GPX4')] <- 1
EntireCohort_pathvar$Mitogene[EntireCohort_pathvar$GENE3 %in% c('PYCR1', 'PISD', 'SLC25A19', 'BCS1L', 'ATAD3A', 'PDHA1', 'PAM16', 'GPX4')] <- 1
EntireCohort_pathvar$Mitogene[EntireCohort_pathvar$GENE4 %in% c('PYCR1', 'PISD', 'SLC25A19', 'BCS1L', 'ATAD3A', 'PDHA1', 'PAM16', 'GPX4')] <- 1
EntireCohort_pathvar$Mitogene[EntireCohort_pathvar$GENE5 %in% c('PYCR1', 'PISD', 'SLC25A19', 'BCS1L', 'ATAD3A', 'PDHA1', 'PAM16', 'GPX4')] <- 1
EntireCohort_pathvar$Mitogene[EntireCohort_pathvar$GENE6 %in% c('PYCR1', 'PISD', 'SLC25A19', 'BCS1L', 'ATAD3A', 'PDHA1', 'PAM16', 'GPX4')] <- 1
EntireCohort_pathvar$Mitogene[EntireCohort_pathvar$GENE7 %in% c('PYCR1', 'PISD', 'SLC25A19', 'BCS1L', 'ATAD3A', 'PDHA1', 'PAM16', 'GPX4')] <- 1
EntireCohort_pathvar$Mitogene[EntireCohort_pathvar$GENE8 %in% c('PYCR1', 'PISD', 'SLC25A19', 'BCS1L', 'ATAD3A', 'PDHA1', 'PAM16', 'GPX4')] <- 1
EntireCohort_pathvar$Mitogene[EntireCohort_pathvar$GENE9 %in% c('PYCR1', 'PISD', 'SLC25A19', 'BCS1L', 'ATAD3A', 'PDHA1', 'PAM16', 'GPX4')] <- 1


EntireCohort_pathvar$Mitogene[is.na(EntireCohort_pathvar$Mitogene)] <- 0
EntireCohort_pathvar$Mitogene <- as.numeric(EntireCohort_pathvar$Mitogene)
unique(EntireCohort_pathvar$Mitogene)
# 0 1

table(EntireCohort_pathvar$Mitogene)
#  0      1 
# 463590   6209 




# now to make a column for each of the Mitogenes
### 'PYCR1'
EntireCohort_pathvar$PYCR1 <- NA
EntireCohort_pathvar$PYCR1[EntireCohort_pathvar$GENE1 %in% c("PYCR1")] <- 1
EntireCohort_pathvar$PYCR1[EntireCohort_pathvar$GENE2 %in% c("PYCR1")] <- 1
EntireCohort_pathvar$PYCR1[EntireCohort_pathvar$GENE3 %in% c("PYCR1")] <- 1
EntireCohort_pathvar$PYCR1[EntireCohort_pathvar$GENE4 %in% c("PYCR1")] <- 1
EntireCohort_pathvar$PYCR1[EntireCohort_pathvar$GENE5 %in% c("PYCR1")] <- 1
EntireCohort_pathvar$PYCR1[EntireCohort_pathvar$GENE6 %in% c("PYCR1")] <- 1
EntireCohort_pathvar$PYCR1[EntireCohort_pathvar$GENE7 %in% c("PYCR1")] <- 1
EntireCohort_pathvar$PYCR1[EntireCohort_pathvar$GENE8 %in% c("PYCR1")] <- 1
EntireCohort_pathvar$PYCR1[EntireCohort_pathvar$GENE9 %in% c("PYCR1")] <- 1

EntireCohort_pathvar$PYCR1[is.na(EntireCohort_pathvar$PYCR1)] <- 0
EntireCohort_pathvar$PYCR1 <- as.numeric(EntireCohort_pathvar$PYCR1)
unique(EntireCohort_pathvar$PYCR1)
# 0 1

table(EntireCohort_pathvar$PYCR1)
#  0      1 
# 469352    447 





### 'PISD' 
EntireCohort_pathvar$PISD <- NA
EntireCohort_pathvar$PISD[EntireCohort_pathvar$GENE1 %in% c("PISD")] <- 1
EntireCohort_pathvar$PISD[EntireCohort_pathvar$GENE2 %in% c("PISD")] <- 1
EntireCohort_pathvar$PISD[EntireCohort_pathvar$GENE3 %in% c("PISD")] <- 1
EntireCohort_pathvar$PISD[EntireCohort_pathvar$GENE4 %in% c("PISD")] <- 1
EntireCohort_pathvar$PISD[EntireCohort_pathvar$GENE5 %in% c("PISD")] <- 1
EntireCohort_pathvar$PISD[EntireCohort_pathvar$GENE6 %in% c("PISD")] <- 1
EntireCohort_pathvar$PISD[EntireCohort_pathvar$GENE7 %in% c("PISD")] <- 1
EntireCohort_pathvar$PISD[EntireCohort_pathvar$GENE8 %in% c("PISD")] <- 1
EntireCohort_pathvar$PISD[EntireCohort_pathvar$GENE9 %in% c("PISD")] <- 1

EntireCohort_pathvar$PISD[is.na(EntireCohort_pathvar$PISD)] <- 0
EntireCohort_pathvar$PISD <- as.numeric(EntireCohort_pathvar$PISD)
unique(EntireCohort_pathvar$PISD)
# 0 1

table(EntireCohort_pathvar$PISD)
#  0      1 
# 468006   1793




### 'SLC25A19'
EntireCohort_pathvar$SLC25A19 <- NA
EntireCohort_pathvar$SLC25A19[EntireCohort_pathvar$GENE1 %in% c("SLC25A19")] <- 1
EntireCohort_pathvar$SLC25A19[EntireCohort_pathvar$GENE2 %in% c("SLC25A19")] <- 1
EntireCohort_pathvar$SLC25A19[EntireCohort_pathvar$GENE3 %in% c("SLC25A19")] <- 1
EntireCohort_pathvar$SLC25A19[EntireCohort_pathvar$GENE4 %in% c("SLC25A19")] <- 1
EntireCohort_pathvar$SLC25A19[EntireCohort_pathvar$GENE5 %in% c("SLC25A19")] <- 1
EntireCohort_pathvar$SLC25A19[EntireCohort_pathvar$GENE6 %in% c("SLC25A19")] <- 1
EntireCohort_pathvar$SLC25A19[EntireCohort_pathvar$GENE7 %in% c("SLC25A19")] <- 1
EntireCohort_pathvar$SLC25A19[EntireCohort_pathvar$GENE8 %in% c("SLC25A19")] <- 1
EntireCohort_pathvar$SLC25A19[EntireCohort_pathvar$GENE9 %in% c("SLC25A19")] <- 1

EntireCohort_pathvar$SLC25A19[is.na(EntireCohort_pathvar$SLC25A19)] <- 0
EntireCohort_pathvar$SLC25A19 <- as.numeric(EntireCohort_pathvar$SLC25A19)
unique(EntireCohort_pathvar$SLC25A19)
# 0 1

table(EntireCohort_pathvar$SLC25A19)
#  0      1 
# 469481    318




### 'BCS1L'
EntireCohort_pathvar$BCS1L <- NA
EntireCohort_pathvar$BCS1L[EntireCohort_pathvar$GENE1 %in% c("BCS1L")] <- 1
EntireCohort_pathvar$BCS1L[EntireCohort_pathvar$GENE2 %in% c("BCS1L")] <- 1
EntireCohort_pathvar$BCS1L[EntireCohort_pathvar$GENE3 %in% c("BCS1L")] <- 1
EntireCohort_pathvar$BCS1L[EntireCohort_pathvar$GENE4 %in% c("BCS1L")] <- 1
EntireCohort_pathvar$BCS1L[EntireCohort_pathvar$GENE5 %in% c("BCS1L")] <- 1
EntireCohort_pathvar$BCS1L[EntireCohort_pathvar$GENE6 %in% c("BCS1L")] <- 1
EntireCohort_pathvar$BCS1L[EntireCohort_pathvar$GENE7 %in% c("BCS1L")] <- 1
EntireCohort_pathvar$BCS1L[EntireCohort_pathvar$GENE8 %in% c("BCS1L")] <- 1
EntireCohort_pathvar$BCS1L[EntireCohort_pathvar$GENE9 %in% c("BCS1L")] <- 1

EntireCohort_pathvar$BCS1L[is.na(EntireCohort_pathvar$BCS1L)] <- 0
EntireCohort_pathvar$BCS1L <- as.numeric(EntireCohort_pathvar$BCS1L)
unique(EntireCohort_pathvar$BCS1L)
# 0 1

table(EntireCohort_pathvar$BCS1L)
#   0      1 
# 468347   1452




### 'ATAD3A'
EntireCohort_pathvar$ATAD3A <- NA
EntireCohort_pathvar$ATAD3A[EntireCohort_pathvar$GENE1 %in% c("ATAD3A")] <- 1
EntireCohort_pathvar$ATAD3A[EntireCohort_pathvar$GENE2 %in% c("ATAD3A")] <- 1
EntireCohort_pathvar$ATAD3A[EntireCohort_pathvar$GENE3 %in% c("ATAD3A")] <- 1
EntireCohort_pathvar$ATAD3A[EntireCohort_pathvar$GENE4 %in% c("ATAD3A")] <- 1
EntireCohort_pathvar$ATAD3A[EntireCohort_pathvar$GENE5 %in% c("ATAD3A")] <- 1
EntireCohort_pathvar$ATAD3A[EntireCohort_pathvar$GENE6 %in% c("ATAD3A")] <- 1
EntireCohort_pathvar$ATAD3A[EntireCohort_pathvar$GENE7 %in% c("ATAD3A")] <- 1
EntireCohort_pathvar$ATAD3A[EntireCohort_pathvar$GENE8 %in% c("ATAD3A")] <- 1
EntireCohort_pathvar$ATAD3A[EntireCohort_pathvar$GENE9 %in% c("ATAD3A")] <- 1

EntireCohort_pathvar$ATAD3A[is.na(EntireCohort_pathvar$ATAD3A)] <- 0
EntireCohort_pathvar$ATAD3A <- as.numeric(EntireCohort_pathvar$ATAD3A)
unique(EntireCohort_pathvar$ATAD3A)
# 0 1

table(EntireCohort_pathvar$ATAD3A)
# 0      1 
# 468314   1485




### 'PDHA1' 
EntireCohort_pathvar$PDHA1 <- NA
EntireCohort_pathvar$PDHA1[EntireCohort_pathvar$GENE1 %in% c("PDHA1")] <- 1
EntireCohort_pathvar$PDHA1[EntireCohort_pathvar$GENE2 %in% c("PDHA1")] <- 1
EntireCohort_pathvar$PDHA1[EntireCohort_pathvar$GENE3 %in% c("PDHA1")] <- 1
EntireCohort_pathvar$PDHA1[EntireCohort_pathvar$GENE4 %in% c("PDHA1")] <- 1
EntireCohort_pathvar$PDHA1[EntireCohort_pathvar$GENE5 %in% c("PDHA1")] <- 1
EntireCohort_pathvar$PDHA1[EntireCohort_pathvar$GENE6 %in% c("PDHA1")] <- 1
EntireCohort_pathvar$PDHA1[EntireCohort_pathvar$GENE7 %in% c("PDHA1")] <- 1
EntireCohort_pathvar$PDHA1[EntireCohort_pathvar$GENE8 %in% c("PDHA1")] <- 1
EntireCohort_pathvar$PDHA1[EntireCohort_pathvar$GENE9 %in% c("PDHA1")] <- 1

EntireCohort_pathvar$PDHA1[is.na(EntireCohort_pathvar$PDHA1)] <- 0
EntireCohort_pathvar$PDHA1 <- as.numeric(EntireCohort_pathvar$PDHA1)
unique(EntireCohort_pathvar$PDHA1)
# 0 1

table(EntireCohort_pathvar$PDHA1)
#  0      1 
# 469597    202 




### 'PAM16' 
EntireCohort_pathvar$PAM16 <- NA
EntireCohort_pathvar$PAM16[EntireCohort_pathvar$GENE1 %in% c("PAM16")] <- 1
EntireCohort_pathvar$PAM16[EntireCohort_pathvar$GENE2 %in% c("PAM16")] <- 1
EntireCohort_pathvar$PAM16[EntireCohort_pathvar$GENE3 %in% c("PAM16")] <- 1
EntireCohort_pathvar$PAM16[EntireCohort_pathvar$GENE4 %in% c("PAM16")] <- 1
EntireCohort_pathvar$PAM16[EntireCohort_pathvar$GENE5 %in% c("PAM16")] <- 1
EntireCohort_pathvar$PAM16[EntireCohort_pathvar$GENE6 %in% c("PAM16")] <- 1
EntireCohort_pathvar$PAM16[EntireCohort_pathvar$GENE7 %in% c("PAM16")] <- 1
EntireCohort_pathvar$PAM16[EntireCohort_pathvar$GENE8 %in% c("PAM16")] <- 1
EntireCohort_pathvar$PAM16[EntireCohort_pathvar$GENE9 %in% c("PAM16")] <- 1

EntireCohort_pathvar$PAM16[is.na(EntireCohort_pathvar$PAM16)] <- 0
EntireCohort_pathvar$PAM16 <- as.numeric(EntireCohort_pathvar$PAM16)
unique(EntireCohort_pathvar$PAM16)
# 0 1

table(EntireCohort_pathvar$PAM16)
#   0      1 
# 469517    282 




### 'GPX4'  
EntireCohort_pathvar$GPX4 <- NA
EntireCohort_pathvar$GPX4[EntireCohort_pathvar$GENE1 %in% c("GPX4")] <- 1
EntireCohort_pathvar$GPX4[EntireCohort_pathvar$GENE2 %in% c("GPX4")] <- 1
EntireCohort_pathvar$GPX4[EntireCohort_pathvar$GENE3 %in% c("GPX4")] <- 1
EntireCohort_pathvar$GPX4[EntireCohort_pathvar$GENE4 %in% c("GPX4")] <- 1
EntireCohort_pathvar$GPX4[EntireCohort_pathvar$GENE5 %in% c("GPX4")] <- 1
EntireCohort_pathvar$GPX4[EntireCohort_pathvar$GENE6 %in% c("GPX4")] <- 1
EntireCohort_pathvar$GPX4[EntireCohort_pathvar$GENE7 %in% c("GPX4")] <- 1
EntireCohort_pathvar$GPX4[EntireCohort_pathvar$GENE8 %in% c("GPX4")] <- 1
EntireCohort_pathvar$GPX4[EntireCohort_pathvar$GENE9 %in% c("GPX4")] <- 1

EntireCohort_pathvar$GPX4[is.na(EntireCohort_pathvar$GPX4)] <- 0
EntireCohort_pathvar$GPX4 <- as.numeric(EntireCohort_pathvar$GPX4)
unique(EntireCohort_pathvar$GPX4)
# 0 1

table(EntireCohort_pathvar$GPX4)
#  0      1 
# 469531    268




### 'AIFM1' 
# no variants present





# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)








# RASgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

# go thr all the genes to check that they're present
grep("FGD1", EntireCohort_pathvar$GENE)
# all seem to be present, except for the following:
# FGFR1 = 0 variants
# RALA = 0 variants
# MRAS = 0 variants

# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$RASgene <- NA
EntireCohort_pathvar$RASgene[EntireCohort_pathvar$GENE1 %in% c('PHEX', 'FGFR3', 'RASA2', 'SHOC2', 'SPRED2', 'FGD1')] <- 1
EntireCohort_pathvar$RASgene[EntireCohort_pathvar$GENE2 %in% c('PHEX', 'FGFR3', 'RASA2', 'SHOC2', 'SPRED2', 'FGD1')] <- 1
EntireCohort_pathvar$RASgene[EntireCohort_pathvar$GENE3 %in% c('PHEX', 'FGFR3', 'RASA2', 'SHOC2', 'SPRED2', 'FGD1')] <- 1
EntireCohort_pathvar$RASgene[EntireCohort_pathvar$GENE4 %in% c('PHEX', 'FGFR3', 'RASA2', 'SHOC2', 'SPRED2', 'FGD1')] <- 1
EntireCohort_pathvar$RASgene[EntireCohort_pathvar$GENE5 %in% c('PHEX', 'FGFR3', 'RASA2', 'SHOC2', 'SPRED2', 'FGD1')] <- 1
EntireCohort_pathvar$RASgene[EntireCohort_pathvar$GENE6 %in% c('PHEX', 'FGFR3', 'RASA2', 'SHOC2', 'SPRED2', 'FGD1')] <- 1
EntireCohort_pathvar$RASgene[EntireCohort_pathvar$GENE7 %in% c('PHEX', 'FGFR3', 'RASA2', 'SHOC2', 'SPRED2', 'FGD1')] <- 1
EntireCohort_pathvar$RASgene[EntireCohort_pathvar$GENE8 %in% c('PHEX', 'FGFR3', 'RASA2', 'SHOC2', 'SPRED2', 'FGD1')] <- 1
EntireCohort_pathvar$RASgene[EntireCohort_pathvar$GENE9 %in% c('PHEX', 'FGFR3', 'RASA2', 'SHOC2', 'SPRED2', 'FGD1')] <- 1


EntireCohort_pathvar$RASgene[is.na(EntireCohort_pathvar$RASgene)] <- 0
EntireCohort_pathvar$RASgene <- as.numeric(EntireCohort_pathvar$RASgene)
unique(EntireCohort_pathvar$RASgene)
# 0 1

table(EntireCohort_pathvar$RASgene)
#      0      1 
# 467585   2214




# now to make a column for each of the RASgenes
### 'PHEX'
EntireCohort_pathvar$PHEX <- NA
EntireCohort_pathvar$PHEX[EntireCohort_pathvar$GENE1 %in% c("PHEX")] <- 1
EntireCohort_pathvar$PHEX[EntireCohort_pathvar$GENE2 %in% c("PHEX")] <- 1
EntireCohort_pathvar$PHEX[EntireCohort_pathvar$GENE3 %in% c("PHEX")] <- 1
EntireCohort_pathvar$PHEX[EntireCohort_pathvar$GENE4 %in% c("PHEX")] <- 1
EntireCohort_pathvar$PHEX[EntireCohort_pathvar$GENE5 %in% c("PHEX")] <- 1
EntireCohort_pathvar$PHEX[EntireCohort_pathvar$GENE6 %in% c("PHEX")] <- 1
EntireCohort_pathvar$PHEX[EntireCohort_pathvar$GENE7 %in% c("PHEX")] <- 1
EntireCohort_pathvar$PHEX[EntireCohort_pathvar$GENE8 %in% c("PHEX")] <- 1
EntireCohort_pathvar$PHEX[EntireCohort_pathvar$GENE9 %in% c("PHEX")] <- 1

EntireCohort_pathvar$PHEX[is.na(EntireCohort_pathvar$PHEX)] <- 0
EntireCohort_pathvar$PHEX <- as.numeric(EntireCohort_pathvar$PHEX)
unique(EntireCohort_pathvar$PHEX)
# 0 1

table(EntireCohort_pathvar$PHEX)
# 0      1 
# 469302    497





### 'FGFR1'
# no variants present




### 'FGFR3'
EntireCohort_pathvar$FGFR3 <- NA
EntireCohort_pathvar$FGFR3[EntireCohort_pathvar$GENE1 %in% c("FGFR3")] <- 1
EntireCohort_pathvar$FGFR3[EntireCohort_pathvar$GENE2 %in% c("FGFR3")] <- 1
EntireCohort_pathvar$FGFR3[EntireCohort_pathvar$GENE3 %in% c("FGFR3")] <- 1
EntireCohort_pathvar$FGFR3[EntireCohort_pathvar$GENE4 %in% c("FGFR3")] <- 1
EntireCohort_pathvar$FGFR3[EntireCohort_pathvar$GENE5 %in% c("FGFR3")] <- 1
EntireCohort_pathvar$FGFR3[EntireCohort_pathvar$GENE6 %in% c("FGFR3")] <- 1
EntireCohort_pathvar$FGFR3[EntireCohort_pathvar$GENE7 %in% c("FGFR3")] <- 1
EntireCohort_pathvar$FGFR3[EntireCohort_pathvar$GENE8 %in% c("FGFR3")] <- 1
EntireCohort_pathvar$FGFR3[EntireCohort_pathvar$GENE9 %in% c("FGFR3")] <- 1

EntireCohort_pathvar$FGFR3[is.na(EntireCohort_pathvar$FGFR3)] <- 0
EntireCohort_pathvar$FGFR3 <- as.numeric(EntireCohort_pathvar$FGFR3)
unique(EntireCohort_pathvar$FGFR3)
# 0 1

table(EntireCohort_pathvar$FGFR3)
# 0      1 
# 469785     14




### 'RASA2'
EntireCohort_pathvar$RASA2 <- NA
EntireCohort_pathvar$RASA2[EntireCohort_pathvar$GENE1 %in% c("RASA2")] <- 1
EntireCohort_pathvar$RASA2[EntireCohort_pathvar$GENE2 %in% c("RASA2")] <- 1
EntireCohort_pathvar$RASA2[EntireCohort_pathvar$GENE3 %in% c("RASA2")] <- 1
EntireCohort_pathvar$RASA2[EntireCohort_pathvar$GENE4 %in% c("RASA2")] <- 1
EntireCohort_pathvar$RASA2[EntireCohort_pathvar$GENE5 %in% c("RASA2")] <- 1
EntireCohort_pathvar$RASA2[EntireCohort_pathvar$GENE6 %in% c("RASA2")] <- 1
EntireCohort_pathvar$RASA2[EntireCohort_pathvar$GENE7 %in% c("RASA2")] <- 1
EntireCohort_pathvar$RASA2[EntireCohort_pathvar$GENE8 %in% c("RASA2")] <- 1
EntireCohort_pathvar$RASA2[EntireCohort_pathvar$GENE9 %in% c("RASA2")] <- 1

EntireCohort_pathvar$RASA2[is.na(EntireCohort_pathvar$RASA2)] <- 0
EntireCohort_pathvar$RASA2 <- as.numeric(EntireCohort_pathvar$RASA2)
unique(EntireCohort_pathvar$RASA2)
# 0 1

table(EntireCohort_pathvar$RASA2)
# 0      1 
# 468894    905




### 'RALA'
# no variants present




### 'SHOC2'
EntireCohort_pathvar$SHOC2 <- NA
EntireCohort_pathvar$SHOC2[EntireCohort_pathvar$GENE1 %in% c("SHOC2")] <- 1
EntireCohort_pathvar$SHOC2[EntireCohort_pathvar$GENE2 %in% c("SHOC2")] <- 1
EntireCohort_pathvar$SHOC2[EntireCohort_pathvar$GENE3 %in% c("SHOC2")] <- 1
EntireCohort_pathvar$SHOC2[EntireCohort_pathvar$GENE4 %in% c("SHOC2")] <- 1
EntireCohort_pathvar$SHOC2[EntireCohort_pathvar$GENE5 %in% c("SHOC2")] <- 1
EntireCohort_pathvar$SHOC2[EntireCohort_pathvar$GENE6 %in% c("SHOC2")] <- 1
EntireCohort_pathvar$SHOC2[EntireCohort_pathvar$GENE7 %in% c("SHOC2")] <- 1
EntireCohort_pathvar$SHOC2[EntireCohort_pathvar$GENE8 %in% c("SHOC2")] <- 1
EntireCohort_pathvar$SHOC2[EntireCohort_pathvar$GENE9 %in% c("SHOC2")] <- 1

EntireCohort_pathvar$SHOC2[is.na(EntireCohort_pathvar$SHOC2)] <- 0
EntireCohort_pathvar$SHOC2 <- as.numeric(EntireCohort_pathvar$SHOC2)
unique(EntireCohort_pathvar$SHOC2)
# 0 1

table(EntireCohort_pathvar$SHOC2)
# 0      1 
# 469793      6




### 'SPRED2'
EntireCohort_pathvar$SPRED2 <- NA
EntireCohort_pathvar$SPRED2[EntireCohort_pathvar$GENE1 %in% c("SPRED2")] <- 1
EntireCohort_pathvar$SPRED2[EntireCohort_pathvar$GENE2 %in% c("SPRED2")] <- 1
EntireCohort_pathvar$SPRED2[EntireCohort_pathvar$GENE3 %in% c("SPRED2")] <- 1
EntireCohort_pathvar$SPRED2[EntireCohort_pathvar$GENE4 %in% c("SPRED2")] <- 1
EntireCohort_pathvar$SPRED2[EntireCohort_pathvar$GENE5 %in% c("SPRED2")] <- 1
EntireCohort_pathvar$SPRED2[EntireCohort_pathvar$GENE6 %in% c("SPRED2")] <- 1
EntireCohort_pathvar$SPRED2[EntireCohort_pathvar$GENE7 %in% c("SPRED2")] <- 1
EntireCohort_pathvar$SPRED2[EntireCohort_pathvar$GENE8 %in% c("SPRED2")] <- 1
EntireCohort_pathvar$SPRED2[EntireCohort_pathvar$GENE9 %in% c("SPRED2")] <- 1

EntireCohort_pathvar$SPRED2[is.na(EntireCohort_pathvar$SPRED2)] <- 0
EntireCohort_pathvar$SPRED2 <- as.numeric(EntireCohort_pathvar$SPRED2)
unique(EntireCohort_pathvar$SPRED2)
# 0 1

table(EntireCohort_pathvar$SPRED2)
# 0      1 
# 469261    538




### 'MRAS'
# no variants present




### 'FGD1'
EntireCohort_pathvar$FGD1 <- NA
EntireCohort_pathvar$FGD1[EntireCohort_pathvar$GENE1 %in% c("FGD1")] <- 1
EntireCohort_pathvar$FGD1[EntireCohort_pathvar$GENE2 %in% c("FGD1")] <- 1
EntireCohort_pathvar$FGD1[EntireCohort_pathvar$GENE3 %in% c("FGD1")] <- 1
EntireCohort_pathvar$FGD1[EntireCohort_pathvar$GENE4 %in% c("FGD1")] <- 1
EntireCohort_pathvar$FGD1[EntireCohort_pathvar$GENE5 %in% c("FGD1")] <- 1
EntireCohort_pathvar$FGD1[EntireCohort_pathvar$GENE6 %in% c("FGD1")] <- 1
EntireCohort_pathvar$FGD1[EntireCohort_pathvar$GENE7 %in% c("FGD1")] <- 1
EntireCohort_pathvar$FGD1[EntireCohort_pathvar$GENE8 %in% c("FGD1")] <- 1
EntireCohort_pathvar$FGD1[EntireCohort_pathvar$GENE9 %in% c("FGD1")] <- 1

EntireCohort_pathvar$FGD1[is.na(EntireCohort_pathvar$FGD1)] <- 0
EntireCohort_pathvar$FGD1 <- as.numeric(EntireCohort_pathvar$FGD1)
unique(EntireCohort_pathvar$FGD1)
# 0 1

table(EntireCohort_pathvar$FGD1)
# 0      1 
# 469539    260







# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)







# RNAPgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

# go thr all the genes to check that they're present
grep("RRP7A", EntireCohort_pathvar$GENE)
# all seem to be present, except for:
# RNU4ATAC = 0
# YARS1 = 0



# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$RNAPgene <- NA
EntireCohort_pathvar$RNAPgene[EntireCohort_pathvar$GENE1 %in% c('CRIPT', 'WDR73', 'TOE1', 'THOC6', 'EXOSC3', 'PUF60', 'PPIL1', 'CDC40', 'CLP1', 'AARS1', 'CTU2', 'GON7', 'LAGE3', 'OSGEP', 'POP1', 'SEPSECS', 'THG1L', 'TP53RK', 'TPRKB', 'TRMT1', 'TRMT10A', 'TSEN15', 'TSEN2', 'TSEN54', 'WDR4', 'YRDC', 'RRP7A')] <- 1
EntireCohort_pathvar$RNAPgene[EntireCohort_pathvar$GENE2 %in% c('CRIPT', 'WDR73', 'TOE1', 'THOC6', 'EXOSC3', 'PUF60', 'PPIL1', 'CDC40', 'CLP1', 'AARS1', 'CTU2', 'GON7', 'LAGE3', 'OSGEP', 'POP1', 'SEPSECS', 'THG1L', 'TP53RK', 'TPRKB', 'TRMT1', 'TRMT10A', 'TSEN15', 'TSEN2', 'TSEN54', 'WDR4', 'YRDC', 'RRP7A')] <- 1
EntireCohort_pathvar$RNAPgene[EntireCohort_pathvar$GENE3 %in% c('CRIPT', 'WDR73', 'TOE1', 'THOC6', 'EXOSC3', 'PUF60', 'PPIL1', 'CDC40', 'CLP1', 'AARS1', 'CTU2', 'GON7', 'LAGE3', 'OSGEP', 'POP1', 'SEPSECS', 'THG1L', 'TP53RK', 'TPRKB', 'TRMT1', 'TRMT10A', 'TSEN15', 'TSEN2', 'TSEN54', 'WDR4', 'YRDC', 'RRP7A')] <- 1
EntireCohort_pathvar$RNAPgene[EntireCohort_pathvar$GENE4 %in% c('CRIPT', 'WDR73', 'TOE1', 'THOC6', 'EXOSC3', 'PUF60', 'PPIL1', 'CDC40', 'CLP1', 'AARS1', 'CTU2', 'GON7', 'LAGE3', 'OSGEP', 'POP1', 'SEPSECS', 'THG1L', 'TP53RK', 'TPRKB', 'TRMT1', 'TRMT10A', 'TSEN15', 'TSEN2', 'TSEN54', 'WDR4', 'YRDC', 'RRP7A')] <- 1
EntireCohort_pathvar$RNAPgene[EntireCohort_pathvar$GENE5 %in% c('CRIPT', 'WDR73', 'TOE1', 'THOC6', 'EXOSC3', 'PUF60', 'PPIL1', 'CDC40', 'CLP1', 'AARS1', 'CTU2', 'GON7', 'LAGE3', 'OSGEP', 'POP1', 'SEPSECS', 'THG1L', 'TP53RK', 'TPRKB', 'TRMT1', 'TRMT10A', 'TSEN15', 'TSEN2', 'TSEN54', 'WDR4', 'YRDC', 'RRP7A')] <- 1
EntireCohort_pathvar$RNAPgene[EntireCohort_pathvar$GENE6 %in% c('CRIPT', 'WDR73', 'TOE1', 'THOC6', 'EXOSC3', 'PUF60', 'PPIL1', 'CDC40', 'CLP1', 'AARS1', 'CTU2', 'GON7', 'LAGE3', 'OSGEP', 'POP1', 'SEPSECS', 'THG1L', 'TP53RK', 'TPRKB', 'TRMT1', 'TRMT10A', 'TSEN15', 'TSEN2', 'TSEN54', 'WDR4', 'YRDC', 'RRP7A')] <- 1
EntireCohort_pathvar$RNAPgene[EntireCohort_pathvar$GENE7 %in% c('CRIPT', 'WDR73', 'TOE1', 'THOC6', 'EXOSC3', 'PUF60', 'PPIL1', 'CDC40', 'CLP1', 'AARS1', 'CTU2', 'GON7', 'LAGE3', 'OSGEP', 'POP1', 'SEPSECS', 'THG1L', 'TP53RK', 'TPRKB', 'TRMT1', 'TRMT10A', 'TSEN15', 'TSEN2', 'TSEN54', 'WDR4', 'YRDC', 'RRP7A')] <- 1
EntireCohort_pathvar$RNAPgene[EntireCohort_pathvar$GENE8 %in% c('CRIPT', 'WDR73', 'TOE1', 'THOC6', 'EXOSC3', 'PUF60', 'PPIL1', 'CDC40', 'CLP1', 'AARS1', 'CTU2', 'GON7', 'LAGE3', 'OSGEP', 'POP1', 'SEPSECS', 'THG1L', 'TP53RK', 'TPRKB', 'TRMT1', 'TRMT10A', 'TSEN15', 'TSEN2', 'TSEN54', 'WDR4', 'YRDC', 'RRP7A')] <- 1
EntireCohort_pathvar$RNAPgene[EntireCohort_pathvar$GENE9 %in% c('CRIPT', 'WDR73', 'TOE1', 'THOC6', 'EXOSC3', 'PUF60', 'PPIL1', 'CDC40', 'CLP1', 'AARS1', 'CTU2', 'GON7', 'LAGE3', 'OSGEP', 'POP1', 'SEPSECS', 'THG1L', 'TP53RK', 'TPRKB', 'TRMT1', 'TRMT10A', 'TSEN15', 'TSEN2', 'TSEN54', 'WDR4', 'YRDC', 'RRP7A')] <- 1


EntireCohort_pathvar$RNAPgene[is.na(EntireCohort_pathvar$RNAPgene)] <- 0
EntireCohort_pathvar$RNAPgene <- as.numeric(EntireCohort_pathvar$RNAPgene)
unique(EntireCohort_pathvar$RNAPgene)
# 1 0

table(EntireCohort_pathvar$RNAPgene)
#  0      1 
# 444092  25707




# now to make a column for each of the RNAPgenes
### 'CRIPT', 
EntireCohort_pathvar$CRIPT <- NA
EntireCohort_pathvar$CRIPT[EntireCohort_pathvar$GENE1 %in% c("CRIPT")] <- 1
EntireCohort_pathvar$CRIPT[EntireCohort_pathvar$GENE2 %in% c("CRIPT")] <- 1
EntireCohort_pathvar$CRIPT[EntireCohort_pathvar$GENE3 %in% c("CRIPT")] <- 1
EntireCohort_pathvar$CRIPT[EntireCohort_pathvar$GENE4 %in% c("CRIPT")] <- 1
EntireCohort_pathvar$CRIPT[EntireCohort_pathvar$GENE5 %in% c("CRIPT")] <- 1
EntireCohort_pathvar$CRIPT[EntireCohort_pathvar$GENE6 %in% c("CRIPT")] <- 1
EntireCohort_pathvar$CRIPT[EntireCohort_pathvar$GENE7 %in% c("CRIPT")] <- 1
EntireCohort_pathvar$CRIPT[EntireCohort_pathvar$GENE8 %in% c("CRIPT")] <- 1
EntireCohort_pathvar$CRIPT[EntireCohort_pathvar$GENE9 %in% c("CRIPT")] <- 1

EntireCohort_pathvar$CRIPT[is.na(EntireCohort_pathvar$CRIPT)] <- 0
EntireCohort_pathvar$CRIPT <- as.numeric(EntireCohort_pathvar$CRIPT)
unique(EntireCohort_pathvar$CRIPT)
# 0 1

table(EntireCohort_pathvar$CRIPT)
# 0      1 
# 469392    407





### 'WDR73', 
EntireCohort_pathvar$WDR73 <- NA
EntireCohort_pathvar$WDR73[EntireCohort_pathvar$GENE1 %in% c("WDR73")] <- 1
EntireCohort_pathvar$WDR73[EntireCohort_pathvar$GENE2 %in% c("WDR73")] <- 1
EntireCohort_pathvar$WDR73[EntireCohort_pathvar$GENE3 %in% c("WDR73")] <- 1
EntireCohort_pathvar$WDR73[EntireCohort_pathvar$GENE4 %in% c("WDR73")] <- 1
EntireCohort_pathvar$WDR73[EntireCohort_pathvar$GENE5 %in% c("WDR73")] <- 1
EntireCohort_pathvar$WDR73[EntireCohort_pathvar$GENE6 %in% c("WDR73")] <- 1
EntireCohort_pathvar$WDR73[EntireCohort_pathvar$GENE7 %in% c("WDR73")] <- 1
EntireCohort_pathvar$WDR73[EntireCohort_pathvar$GENE8 %in% c("WDR73")] <- 1
EntireCohort_pathvar$WDR73[EntireCohort_pathvar$GENE9 %in% c("WDR73")] <- 1

EntireCohort_pathvar$WDR73[is.na(EntireCohort_pathvar$WDR73)] <- 0
EntireCohort_pathvar$WDR73 <- as.numeric(EntireCohort_pathvar$WDR73)
unique(EntireCohort_pathvar$WDR73)
# 0 1

table(EntireCohort_pathvar$WDR73)
#  0      1 
# 469206    593




### 'TOE1', 
EntireCohort_pathvar$TOE1 <- NA
EntireCohort_pathvar$TOE1[EntireCohort_pathvar$GENE1 %in% c("TOE1")] <- 1
EntireCohort_pathvar$TOE1[EntireCohort_pathvar$GENE2 %in% c("TOE1")] <- 1
EntireCohort_pathvar$TOE1[EntireCohort_pathvar$GENE3 %in% c("TOE1")] <- 1
EntireCohort_pathvar$TOE1[EntireCohort_pathvar$GENE4 %in% c("TOE1")] <- 1
EntireCohort_pathvar$TOE1[EntireCohort_pathvar$GENE5 %in% c("TOE1")] <- 1
EntireCohort_pathvar$TOE1[EntireCohort_pathvar$GENE6 %in% c("TOE1")] <- 1
EntireCohort_pathvar$TOE1[EntireCohort_pathvar$GENE7 %in% c("TOE1")] <- 1
EntireCohort_pathvar$TOE1[EntireCohort_pathvar$GENE8 %in% c("TOE1")] <- 1
EntireCohort_pathvar$TOE1[EntireCohort_pathvar$GENE9 %in% c("TOE1")] <- 1

EntireCohort_pathvar$TOE1[is.na(EntireCohort_pathvar$TOE1)] <- 0
EntireCohort_pathvar$TOE1 <- as.numeric(EntireCohort_pathvar$TOE1)
unique(EntireCohort_pathvar$TOE1)
# 0 1

table(EntireCohort_pathvar$TOE1)
#  0      1 
# 468823    976 




### 'THOC6', 
EntireCohort_pathvar$THOC6 <- NA
EntireCohort_pathvar$THOC6[EntireCohort_pathvar$GENE1 %in% c("THOC6")] <- 1
EntireCohort_pathvar$THOC6[EntireCohort_pathvar$GENE2 %in% c("THOC6")] <- 1
EntireCohort_pathvar$THOC6[EntireCohort_pathvar$GENE3 %in% c("THOC6")] <- 1
EntireCohort_pathvar$THOC6[EntireCohort_pathvar$GENE4 %in% c("THOC6")] <- 1
EntireCohort_pathvar$THOC6[EntireCohort_pathvar$GENE5 %in% c("THOC6")] <- 1
EntireCohort_pathvar$THOC6[EntireCohort_pathvar$GENE6 %in% c("THOC6")] <- 1
EntireCohort_pathvar$THOC6[EntireCohort_pathvar$GENE7 %in% c("THOC6")] <- 1
EntireCohort_pathvar$THOC6[EntireCohort_pathvar$GENE8 %in% c("THOC6")] <- 1
EntireCohort_pathvar$THOC6[EntireCohort_pathvar$GENE9 %in% c("THOC6")] <- 1

EntireCohort_pathvar$THOC6[is.na(EntireCohort_pathvar$THOC6)] <- 0
EntireCohort_pathvar$THOC6 <- as.numeric(EntireCohort_pathvar$THOC6)
unique(EntireCohort_pathvar$THOC6)
# 0 1

table(EntireCohort_pathvar$THOC6)
#  0      1 
# 468639   1160 




### 'EXOSC3', 
EntireCohort_pathvar$EXOSC3 <- NA
EntireCohort_pathvar$EXOSC3[EntireCohort_pathvar$GENE1 %in% c("EXOSC3")] <- 1
EntireCohort_pathvar$EXOSC3[EntireCohort_pathvar$GENE2 %in% c("EXOSC3")] <- 1
EntireCohort_pathvar$EXOSC3[EntireCohort_pathvar$GENE3 %in% c("EXOSC3")] <- 1
EntireCohort_pathvar$EXOSC3[EntireCohort_pathvar$GENE4 %in% c("EXOSC3")] <- 1
EntireCohort_pathvar$EXOSC3[EntireCohort_pathvar$GENE5 %in% c("EXOSC3")] <- 1
EntireCohort_pathvar$EXOSC3[EntireCohort_pathvar$GENE6 %in% c("EXOSC3")] <- 1
EntireCohort_pathvar$EXOSC3[EntireCohort_pathvar$GENE7 %in% c("EXOSC3")] <- 1
EntireCohort_pathvar$EXOSC3[EntireCohort_pathvar$GENE8 %in% c("EXOSC3")] <- 1
EntireCohort_pathvar$EXOSC3[EntireCohort_pathvar$GENE9 %in% c("EXOSC3")] <- 1

EntireCohort_pathvar$EXOSC3[is.na(EntireCohort_pathvar$EXOSC3)] <- 0
EntireCohort_pathvar$EXOSC3 <- as.numeric(EntireCohort_pathvar$EXOSC3)
unique(EntireCohort_pathvar$EXOSC3)
# 0 1

table(EntireCohort_pathvar$EXOSC3)
# 0      1 
# 468119   1680




### 'PUF60', 
EntireCohort_pathvar$PUF60 <- NA
EntireCohort_pathvar$PUF60[EntireCohort_pathvar$GENE1 %in% c("PUF60")] <- 1
EntireCohort_pathvar$PUF60[EntireCohort_pathvar$GENE2 %in% c("PUF60")] <- 1
EntireCohort_pathvar$PUF60[EntireCohort_pathvar$GENE3 %in% c("PUF60")] <- 1
EntireCohort_pathvar$PUF60[EntireCohort_pathvar$GENE4 %in% c("PUF60")] <- 1
EntireCohort_pathvar$PUF60[EntireCohort_pathvar$GENE5 %in% c("PUF60")] <- 1
EntireCohort_pathvar$PUF60[EntireCohort_pathvar$GENE6 %in% c("PUF60")] <- 1
EntireCohort_pathvar$PUF60[EntireCohort_pathvar$GENE7 %in% c("PUF60")] <- 1
EntireCohort_pathvar$PUF60[EntireCohort_pathvar$GENE8 %in% c("PUF60")] <- 1
EntireCohort_pathvar$PUF60[EntireCohort_pathvar$GENE9 %in% c("PUF60")] <- 1

EntireCohort_pathvar$PUF60[is.na(EntireCohort_pathvar$PUF60)] <- 0
EntireCohort_pathvar$PUF60 <- as.numeric(EntireCohort_pathvar$PUF60)
unique(EntireCohort_pathvar$PUF60)
# 0 1

table(EntireCohort_pathvar$PUF60)
# 0      1 
# 469696    103




### 'PPIL1', 
EntireCohort_pathvar$PPIL1 <- NA
EntireCohort_pathvar$PPIL1[EntireCohort_pathvar$GENE1 %in% c("PPIL1")] <- 1
EntireCohort_pathvar$PPIL1[EntireCohort_pathvar$GENE2 %in% c("PPIL1")] <- 1
EntireCohort_pathvar$PPIL1[EntireCohort_pathvar$GENE3 %in% c("PPIL1")] <- 1
EntireCohort_pathvar$PPIL1[EntireCohort_pathvar$GENE4 %in% c("PPIL1")] <- 1
EntireCohort_pathvar$PPIL1[EntireCohort_pathvar$GENE5 %in% c("PPIL1")] <- 1
EntireCohort_pathvar$PPIL1[EntireCohort_pathvar$GENE6 %in% c("PPIL1")] <- 1
EntireCohort_pathvar$PPIL1[EntireCohort_pathvar$GENE7 %in% c("PPIL1")] <- 1
EntireCohort_pathvar$PPIL1[EntireCohort_pathvar$GENE8 %in% c("PPIL1")] <- 1
EntireCohort_pathvar$PPIL1[EntireCohort_pathvar$GENE9 %in% c("PPIL1")] <- 1

EntireCohort_pathvar$PPIL1[is.na(EntireCohort_pathvar$PPIL1)] <- 0
EntireCohort_pathvar$PPIL1 <- as.numeric(EntireCohort_pathvar$PPIL1)
unique(EntireCohort_pathvar$PPIL1)
# 0 1

table(EntireCohort_pathvar$PPIL1)
#   0      1 
# 469625    174




### 'RNU4ATAC', 
# no variants present




### 'CDC40', 
EntireCohort_pathvar$CDC40 <- NA
EntireCohort_pathvar$CDC40[EntireCohort_pathvar$GENE1 %in% c("CDC40")] <- 1
EntireCohort_pathvar$CDC40[EntireCohort_pathvar$GENE2 %in% c("CDC40")] <- 1
EntireCohort_pathvar$CDC40[EntireCohort_pathvar$GENE3 %in% c("CDC40")] <- 1
EntireCohort_pathvar$CDC40[EntireCohort_pathvar$GENE4 %in% c("CDC40")] <- 1
EntireCohort_pathvar$CDC40[EntireCohort_pathvar$GENE5 %in% c("CDC40")] <- 1
EntireCohort_pathvar$CDC40[EntireCohort_pathvar$GENE6 %in% c("CDC40")] <- 1
EntireCohort_pathvar$CDC40[EntireCohort_pathvar$GENE7 %in% c("CDC40")] <- 1
EntireCohort_pathvar$CDC40[EntireCohort_pathvar$GENE8 %in% c("CDC40")] <- 1
EntireCohort_pathvar$CDC40[EntireCohort_pathvar$GENE9 %in% c("CDC40")] <- 1

EntireCohort_pathvar$CDC40[is.na(EntireCohort_pathvar$CDC40)] <- 0
EntireCohort_pathvar$CDC40 <- as.numeric(EntireCohort_pathvar$CDC40)
unique(EntireCohort_pathvar$CDC40)
# 0 1

table(EntireCohort_pathvar$CDC40)
# 0      1 
# 469351    448





### 'CLP1',     
EntireCohort_pathvar$CLP1 <- NA
EntireCohort_pathvar$CLP1[EntireCohort_pathvar$GENE1 %in% c("CLP1")] <- 1
EntireCohort_pathvar$CLP1[EntireCohort_pathvar$GENE2 %in% c("CLP1")] <- 1
EntireCohort_pathvar$CLP1[EntireCohort_pathvar$GENE3 %in% c("CLP1")] <- 1
EntireCohort_pathvar$CLP1[EntireCohort_pathvar$GENE4 %in% c("CLP1")] <- 1
EntireCohort_pathvar$CLP1[EntireCohort_pathvar$GENE5 %in% c("CLP1")] <- 1
EntireCohort_pathvar$CLP1[EntireCohort_pathvar$GENE6 %in% c("CLP1")] <- 1
EntireCohort_pathvar$CLP1[EntireCohort_pathvar$GENE7 %in% c("CLP1")] <- 1
EntireCohort_pathvar$CLP1[EntireCohort_pathvar$GENE8 %in% c("CLP1")] <- 1
EntireCohort_pathvar$CLP1[EntireCohort_pathvar$GENE9 %in% c("CLP1")] <- 1

EntireCohort_pathvar$CLP1[is.na(EntireCohort_pathvar$CLP1)] <- 0
EntireCohort_pathvar$CLP1 <- as.numeric(EntireCohort_pathvar$CLP1)
unique(EntireCohort_pathvar$CLP1)
# 0 1

table(EntireCohort_pathvar$CLP1)
#  0      1 
# 469246    553




### 'AARS1', 
EntireCohort_pathvar$AARS1 <- NA
EntireCohort_pathvar$AARS1[EntireCohort_pathvar$GENE1 %in% c("AARS1")] <- 1
EntireCohort_pathvar$AARS1[EntireCohort_pathvar$GENE2 %in% c("AARS1")] <- 1
EntireCohort_pathvar$AARS1[EntireCohort_pathvar$GENE3 %in% c("AARS1")] <- 1
EntireCohort_pathvar$AARS1[EntireCohort_pathvar$GENE4 %in% c("AARS1")] <- 1
EntireCohort_pathvar$AARS1[EntireCohort_pathvar$GENE5 %in% c("AARS1")] <- 1
EntireCohort_pathvar$AARS1[EntireCohort_pathvar$GENE6 %in% c("AARS1")] <- 1
EntireCohort_pathvar$AARS1[EntireCohort_pathvar$GENE7 %in% c("AARS1")] <- 1
EntireCohort_pathvar$AARS1[EntireCohort_pathvar$GENE8 %in% c("AARS1")] <- 1
EntireCohort_pathvar$AARS1[EntireCohort_pathvar$GENE9 %in% c("AARS1")] <- 1

EntireCohort_pathvar$AARS1[is.na(EntireCohort_pathvar$AARS1)] <- 0
EntireCohort_pathvar$AARS1 <- as.numeric(EntireCohort_pathvar$AARS1)
unique(EntireCohort_pathvar$AARS1)
# 0 1

table(EntireCohort_pathvar$AARS1)
# 0      1 
# 465558   4241




### 'CTU2', 
EntireCohort_pathvar$CTU2 <- NA
EntireCohort_pathvar$CTU2[EntireCohort_pathvar$GENE1 %in% c("CTU2")] <- 1
EntireCohort_pathvar$CTU2[EntireCohort_pathvar$GENE2 %in% c("CTU2")] <- 1
EntireCohort_pathvar$CTU2[EntireCohort_pathvar$GENE3 %in% c("CTU2")] <- 1
EntireCohort_pathvar$CTU2[EntireCohort_pathvar$GENE4 %in% c("CTU2")] <- 1
EntireCohort_pathvar$CTU2[EntireCohort_pathvar$GENE5 %in% c("CTU2")] <- 1
EntireCohort_pathvar$CTU2[EntireCohort_pathvar$GENE6 %in% c("CTU2")] <- 1
EntireCohort_pathvar$CTU2[EntireCohort_pathvar$GENE7 %in% c("CTU2")] <- 1
EntireCohort_pathvar$CTU2[EntireCohort_pathvar$GENE8 %in% c("CTU2")] <- 1
EntireCohort_pathvar$CTU2[EntireCohort_pathvar$GENE9 %in% c("CTU2")] <- 1

EntireCohort_pathvar$CTU2[is.na(EntireCohort_pathvar$CTU2)] <- 0
EntireCohort_pathvar$CTU2 <- as.numeric(EntireCohort_pathvar$CTU2)
unique(EntireCohort_pathvar$CTU2)
# 0 1

table(EntireCohort_pathvar$CTU2)
# 0      1 
# 465728   4071




### 'GON7',   
EntireCohort_pathvar$GON7 <- NA
EntireCohort_pathvar$GON7[EntireCohort_pathvar$GENE1 %in% c("GON7")] <- 1
EntireCohort_pathvar$GON7[EntireCohort_pathvar$GENE2 %in% c("GON7")] <- 1
EntireCohort_pathvar$GON7[EntireCohort_pathvar$GENE3 %in% c("GON7")] <- 1
EntireCohort_pathvar$GON7[EntireCohort_pathvar$GENE4 %in% c("GON7")] <- 1
EntireCohort_pathvar$GON7[EntireCohort_pathvar$GENE5 %in% c("GON7")] <- 1
EntireCohort_pathvar$GON7[EntireCohort_pathvar$GENE6 %in% c("GON7")] <- 1
EntireCohort_pathvar$GON7[EntireCohort_pathvar$GENE7 %in% c("GON7")] <- 1
EntireCohort_pathvar$GON7[EntireCohort_pathvar$GENE8 %in% c("GON7")] <- 1
EntireCohort_pathvar$GON7[EntireCohort_pathvar$GENE9 %in% c("GON7")] <- 1

EntireCohort_pathvar$GON7[is.na(EntireCohort_pathvar$GON7)] <- 0
EntireCohort_pathvar$GON7 <- as.numeric(EntireCohort_pathvar$GON7)
unique(EntireCohort_pathvar$GON7)
# 0 1

table(EntireCohort_pathvar$GON7)
# 0      1 
# 469647    152




### 'LAGE3', 
EntireCohort_pathvar$LAGE3 <- NA
EntireCohort_pathvar$LAGE3[EntireCohort_pathvar$GENE1 %in% c("LAGE3")] <- 1
EntireCohort_pathvar$LAGE3[EntireCohort_pathvar$GENE2 %in% c("LAGE3")] <- 1
EntireCohort_pathvar$LAGE3[EntireCohort_pathvar$GENE3 %in% c("LAGE3")] <- 1
EntireCohort_pathvar$LAGE3[EntireCohort_pathvar$GENE4 %in% c("LAGE3")] <- 1
EntireCohort_pathvar$LAGE3[EntireCohort_pathvar$GENE5 %in% c("LAGE3")] <- 1
EntireCohort_pathvar$LAGE3[EntireCohort_pathvar$GENE6 %in% c("LAGE3")] <- 1
EntireCohort_pathvar$LAGE3[EntireCohort_pathvar$GENE7 %in% c("LAGE3")] <- 1
EntireCohort_pathvar$LAGE3[EntireCohort_pathvar$GENE8 %in% c("LAGE3")] <- 1
EntireCohort_pathvar$LAGE3[EntireCohort_pathvar$GENE9 %in% c("LAGE3")] <- 1

EntireCohort_pathvar$LAGE3[is.na(EntireCohort_pathvar$LAGE3)] <- 0
EntireCohort_pathvar$LAGE3 <- as.numeric(EntireCohort_pathvar$LAGE3)
unique(EntireCohort_pathvar$LAGE3)
# 0 1

table(EntireCohort_pathvar$LAGE3)
#  0      1 
# 469797      2




### 'OSGEP', 
EntireCohort_pathvar$OSGEP <- NA
EntireCohort_pathvar$OSGEP[EntireCohort_pathvar$GENE1 %in% c("OSGEP")] <- 1
EntireCohort_pathvar$OSGEP[EntireCohort_pathvar$GENE2 %in% c("OSGEP")] <- 1
EntireCohort_pathvar$OSGEP[EntireCohort_pathvar$GENE3 %in% c("OSGEP")] <- 1
EntireCohort_pathvar$OSGEP[EntireCohort_pathvar$GENE4 %in% c("OSGEP")] <- 1
EntireCohort_pathvar$OSGEP[EntireCohort_pathvar$GENE5 %in% c("OSGEP")] <- 1
EntireCohort_pathvar$OSGEP[EntireCohort_pathvar$GENE6 %in% c("OSGEP")] <- 1
EntireCohort_pathvar$OSGEP[EntireCohort_pathvar$GENE7 %in% c("OSGEP")] <- 1
EntireCohort_pathvar$OSGEP[EntireCohort_pathvar$GENE8 %in% c("OSGEP")] <- 1
EntireCohort_pathvar$OSGEP[EntireCohort_pathvar$GENE9 %in% c("OSGEP")] <- 1

EntireCohort_pathvar$OSGEP[is.na(EntireCohort_pathvar$OSGEP)] <- 0
EntireCohort_pathvar$OSGEP <- as.numeric(EntireCohort_pathvar$OSGEP)
unique(EntireCohort_pathvar$OSGEP)
# 0 1

table(EntireCohort_pathvar$OSGEP)
# 0      1 
# 467944   1855





### 'POP1', 
EntireCohort_pathvar$POP1 <- NA
EntireCohort_pathvar$POP1[EntireCohort_pathvar$GENE1 %in% c("POP1")] <- 1
EntireCohort_pathvar$POP1[EntireCohort_pathvar$GENE2 %in% c("POP1")] <- 1
EntireCohort_pathvar$POP1[EntireCohort_pathvar$GENE3 %in% c("POP1")] <- 1
EntireCohort_pathvar$POP1[EntireCohort_pathvar$GENE4 %in% c("POP1")] <- 1
EntireCohort_pathvar$POP1[EntireCohort_pathvar$GENE5 %in% c("POP1")] <- 1
EntireCohort_pathvar$POP1[EntireCohort_pathvar$GENE6 %in% c("POP1")] <- 1
EntireCohort_pathvar$POP1[EntireCohort_pathvar$GENE7 %in% c("POP1")] <- 1
EntireCohort_pathvar$POP1[EntireCohort_pathvar$GENE8 %in% c("POP1")] <- 1
EntireCohort_pathvar$POP1[EntireCohort_pathvar$GENE9 %in% c("POP1")] <- 1

EntireCohort_pathvar$POP1[is.na(EntireCohort_pathvar$POP1)] <- 0
EntireCohort_pathvar$POP1 <- as.numeric(EntireCohort_pathvar$POP1)
unique(EntireCohort_pathvar$POP1)
# 0 1

table(EntireCohort_pathvar$POP1)
# 0      1 
# 468175   1624





### 'SEPSECS', 
EntireCohort_pathvar$SEPSECS <- NA
EntireCohort_pathvar$SEPSECS[EntireCohort_pathvar$GENE1 %in% c("SEPSECS")] <- 1
EntireCohort_pathvar$SEPSECS[EntireCohort_pathvar$GENE2 %in% c("SEPSECS")] <- 1
EntireCohort_pathvar$SEPSECS[EntireCohort_pathvar$GENE3 %in% c("SEPSECS")] <- 1
EntireCohort_pathvar$SEPSECS[EntireCohort_pathvar$GENE4 %in% c("SEPSECS")] <- 1
EntireCohort_pathvar$SEPSECS[EntireCohort_pathvar$GENE5 %in% c("SEPSECS")] <- 1
EntireCohort_pathvar$SEPSECS[EntireCohort_pathvar$GENE6 %in% c("SEPSECS")] <- 1
EntireCohort_pathvar$SEPSECS[EntireCohort_pathvar$GENE7 %in% c("SEPSECS")] <- 1
EntireCohort_pathvar$SEPSECS[EntireCohort_pathvar$GENE8 %in% c("SEPSECS")] <- 1
EntireCohort_pathvar$SEPSECS[EntireCohort_pathvar$GENE9 %in% c("SEPSECS")] <- 1

EntireCohort_pathvar$SEPSECS[is.na(EntireCohort_pathvar$SEPSECS)] <- 0
EntireCohort_pathvar$SEPSECS <- as.numeric(EntireCohort_pathvar$SEPSECS)
unique(EntireCohort_pathvar$SEPSECS)
# 0 1

table(EntireCohort_pathvar$SEPSECS)
#  0      1 
# 468846    953




### 'THG1L',
EntireCohort_pathvar$THG1L <- NA
EntireCohort_pathvar$THG1L[EntireCohort_pathvar$GENE1 %in% c("THG1L")] <- 1
EntireCohort_pathvar$THG1L[EntireCohort_pathvar$GENE2 %in% c("THG1L")] <- 1
EntireCohort_pathvar$THG1L[EntireCohort_pathvar$GENE3 %in% c("THG1L")] <- 1
EntireCohort_pathvar$THG1L[EntireCohort_pathvar$GENE4 %in% c("THG1L")] <- 1
EntireCohort_pathvar$THG1L[EntireCohort_pathvar$GENE5 %in% c("THG1L")] <- 1
EntireCohort_pathvar$THG1L[EntireCohort_pathvar$GENE6 %in% c("THG1L")] <- 1
EntireCohort_pathvar$THG1L[EntireCohort_pathvar$GENE7 %in% c("THG1L")] <- 1
EntireCohort_pathvar$THG1L[EntireCohort_pathvar$GENE8 %in% c("THG1L")] <- 1
EntireCohort_pathvar$THG1L[EntireCohort_pathvar$GENE9 %in% c("THG1L")] <- 1

EntireCohort_pathvar$THG1L[is.na(EntireCohort_pathvar$THG1L)] <- 0
EntireCohort_pathvar$THG1L <- as.numeric(EntireCohort_pathvar$THG1L)
unique(EntireCohort_pathvar$THG1L)
# 1 0

table(EntireCohort_pathvar$THG1L)
#  0      1 
# 468543   1256 




### 'TP53RK', 
EntireCohort_pathvar$TP53RK <- NA
EntireCohort_pathvar$TP53RK[EntireCohort_pathvar$GENE1 %in% c("TP53RK")] <- 1
EntireCohort_pathvar$TP53RK[EntireCohort_pathvar$GENE2 %in% c("TP53RK")] <- 1
EntireCohort_pathvar$TP53RK[EntireCohort_pathvar$GENE3 %in% c("TP53RK")] <- 1
EntireCohort_pathvar$TP53RK[EntireCohort_pathvar$GENE4 %in% c("TP53RK")] <- 1
EntireCohort_pathvar$TP53RK[EntireCohort_pathvar$GENE5 %in% c("TP53RK")] <- 1
EntireCohort_pathvar$TP53RK[EntireCohort_pathvar$GENE6 %in% c("TP53RK")] <- 1
EntireCohort_pathvar$TP53RK[EntireCohort_pathvar$GENE7 %in% c("TP53RK")] <- 1
EntireCohort_pathvar$TP53RK[EntireCohort_pathvar$GENE8 %in% c("TP53RK")] <- 1
EntireCohort_pathvar$TP53RK[EntireCohort_pathvar$GENE9 %in% c("TP53RK")] <- 1

EntireCohort_pathvar$TP53RK[is.na(EntireCohort_pathvar$TP53RK)] <- 0
EntireCohort_pathvar$TP53RK <- as.numeric(EntireCohort_pathvar$TP53RK)
unique(EntireCohort_pathvar$TP53RK)
# 0 1

table(EntireCohort_pathvar$TP53RK)
# 0      1 
# 469357    442




### 'TPRKB', 
EntireCohort_pathvar$TPRKB <- NA
EntireCohort_pathvar$TPRKB[EntireCohort_pathvar$GENE1 %in% c("TPRKB")] <- 1
EntireCohort_pathvar$TPRKB[EntireCohort_pathvar$GENE2 %in% c("TPRKB")] <- 1
EntireCohort_pathvar$TPRKB[EntireCohort_pathvar$GENE3 %in% c("TPRKB")] <- 1
EntireCohort_pathvar$TPRKB[EntireCohort_pathvar$GENE4 %in% c("TPRKB")] <- 1
EntireCohort_pathvar$TPRKB[EntireCohort_pathvar$GENE5 %in% c("TPRKB")] <- 1
EntireCohort_pathvar$TPRKB[EntireCohort_pathvar$GENE6 %in% c("TPRKB")] <- 1
EntireCohort_pathvar$TPRKB[EntireCohort_pathvar$GENE7 %in% c("TPRKB")] <- 1
EntireCohort_pathvar$TPRKB[EntireCohort_pathvar$GENE8 %in% c("TPRKB")] <- 1
EntireCohort_pathvar$TPRKB[EntireCohort_pathvar$GENE9 %in% c("TPRKB")] <- 1

EntireCohort_pathvar$TPRKB[is.na(EntireCohort_pathvar$TPRKB)] <- 0
EntireCohort_pathvar$TPRKB <- as.numeric(EntireCohort_pathvar$TPRKB)
unique(EntireCohort_pathvar$TPRKB)
# 0 1

table(EntireCohort_pathvar$TPRKB)
#   0      1 
# 469716     83




### 'TRMT1', 
EntireCohort_pathvar$TRMT1 <- NA
EntireCohort_pathvar$TRMT1[EntireCohort_pathvar$GENE1 %in% c("TRMT1")] <- 1
EntireCohort_pathvar$TRMT1[EntireCohort_pathvar$GENE2 %in% c("TRMT1")] <- 1
EntireCohort_pathvar$TRMT1[EntireCohort_pathvar$GENE3 %in% c("TRMT1")] <- 1
EntireCohort_pathvar$TRMT1[EntireCohort_pathvar$GENE4 %in% c("TRMT1")] <- 1
EntireCohort_pathvar$TRMT1[EntireCohort_pathvar$GENE5 %in% c("TRMT1")] <- 1
EntireCohort_pathvar$TRMT1[EntireCohort_pathvar$GENE6 %in% c("TRMT1")] <- 1
EntireCohort_pathvar$TRMT1[EntireCohort_pathvar$GENE7 %in% c("TRMT1")] <- 1
EntireCohort_pathvar$TRMT1[EntireCohort_pathvar$GENE8 %in% c("TRMT1")] <- 1
EntireCohort_pathvar$TRMT1[EntireCohort_pathvar$GENE9 %in% c("TRMT1")] <- 1

EntireCohort_pathvar$TRMT1[is.na(EntireCohort_pathvar$TRMT1)] <- 0
EntireCohort_pathvar$TRMT1 <- as.numeric(EntireCohort_pathvar$TRMT1)
unique(EntireCohort_pathvar$TRMT1)
# 0 1

table(EntireCohort_pathvar$TRMT1)
#  0      1 
# 468636   1163 




### 'TRMT10A', 
EntireCohort_pathvar$TRMT10A <- NA
EntireCohort_pathvar$TRMT10A[EntireCohort_pathvar$GENE1 %in% c("TRMT10A")] <- 1
EntireCohort_pathvar$TRMT10A[EntireCohort_pathvar$GENE2 %in% c("TRMT10A")] <- 1
EntireCohort_pathvar$TRMT10A[EntireCohort_pathvar$GENE3 %in% c("TRMT10A")] <- 1
EntireCohort_pathvar$TRMT10A[EntireCohort_pathvar$GENE4 %in% c("TRMT10A")] <- 1
EntireCohort_pathvar$TRMT10A[EntireCohort_pathvar$GENE5 %in% c("TRMT10A")] <- 1
EntireCohort_pathvar$TRMT10A[EntireCohort_pathvar$GENE6 %in% c("TRMT10A")] <- 1
EntireCohort_pathvar$TRMT10A[EntireCohort_pathvar$GENE7 %in% c("TRMT10A")] <- 1
EntireCohort_pathvar$TRMT10A[EntireCohort_pathvar$GENE8 %in% c("TRMT10A")] <- 1
EntireCohort_pathvar$TRMT10A[EntireCohort_pathvar$GENE9 %in% c("TRMT10A")] <- 1

EntireCohort_pathvar$TRMT10A[is.na(EntireCohort_pathvar$TRMT10A)] <- 0
EntireCohort_pathvar$TRMT10A <- as.numeric(EntireCohort_pathvar$TRMT10A)
unique(EntireCohort_pathvar$TRMT10A)
# 0 1

table(EntireCohort_pathvar$TRMT10A)
#  0      1 
# 468608   1191






### 'TSEN15', 
EntireCohort_pathvar$TSEN15 <- NA
EntireCohort_pathvar$TSEN15[EntireCohort_pathvar$GENE1 %in% c("TSEN15")] <- 1
EntireCohort_pathvar$TSEN15[EntireCohort_pathvar$GENE2 %in% c("TSEN15")] <- 1
EntireCohort_pathvar$TSEN15[EntireCohort_pathvar$GENE3 %in% c("TSEN15")] <- 1
EntireCohort_pathvar$TSEN15[EntireCohort_pathvar$GENE4 %in% c("TSEN15")] <- 1
EntireCohort_pathvar$TSEN15[EntireCohort_pathvar$GENE5 %in% c("TSEN15")] <- 1
EntireCohort_pathvar$TSEN15[EntireCohort_pathvar$GENE6 %in% c("TSEN15")] <- 1
EntireCohort_pathvar$TSEN15[EntireCohort_pathvar$GENE7 %in% c("TSEN15")] <- 1
EntireCohort_pathvar$TSEN15[EntireCohort_pathvar$GENE8 %in% c("TSEN15")] <- 1
EntireCohort_pathvar$TSEN15[EntireCohort_pathvar$GENE9 %in% c("TSEN15")] <- 1

EntireCohort_pathvar$TSEN15[is.na(EntireCohort_pathvar$TSEN15)] <- 0
EntireCohort_pathvar$TSEN15 <- as.numeric(EntireCohort_pathvar$TSEN15)
unique(EntireCohort_pathvar$TSEN15)
# 0 1

table(EntireCohort_pathvar$TSEN15)
#  0      1 
# 469637    162 





### 'TSEN2', 
EntireCohort_pathvar$TSEN2 <- NA
EntireCohort_pathvar$TSEN2[EntireCohort_pathvar$GENE1 %in% c("TSEN2")] <- 1
EntireCohort_pathvar$TSEN2[EntireCohort_pathvar$GENE2 %in% c("TSEN2")] <- 1
EntireCohort_pathvar$TSEN2[EntireCohort_pathvar$GENE3 %in% c("TSEN2")] <- 1
EntireCohort_pathvar$TSEN2[EntireCohort_pathvar$GENE4 %in% c("TSEN2")] <- 1
EntireCohort_pathvar$TSEN2[EntireCohort_pathvar$GENE5 %in% c("TSEN2")] <- 1
EntireCohort_pathvar$TSEN2[EntireCohort_pathvar$GENE6 %in% c("TSEN2")] <- 1
EntireCohort_pathvar$TSEN2[EntireCohort_pathvar$GENE7 %in% c("TSEN2")] <- 1
EntireCohort_pathvar$TSEN2[EntireCohort_pathvar$GENE8 %in% c("TSEN2")] <- 1
EntireCohort_pathvar$TSEN2[EntireCohort_pathvar$GENE9 %in% c("TSEN2")] <- 1

EntireCohort_pathvar$TSEN2[is.na(EntireCohort_pathvar$TSEN2)] <- 0
EntireCohort_pathvar$TSEN2 <- as.numeric(EntireCohort_pathvar$TSEN2)
unique(EntireCohort_pathvar$TSEN2)
# 0 1

table(EntireCohort_pathvar$TSEN2)
#  0      1 
# 469061    738 




### 'TSEN54',  
EntireCohort_pathvar$TSEN54 <- NA
EntireCohort_pathvar$TSEN54[EntireCohort_pathvar$GENE1 %in% c("TSEN54")] <- 1
EntireCohort_pathvar$TSEN54[EntireCohort_pathvar$GENE2 %in% c("TSEN54")] <- 1
EntireCohort_pathvar$TSEN54[EntireCohort_pathvar$GENE3 %in% c("TSEN54")] <- 1
EntireCohort_pathvar$TSEN54[EntireCohort_pathvar$GENE4 %in% c("TSEN54")] <- 1
EntireCohort_pathvar$TSEN54[EntireCohort_pathvar$GENE5 %in% c("TSEN54")] <- 1
EntireCohort_pathvar$TSEN54[EntireCohort_pathvar$GENE6 %in% c("TSEN54")] <- 1
EntireCohort_pathvar$TSEN54[EntireCohort_pathvar$GENE7 %in% c("TSEN54")] <- 1
EntireCohort_pathvar$TSEN54[EntireCohort_pathvar$GENE8 %in% c("TSEN54")] <- 1
EntireCohort_pathvar$TSEN54[EntireCohort_pathvar$GENE9 %in% c("TSEN54")] <- 1

EntireCohort_pathvar$TSEN54[is.na(EntireCohort_pathvar$TSEN54)] <- 0
EntireCohort_pathvar$TSEN54 <- as.numeric(EntireCohort_pathvar$TSEN54)
unique(EntireCohort_pathvar$TSEN54)
# 0 1

table(EntireCohort_pathvar$TSEN54)
#  0      1 
# 469157    642 




### 'WDR4',         
EntireCohort_pathvar$WDR4 <- NA
EntireCohort_pathvar$WDR4[EntireCohort_pathvar$GENE1 %in% c("WDR4")] <- 1
EntireCohort_pathvar$WDR4[EntireCohort_pathvar$GENE2 %in% c("WDR4")] <- 1
EntireCohort_pathvar$WDR4[EntireCohort_pathvar$GENE3 %in% c("WDR4")] <- 1
EntireCohort_pathvar$WDR4[EntireCohort_pathvar$GENE4 %in% c("WDR4")] <- 1
EntireCohort_pathvar$WDR4[EntireCohort_pathvar$GENE5 %in% c("WDR4")] <- 1
EntireCohort_pathvar$WDR4[EntireCohort_pathvar$GENE6 %in% c("WDR4")] <- 1
EntireCohort_pathvar$WDR4[EntireCohort_pathvar$GENE7 %in% c("WDR4")] <- 1
EntireCohort_pathvar$WDR4[EntireCohort_pathvar$GENE8 %in% c("WDR4")] <- 1
EntireCohort_pathvar$WDR4[EntireCohort_pathvar$GENE9 %in% c("WDR4")] <- 1

EntireCohort_pathvar$WDR4[is.na(EntireCohort_pathvar$WDR4)] <- 0
EntireCohort_pathvar$WDR4 <- as.numeric(EntireCohort_pathvar$WDR4)
unique(EntireCohort_pathvar$WDR4)
# 0 1

table(EntireCohort_pathvar$WDR4)
# 0      1 
# 469485    314




### 'YARS1', 
# no variants present




### 'YRDC' 
EntireCohort_pathvar$YRDC <- NA
EntireCohort_pathvar$YRDC[EntireCohort_pathvar$GENE1 %in% c("YRDC")] <- 1
EntireCohort_pathvar$YRDC[EntireCohort_pathvar$GENE2 %in% c("YRDC")] <- 1
EntireCohort_pathvar$YRDC[EntireCohort_pathvar$GENE3 %in% c("YRDC")] <- 1
EntireCohort_pathvar$YRDC[EntireCohort_pathvar$GENE4 %in% c("YRDC")] <- 1
EntireCohort_pathvar$YRDC[EntireCohort_pathvar$GENE5 %in% c("YRDC")] <- 1
EntireCohort_pathvar$YRDC[EntireCohort_pathvar$GENE6 %in% c("YRDC")] <- 1
EntireCohort_pathvar$YRDC[EntireCohort_pathvar$GENE7 %in% c("YRDC")] <- 1
EntireCohort_pathvar$YRDC[EntireCohort_pathvar$GENE8 %in% c("YRDC")] <- 1
EntireCohort_pathvar$YRDC[EntireCohort_pathvar$GENE9 %in% c("YRDC")] <- 1

EntireCohort_pathvar$YRDC[is.na(EntireCohort_pathvar$YRDC)] <- 0
EntireCohort_pathvar$YRDC <- as.numeric(EntireCohort_pathvar$YRDC)
unique(EntireCohort_pathvar$YRDC)
# 0 1

table(EntireCohort_pathvar$YRDC)
# 0      1 
# 469349    450 




### 'RRP7A'  
EntireCohort_pathvar$RRP7A <- NA
EntireCohort_pathvar$RRP7A[EntireCohort_pathvar$GENE1 %in% c("RRP7A")] <- 1
EntireCohort_pathvar$RRP7A[EntireCohort_pathvar$GENE2 %in% c("RRP7A")] <- 1
EntireCohort_pathvar$RRP7A[EntireCohort_pathvar$GENE3 %in% c("RRP7A")] <- 1
EntireCohort_pathvar$RRP7A[EntireCohort_pathvar$GENE4 %in% c("RRP7A")] <- 1
EntireCohort_pathvar$RRP7A[EntireCohort_pathvar$GENE5 %in% c("RRP7A")] <- 1
EntireCohort_pathvar$RRP7A[EntireCohort_pathvar$GENE6 %in% c("RRP7A")] <- 1
EntireCohort_pathvar$RRP7A[EntireCohort_pathvar$GENE7 %in% c("RRP7A")] <- 1
EntireCohort_pathvar$RRP7A[EntireCohort_pathvar$GENE8 %in% c("RRP7A")] <- 1
EntireCohort_pathvar$RRP7A[EntireCohort_pathvar$GENE9 %in% c("RRP7A")] <- 1

EntireCohort_pathvar$RRP7A[is.na(EntireCohort_pathvar$RRP7A)] <- 0
EntireCohort_pathvar$RRP7A <- as.numeric(EntireCohort_pathvar$RRP7A)
unique(EntireCohort_pathvar$RRP7A)
# 0 1

table(EntireCohort_pathvar$RRP7A)
#  0      1 
# 468853    946 






# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)




# TFgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

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












# TGFbBMPgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

# go thr all the genes to check that they're present
grep("LTBP3", EntireCohort_pathvar$GENE)
# all seem to be present except for:
# BGN = 0
# FBN1 = 0


# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$TGFbBMPgene <- NA
EntireCohort_pathvar$TGFbBMPgene[EntireCohort_pathvar$GENE1 %in% c('SLC39A13', 'ADAMTSL2', 'SOX9', 'DDRGK1', 'LEMD3', 'BMP2', 'BMPR1B', 'GDF5', 'BMP1', 'LTBP3')] <- 1
EntireCohort_pathvar$TGFbBMPgene[EntireCohort_pathvar$GENE2 %in% c('SLC39A13', 'ADAMTSL2', 'SOX9', 'DDRGK1', 'LEMD3', 'BMP2', 'BMPR1B', 'GDF5', 'BMP1', 'LTBP3')] <- 1
EntireCohort_pathvar$TGFbBMPgene[EntireCohort_pathvar$GENE3 %in% c('SLC39A13', 'ADAMTSL2', 'SOX9', 'DDRGK1', 'LEMD3', 'BMP2', 'BMPR1B', 'GDF5', 'BMP1', 'LTBP3')] <- 1
EntireCohort_pathvar$TGFbBMPgene[EntireCohort_pathvar$GENE4 %in% c('SLC39A13', 'ADAMTSL2', 'SOX9', 'DDRGK1', 'LEMD3', 'BMP2', 'BMPR1B', 'GDF5', 'BMP1', 'LTBP3')] <- 1
EntireCohort_pathvar$TGFbBMPgene[EntireCohort_pathvar$GENE5 %in% c('SLC39A13', 'ADAMTSL2', 'SOX9', 'DDRGK1', 'LEMD3', 'BMP2', 'BMPR1B', 'GDF5', 'BMP1', 'LTBP3')] <- 1
EntireCohort_pathvar$TGFbBMPgene[EntireCohort_pathvar$GENE6 %in% c('SLC39A13', 'ADAMTSL2', 'SOX9', 'DDRGK1', 'LEMD3', 'BMP2', 'BMPR1B', 'GDF5', 'BMP1', 'LTBP3')] <- 1
EntireCohort_pathvar$TGFbBMPgene[EntireCohort_pathvar$GENE7 %in% c('SLC39A13', 'ADAMTSL2', 'SOX9', 'DDRGK1', 'LEMD3', 'BMP2', 'BMPR1B', 'GDF5', 'BMP1', 'LTBP3')] <- 1
EntireCohort_pathvar$TGFbBMPgene[EntireCohort_pathvar$GENE8 %in% c('SLC39A13', 'ADAMTSL2', 'SOX9', 'DDRGK1', 'LEMD3', 'BMP2', 'BMPR1B', 'GDF5', 'BMP1', 'LTBP3')] <- 1
EntireCohort_pathvar$TGFbBMPgene[EntireCohort_pathvar$GENE9 %in% c('SLC39A13', 'ADAMTSL2', 'SOX9', 'DDRGK1', 'LEMD3', 'BMP2', 'BMPR1B', 'GDF5', 'BMP1', 'LTBP3')] <- 1


EntireCohort_pathvar$TGFbBMPgene[is.na(EntireCohort_pathvar$TGFbBMPgene)] <- 0
EntireCohort_pathvar$TGFbBMPgene <- as.numeric(EntireCohort_pathvar$TGFbBMPgene)
unique(EntireCohort_pathvar$TGFbBMPgene)
# 0 1

table(EntireCohort_pathvar$TGFbBMPgene)
#     0      1 
# 456473  13326




# now to make a column for each of the TGFbBMPgenes
### 'SLC39A13'
EntireCohort_pathvar$SLC39A13 <- NA
EntireCohort_pathvar$SLC39A13[EntireCohort_pathvar$GENE1 %in% c("SLC39A13")] <- 1
EntireCohort_pathvar$SLC39A13[EntireCohort_pathvar$GENE2 %in% c("SLC39A13")] <- 1
EntireCohort_pathvar$SLC39A13[EntireCohort_pathvar$GENE3 %in% c("SLC39A13")] <- 1
EntireCohort_pathvar$SLC39A13[EntireCohort_pathvar$GENE4 %in% c("SLC39A13")] <- 1
EntireCohort_pathvar$SLC39A13[EntireCohort_pathvar$GENE5 %in% c("SLC39A13")] <- 1
EntireCohort_pathvar$SLC39A13[EntireCohort_pathvar$GENE6 %in% c("SLC39A13")] <- 1
EntireCohort_pathvar$SLC39A13[EntireCohort_pathvar$GENE7 %in% c("SLC39A13")] <- 1
EntireCohort_pathvar$SLC39A13[EntireCohort_pathvar$GENE8 %in% c("SLC39A13")] <- 1
EntireCohort_pathvar$SLC39A13[EntireCohort_pathvar$GENE9 %in% c("SLC39A13")] <- 1

EntireCohort_pathvar$SLC39A13[is.na(EntireCohort_pathvar$SLC39A13)] <- 0
EntireCohort_pathvar$SLC39A13 <- as.numeric(EntireCohort_pathvar$SLC39A13)
unique(EntireCohort_pathvar$SLC39A13)
# 0 1

table(EntireCohort_pathvar$SLC39A13)
#      0      1 
# 469425    374 





### 'ADAMTSL2'
EntireCohort_pathvar$ADAMTSL2 <- NA
EntireCohort_pathvar$ADAMTSL2[EntireCohort_pathvar$GENE1 %in% c("ADAMTSL2")] <- 1
EntireCohort_pathvar$ADAMTSL2[EntireCohort_pathvar$GENE2 %in% c("ADAMTSL2")] <- 1
EntireCohort_pathvar$ADAMTSL2[EntireCohort_pathvar$GENE3 %in% c("ADAMTSL2")] <- 1
EntireCohort_pathvar$ADAMTSL2[EntireCohort_pathvar$GENE4 %in% c("ADAMTSL2")] <- 1
EntireCohort_pathvar$ADAMTSL2[EntireCohort_pathvar$GENE5 %in% c("ADAMTSL2")] <- 1
EntireCohort_pathvar$ADAMTSL2[EntireCohort_pathvar$GENE6 %in% c("ADAMTSL2")] <- 1
EntireCohort_pathvar$ADAMTSL2[EntireCohort_pathvar$GENE7 %in% c("ADAMTSL2")] <- 1
EntireCohort_pathvar$ADAMTSL2[EntireCohort_pathvar$GENE8 %in% c("ADAMTSL2")] <- 1
EntireCohort_pathvar$ADAMTSL2[EntireCohort_pathvar$GENE9 %in% c("ADAMTSL2")] <- 1

EntireCohort_pathvar$ADAMTSL2[is.na(EntireCohort_pathvar$ADAMTSL2)] <- 0
EntireCohort_pathvar$ADAMTSL2 <- as.numeric(EntireCohort_pathvar$ADAMTSL2)
unique(EntireCohort_pathvar$ADAMTSL2)
# 0 1

table(EntireCohort_pathvar$ADAMTSL2)
#      0      1 
# 468078   1721




### 'SOX9'
EntireCohort_pathvar$SOX9 <- NA
EntireCohort_pathvar$SOX9[EntireCohort_pathvar$GENE1 %in% c("SOX9")] <- 1
EntireCohort_pathvar$SOX9[EntireCohort_pathvar$GENE2 %in% c("SOX9")] <- 1
EntireCohort_pathvar$SOX9[EntireCohort_pathvar$GENE3 %in% c("SOX9")] <- 1
EntireCohort_pathvar$SOX9[EntireCohort_pathvar$GENE4 %in% c("SOX9")] <- 1
EntireCohort_pathvar$SOX9[EntireCohort_pathvar$GENE5 %in% c("SOX9")] <- 1
EntireCohort_pathvar$SOX9[EntireCohort_pathvar$GENE6 %in% c("SOX9")] <- 1
EntireCohort_pathvar$SOX9[EntireCohort_pathvar$GENE7 %in% c("SOX9")] <- 1
EntireCohort_pathvar$SOX9[EntireCohort_pathvar$GENE8 %in% c("SOX9")] <- 1
EntireCohort_pathvar$SOX9[EntireCohort_pathvar$GENE9 %in% c("SOX9")] <- 1

EntireCohort_pathvar$SOX9[is.na(EntireCohort_pathvar$SOX9)] <- 0
EntireCohort_pathvar$SOX9 <- as.numeric(EntireCohort_pathvar$SOX9)
unique(EntireCohort_pathvar$SOX9)
# 0 1

table(EntireCohort_pathvar$SOX9)
#  0      1 
# 469398    401




### 'DDRGK1'
EntireCohort_pathvar$DDRGK1 <- NA
EntireCohort_pathvar$DDRGK1[EntireCohort_pathvar$GENE1 %in% c("DDRGK1")] <- 1
EntireCohort_pathvar$DDRGK1[EntireCohort_pathvar$GENE2 %in% c("DDRGK1")] <- 1
EntireCohort_pathvar$DDRGK1[EntireCohort_pathvar$GENE3 %in% c("DDRGK1")] <- 1
EntireCohort_pathvar$DDRGK1[EntireCohort_pathvar$GENE4 %in% c("DDRGK1")] <- 1
EntireCohort_pathvar$DDRGK1[EntireCohort_pathvar$GENE5 %in% c("DDRGK1")] <- 1
EntireCohort_pathvar$DDRGK1[EntireCohort_pathvar$GENE6 %in% c("DDRGK1")] <- 1
EntireCohort_pathvar$DDRGK1[EntireCohort_pathvar$GENE7 %in% c("DDRGK1")] <- 1
EntireCohort_pathvar$DDRGK1[EntireCohort_pathvar$GENE8 %in% c("DDRGK1")] <- 1
EntireCohort_pathvar$DDRGK1[EntireCohort_pathvar$GENE9 %in% c("DDRGK1")] <- 1

EntireCohort_pathvar$DDRGK1[is.na(EntireCohort_pathvar$DDRGK1)] <- 0
EntireCohort_pathvar$DDRGK1 <- as.numeric(EntireCohort_pathvar$DDRGK1)
unique(EntireCohort_pathvar$DDRGK1)
# 0 1

table(EntireCohort_pathvar$DDRGK1)
#   0      1 
# 469043    756




### 'BGN'
# no variants present




### 'LEMD3'
EntireCohort_pathvar$LEMD3 <- NA
EntireCohort_pathvar$LEMD3[EntireCohort_pathvar$GENE1 %in% c("LEMD3")] <- 1
EntireCohort_pathvar$LEMD3[EntireCohort_pathvar$GENE2 %in% c("LEMD3")] <- 1
EntireCohort_pathvar$LEMD3[EntireCohort_pathvar$GENE3 %in% c("LEMD3")] <- 1
EntireCohort_pathvar$LEMD3[EntireCohort_pathvar$GENE4 %in% c("LEMD3")] <- 1
EntireCohort_pathvar$LEMD3[EntireCohort_pathvar$GENE5 %in% c("LEMD3")] <- 1
EntireCohort_pathvar$LEMD3[EntireCohort_pathvar$GENE6 %in% c("LEMD3")] <- 1
EntireCohort_pathvar$LEMD3[EntireCohort_pathvar$GENE7 %in% c("LEMD3")] <- 1
EntireCohort_pathvar$LEMD3[EntireCohort_pathvar$GENE8 %in% c("LEMD3")] <- 1
EntireCohort_pathvar$LEMD3[EntireCohort_pathvar$GENE9 %in% c("LEMD3")] <- 1

EntireCohort_pathvar$LEMD3[is.na(EntireCohort_pathvar$LEMD3)] <- 0
EntireCohort_pathvar$LEMD3 <- as.numeric(EntireCohort_pathvar$LEMD3)
unique(EntireCohort_pathvar$LEMD3)
# 0 1

table(EntireCohort_pathvar$LEMD3)
# 0      1 
# 469215    584




### 'BMP2'
EntireCohort_pathvar$BMP2 <- NA
EntireCohort_pathvar$BMP2[EntireCohort_pathvar$GENE1 %in% c("BMP2")] <- 1
EntireCohort_pathvar$BMP2[EntireCohort_pathvar$GENE2 %in% c("BMP2")] <- 1
EntireCohort_pathvar$BMP2[EntireCohort_pathvar$GENE3 %in% c("BMP2")] <- 1
EntireCohort_pathvar$BMP2[EntireCohort_pathvar$GENE4 %in% c("BMP2")] <- 1
EntireCohort_pathvar$BMP2[EntireCohort_pathvar$GENE5 %in% c("BMP2")] <- 1
EntireCohort_pathvar$BMP2[EntireCohort_pathvar$GENE6 %in% c("BMP2")] <- 1
EntireCohort_pathvar$BMP2[EntireCohort_pathvar$GENE7 %in% c("BMP2")] <- 1
EntireCohort_pathvar$BMP2[EntireCohort_pathvar$GENE8 %in% c("BMP2")] <- 1
EntireCohort_pathvar$BMP2[EntireCohort_pathvar$GENE9 %in% c("BMP2")] <- 1

EntireCohort_pathvar$BMP2[is.na(EntireCohort_pathvar$BMP2)] <- 0
EntireCohort_pathvar$BMP2 <- as.numeric(EntireCohort_pathvar$BMP2)
unique(EntireCohort_pathvar$BMP2)
# 0 1

table(EntireCohort_pathvar$BMP2)
#  0      1 
# 469114    685 




### 'BMPR1B'
EntireCohort_pathvar$BMPR1B <- NA
EntireCohort_pathvar$BMPR1B[EntireCohort_pathvar$GENE1 %in% c("BMPR1B")] <- 1
EntireCohort_pathvar$BMPR1B[EntireCohort_pathvar$GENE2 %in% c("BMPR1B")] <- 1
EntireCohort_pathvar$BMPR1B[EntireCohort_pathvar$GENE3 %in% c("BMPR1B")] <- 1
EntireCohort_pathvar$BMPR1B[EntireCohort_pathvar$GENE4 %in% c("BMPR1B")] <- 1
EntireCohort_pathvar$BMPR1B[EntireCohort_pathvar$GENE5 %in% c("BMPR1B")] <- 1
EntireCohort_pathvar$BMPR1B[EntireCohort_pathvar$GENE6 %in% c("BMPR1B")] <- 1
EntireCohort_pathvar$BMPR1B[EntireCohort_pathvar$GENE7 %in% c("BMPR1B")] <- 1
EntireCohort_pathvar$BMPR1B[EntireCohort_pathvar$GENE8 %in% c("BMPR1B")] <- 1
EntireCohort_pathvar$BMPR1B[EntireCohort_pathvar$GENE9 %in% c("BMPR1B")] <- 1

EntireCohort_pathvar$BMPR1B[is.na(EntireCohort_pathvar$BMPR1B)] <- 0
EntireCohort_pathvar$BMPR1B <- as.numeric(EntireCohort_pathvar$BMPR1B)
unique(EntireCohort_pathvar$BMPR1B)
# 0 1

table(EntireCohort_pathvar$BMPR1B)
#    0      1 
# 467929   1870




### 'GDF5'
EntireCohort_pathvar$GDF5 <- NA
EntireCohort_pathvar$GDF5[EntireCohort_pathvar$GENE1 %in% c("GDF5")] <- 1
EntireCohort_pathvar$GDF5[EntireCohort_pathvar$GENE2 %in% c("GDF5")] <- 1
EntireCohort_pathvar$GDF5[EntireCohort_pathvar$GENE3 %in% c("GDF5")] <- 1
EntireCohort_pathvar$GDF5[EntireCohort_pathvar$GENE4 %in% c("GDF5")] <- 1
EntireCohort_pathvar$GDF5[EntireCohort_pathvar$GENE5 %in% c("GDF5")] <- 1
EntireCohort_pathvar$GDF5[EntireCohort_pathvar$GENE6 %in% c("GDF5")] <- 1
EntireCohort_pathvar$GDF5[EntireCohort_pathvar$GENE7 %in% c("GDF5")] <- 1
EntireCohort_pathvar$GDF5[EntireCohort_pathvar$GENE8 %in% c("GDF5")] <- 1
EntireCohort_pathvar$GDF5[EntireCohort_pathvar$GENE9 %in% c("GDF5")] <- 1

EntireCohort_pathvar$GDF5[is.na(EntireCohort_pathvar$GDF5)] <- 0
EntireCohort_pathvar$GDF5 <- as.numeric(EntireCohort_pathvar$GDF5)
unique(EntireCohort_pathvar$GDF5)
# 0 1

table(EntireCohort_pathvar$GDF5)
#  0      1 
# 469220    579





### 'BMP1'
EntireCohort_pathvar$BMP1 <- NA
EntireCohort_pathvar$BMP1[EntireCohort_pathvar$GENE1 %in% c("BMP1")] <- 1
EntireCohort_pathvar$BMP1[EntireCohort_pathvar$GENE2 %in% c("BMP1")] <- 1
EntireCohort_pathvar$BMP1[EntireCohort_pathvar$GENE3 %in% c("BMP1")] <- 1
EntireCohort_pathvar$BMP1[EntireCohort_pathvar$GENE4 %in% c("BMP1")] <- 1
EntireCohort_pathvar$BMP1[EntireCohort_pathvar$GENE5 %in% c("BMP1")] <- 1
EntireCohort_pathvar$BMP1[EntireCohort_pathvar$GENE6 %in% c("BMP1")] <- 1
EntireCohort_pathvar$BMP1[EntireCohort_pathvar$GENE7 %in% c("BMP1")] <- 1
EntireCohort_pathvar$BMP1[EntireCohort_pathvar$GENE8 %in% c("BMP1")] <- 1
EntireCohort_pathvar$BMP1[EntireCohort_pathvar$GENE9 %in% c("BMP1")] <- 1

EntireCohort_pathvar$BMP1[is.na(EntireCohort_pathvar$BMP1)] <- 0
EntireCohort_pathvar$BMP1 <- as.numeric(EntireCohort_pathvar$BMP1)
unique(EntireCohort_pathvar$BMP1)
# 0 1

table(EntireCohort_pathvar$BMP1)
#   0      1 
# 465572   4227




### 'FBN1'
# no variants present




### 'LTBP3'
EntireCohort_pathvar$LTBP3 <- NA
EntireCohort_pathvar$LTBP3[EntireCohort_pathvar$GENE1 %in% c("LTBP3")] <- 1
EntireCohort_pathvar$LTBP3[EntireCohort_pathvar$GENE2 %in% c("LTBP3")] <- 1
EntireCohort_pathvar$LTBP3[EntireCohort_pathvar$GENE3 %in% c("LTBP3")] <- 1
EntireCohort_pathvar$LTBP3[EntireCohort_pathvar$GENE4 %in% c("LTBP3")] <- 1
EntireCohort_pathvar$LTBP3[EntireCohort_pathvar$GENE5 %in% c("LTBP3")] <- 1
EntireCohort_pathvar$LTBP3[EntireCohort_pathvar$GENE6 %in% c("LTBP3")] <- 1
EntireCohort_pathvar$LTBP3[EntireCohort_pathvar$GENE7 %in% c("LTBP3")] <- 1
EntireCohort_pathvar$LTBP3[EntireCohort_pathvar$GENE8 %in% c("LTBP3")] <- 1
EntireCohort_pathvar$LTBP3[EntireCohort_pathvar$GENE9 %in% c("LTBP3")] <- 1

EntireCohort_pathvar$LTBP3[is.na(EntireCohort_pathvar$LTBP3)] <- 0
EntireCohort_pathvar$LTBP3 <- as.numeric(EntireCohort_pathvar$LTBP3)
unique(EntireCohort_pathvar$LTBP3)
# 0 1

table(EntireCohort_pathvar$LTBP3)
#  0      1 
# 467488   2311








# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)











# THgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

# go thr all the genes to check that they're present
grep("GNAS", EntireCohort_pathvar$GENE)
# all seem to be present except for these genes
# PRKAR1A = 0
# THRA = 0



# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$THgene <- NA
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE1 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS')] <- 1
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE2 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS')] <- 1
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE3 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS')] <- 1
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE4 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS')] <- 1
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE5 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS')] <- 1
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE6 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS')] <- 1
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE7 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS')] <- 1
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE8 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS')] <- 1
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE9 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS')] <- 1


EntireCohort_pathvar$THgene[is.na(EntireCohort_pathvar$THgene)] <- 0
EntireCohort_pathvar$THgene <- as.numeric(EntireCohort_pathvar$THgene)
unique(EntireCohort_pathvar$THgene)
# 0 1 

table(EntireCohort_pathvar$THgene)
#     0      1 
# 462292   7507



# now to make a column for each of the THgenes
### 'PRKAR1A' = no variants present

### 'SIK3'
EntireCohort_pathvar$SIK3 <- NA
EntireCohort_pathvar$SIK3[EntireCohort_pathvar$GENE1 %in% c("SIK3")] <- 1
EntireCohort_pathvar$SIK3[EntireCohort_pathvar$GENE2 %in% c("SIK3")] <- 1
EntireCohort_pathvar$SIK3[EntireCohort_pathvar$GENE3 %in% c("SIK3")] <- 1
EntireCohort_pathvar$SIK3[EntireCohort_pathvar$GENE4 %in% c("SIK3")] <- 1
EntireCohort_pathvar$SIK3[EntireCohort_pathvar$GENE5 %in% c("SIK3")] <- 1
EntireCohort_pathvar$SIK3[EntireCohort_pathvar$GENE6 %in% c("SIK3")] <- 1
EntireCohort_pathvar$SIK3[EntireCohort_pathvar$GENE7 %in% c("SIK3")] <- 1
EntireCohort_pathvar$SIK3[EntireCohort_pathvar$GENE8 %in% c("SIK3")] <- 1
EntireCohort_pathvar$SIK3[EntireCohort_pathvar$GENE9 %in% c("SIK3")] <- 1

EntireCohort_pathvar$SIK3[is.na(EntireCohort_pathvar$SIK3)] <- 0
EntireCohort_pathvar$SIK3 <- as.numeric(EntireCohort_pathvar$SIK3)
unique(EntireCohort_pathvar$SIK3)
# 0 1 

table(EntireCohort_pathvar$SIK3)
# 0      1 
# 467073   2726




### 'TBCE'
EntireCohort_pathvar$TBCE <- NA
EntireCohort_pathvar$TBCE[EntireCohort_pathvar$GENE1 %in% c("TBCE")] <- 1
EntireCohort_pathvar$TBCE[EntireCohort_pathvar$GENE2 %in% c("TBCE")] <- 1
EntireCohort_pathvar$TBCE[EntireCohort_pathvar$GENE3 %in% c("TBCE")] <- 1
EntireCohort_pathvar$TBCE[EntireCohort_pathvar$GENE4 %in% c("TBCE")] <- 1
EntireCohort_pathvar$TBCE[EntireCohort_pathvar$GENE5 %in% c("TBCE")] <- 1
EntireCohort_pathvar$TBCE[EntireCohort_pathvar$GENE6 %in% c("TBCE")] <- 1
EntireCohort_pathvar$TBCE[EntireCohort_pathvar$GENE7 %in% c("TBCE")] <- 1
EntireCohort_pathvar$TBCE[EntireCohort_pathvar$GENE8 %in% c("TBCE")] <- 1
EntireCohort_pathvar$TBCE[EntireCohort_pathvar$GENE9 %in% c("TBCE")] <- 1

EntireCohort_pathvar$TBCE[is.na(EntireCohort_pathvar$TBCE)] <- 0
EntireCohort_pathvar$TBCE <- as.numeric(EntireCohort_pathvar$TBCE)
unique(EntireCohort_pathvar$TBCE)
# 0 1

table(EntireCohort_pathvar$TBCE)
# 0      1 
# 469224    575




### 'TRIP11'
EntireCohort_pathvar$TRIP11 <- NA
EntireCohort_pathvar$TRIP11[EntireCohort_pathvar$GENE1 %in% c("TRIP11")] <- 1
EntireCohort_pathvar$TRIP11[EntireCohort_pathvar$GENE2 %in% c("TRIP11")] <- 1
EntireCohort_pathvar$TRIP11[EntireCohort_pathvar$GENE3 %in% c("TRIP11")] <- 1
EntireCohort_pathvar$TRIP11[EntireCohort_pathvar$GENE4 %in% c("TRIP11")] <- 1
EntireCohort_pathvar$TRIP11[EntireCohort_pathvar$GENE5 %in% c("TRIP11")] <- 1
EntireCohort_pathvar$TRIP11[EntireCohort_pathvar$GENE6 %in% c("TRIP11")] <- 1
EntireCohort_pathvar$TRIP11[EntireCohort_pathvar$GENE7 %in% c("TRIP11")] <- 1
EntireCohort_pathvar$TRIP11[EntireCohort_pathvar$GENE8 %in% c("TRIP11")] <- 1
EntireCohort_pathvar$TRIP11[EntireCohort_pathvar$GENE9 %in% c("TRIP11")] <- 1

EntireCohort_pathvar$TRIP11[is.na(EntireCohort_pathvar$TRIP11)] <- 0
EntireCohort_pathvar$TRIP11 <- as.numeric(EntireCohort_pathvar$TRIP11)
unique(EntireCohort_pathvar$TRIP11)
# 0 1

table(EntireCohort_pathvar$TRIP11)
#  0      1 
# 467791   2008 




### 'PTHLH'
EntireCohort_pathvar$PTHLH <- NA
EntireCohort_pathvar$PTHLH[EntireCohort_pathvar$GENE1 %in% c("PTHLH")] <- 1
EntireCohort_pathvar$PTHLH[EntireCohort_pathvar$GENE2 %in% c("PTHLH")] <- 1
EntireCohort_pathvar$PTHLH[EntireCohort_pathvar$GENE3 %in% c("PTHLH")] <- 1
EntireCohort_pathvar$PTHLH[EntireCohort_pathvar$GENE4 %in% c("PTHLH")] <- 1
EntireCohort_pathvar$PTHLH[EntireCohort_pathvar$GENE5 %in% c("PTHLH")] <- 1
EntireCohort_pathvar$PTHLH[EntireCohort_pathvar$GENE6 %in% c("PTHLH")] <- 1
EntireCohort_pathvar$PTHLH[EntireCohort_pathvar$GENE7 %in% c("PTHLH")] <- 1
EntireCohort_pathvar$PTHLH[EntireCohort_pathvar$GENE8 %in% c("PTHLH")] <- 1
EntireCohort_pathvar$PTHLH[EntireCohort_pathvar$GENE9 %in% c("PTHLH")] <- 1

EntireCohort_pathvar$PTHLH[is.na(EntireCohort_pathvar$PTHLH)] <- 0
EntireCohort_pathvar$PTHLH <- as.numeric(EntireCohort_pathvar$PTHLH)
unique(EntireCohort_pathvar$PTHLH)
# 0 1

table(EntireCohort_pathvar$PTHLH)
# 0      1 
# 469602    197




### 'PTH1R'
EntireCohort_pathvar$PTH1R <- NA
EntireCohort_pathvar$PTH1R[EntireCohort_pathvar$GENE1 %in% c("PTH1R")] <- 1
EntireCohort_pathvar$PTH1R[EntireCohort_pathvar$GENE2 %in% c("PTH1R")] <- 1
EntireCohort_pathvar$PTH1R[EntireCohort_pathvar$GENE3 %in% c("PTH1R")] <- 1
EntireCohort_pathvar$PTH1R[EntireCohort_pathvar$GENE4 %in% c("PTH1R")] <- 1
EntireCohort_pathvar$PTH1R[EntireCohort_pathvar$GENE5 %in% c("PTH1R")] <- 1
EntireCohort_pathvar$PTH1R[EntireCohort_pathvar$GENE6 %in% c("PTH1R")] <- 1
EntireCohort_pathvar$PTH1R[EntireCohort_pathvar$GENE7 %in% c("PTH1R")] <- 1
EntireCohort_pathvar$PTH1R[EntireCohort_pathvar$GENE8 %in% c("PTH1R")] <- 1
EntireCohort_pathvar$PTH1R[EntireCohort_pathvar$GENE9 %in% c("PTH1R")] <- 1

EntireCohort_pathvar$PTH1R[is.na(EntireCohort_pathvar$PTH1R)] <- 0
EntireCohort_pathvar$PTH1R <- as.numeric(EntireCohort_pathvar$PTH1R)
unique(EntireCohort_pathvar$PTH1R)
# 0 1

table(EntireCohort_pathvar$PTH1R)
#  0      1 
# 469187    612


### 'THRA' = no variants present


### 'GNAS'
EntireCohort_pathvar$GNAS <- NA
EntireCohort_pathvar$GNAS[EntireCohort_pathvar$GENE1 %in% c("GNAS")] <- 1
EntireCohort_pathvar$GNAS[EntireCohort_pathvar$GENE2 %in% c("GNAS")] <- 1
EntireCohort_pathvar$GNAS[EntireCohort_pathvar$GENE3 %in% c("GNAS")] <- 1
EntireCohort_pathvar$GNAS[EntireCohort_pathvar$GENE4 %in% c("GNAS")] <- 1
EntireCohort_pathvar$GNAS[EntireCohort_pathvar$GENE5 %in% c("GNAS")] <- 1
EntireCohort_pathvar$GNAS[EntireCohort_pathvar$GENE6 %in% c("GNAS")] <- 1
EntireCohort_pathvar$GNAS[EntireCohort_pathvar$GENE7 %in% c("GNAS")] <- 1
EntireCohort_pathvar$GNAS[EntireCohort_pathvar$GENE8 %in% c("GNAS")] <- 1
EntireCohort_pathvar$GNAS[EntireCohort_pathvar$GENE9 %in% c("GNAS")] <- 1

EntireCohort_pathvar$GNAS[is.na(EntireCohort_pathvar$GNAS)] <- 0
EntireCohort_pathvar$GNAS <- as.numeric(EntireCohort_pathvar$GNAS)
unique(EntireCohort_pathvar$GNAS)
# 0 1

table(EntireCohort_pathvar$GNAS)
# 0      1 
# 468358   1441




# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)




# WNTSHHgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

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




# checking for SS/M diagnoses --------------------------------------------------

# Q02 Microcephaly
# 
# Q04 Other congenital malformations of brain
# 
# Q77 Osteochondrodysplasia with defects of growth of tubular bones and spine
# 
# Q78 Other osteochondrodysplasias
# 
# Q79.6 Ehlers-Danlos syndrome
# Q79.8 Other congenital malformations of musculoskeletal system
# Q79.9 Congenital malformation of musculoskeletal system, unspecified
# 
# Q87.1 Congenital malformation syndromes predominantly associated with short stature
# Q87.2 Congenital malformation syndromes predominantly involving limbs
# 
# Q96 Turner's syndrome

# are any of these people carriers of vars

EntireCohort_Diagnoses_ICD10$eid <- as.character(EntireCohort_Diagnoses_ICD10$eid)




# CAgene [non-LoF vars present]
CAgene_1 <- EntireCohort_pathvar[EntireCohort_pathvar$CAgene == "1",]
test1 <- merge(EntireCohort_Diagnoses_ICD10, CAgene_1, by = c("eid", "eid"))
#  40 obs

test1$eid
# [1] 1319728 1354999 1395783 1666711 1709285 1772972 1857093 1865458 1879176 1920283 1987949 2115376 2519776 2626422 2654573 2805253 2847418
# [18] 3079809 3135104 3456772 3597233 3634548 3757282 3922237 4027505 4312376 4449879 4473170 4621506 4627044 4666537 4837383 4924862 4987061
# [35] 5274452 5332220 5413923 5509413 5761678 5938652
test1$p41270
## relevant:
# [4] Q87.1 Congenital malformation syndromes predominantly associated with short stature
# [9] Q87.1 Congenital malformation syndromes predominantly associated with short stature
# [18] Q87.1 Congenital malformation syndromes predominantly associated with short stature
# [29] Q87.1 Congenital malformation syndromes predominantly associated with short stature

## relevant eids:
# [4] 1666711 = PTPN11 (GOF/DN dom), 12:112450362:A:G_G (0 cancer)
# [9] 1879176 = SOS1 (GOF dom), 2:39051202:A:G_G. RECQL4 (LOF rec), 8:144514449:A:G_G (0 cancer)
# [18] 3079809 = PTPN11 (GOF/DN dom), 12:112450408:G:C_C (cns cancer)
# [29] 4621506 = SOS1 (GOF dom), 2:39022773:C:T_T (0 cancer)






# check vars associated with GoF/DN ---------------------------------------

#  BRAF, CBL, KRAS, NRAS, PTPN11, RIT1, RRAS2, SMAD4, SOS1, and THRB

test1 <- EntireCohort_pathvar[EntireCohort_pathvar$THRB == 1,]
test1$GENE
test1$VAR

## BRAF 
# 7:140801550:G:A_A

## CBL 
# 11:119278541:G:A_A = x15
# 11:119278181:T:A_A = x1
# 11:119278181:T:C_C = x3                                                                                  

## KRAS
# 12:25245372:T:C_C = x2
# 12:25227345:C:A_A = x2
# 12:25209904:T:C_C = x1


## NRAS 
# 1:114716126:C:T_T = x12
# 1:114716127:C:T_T = x1


## PTPN11
# 12:112455968:A:G_G = x6
# 12:112450346:A:G_G = x4
# 12:112477719:A:G_G = x5
# 12:112450368:A:G_G = x5    
# 12:112489084:G:A_A = x1
# 12:112472981:G:A_A = x10
# 12:112489105:A:G_G = x11
# 12:112472961:G:T_T = x1
# 12:112489068:C:T_T = x5
# 12:112450389:A:G_G = x2
# 12:112453279:G:C_C = x1
# 12:112473031:A:G_G = x1
# 12:112489048:C:A_A = x1
# 12:112488466:C:T_T = x6
# 12:112450361:G:A_A = x1
# 12:112489048:C:T_T = x2
# 12:112450362:A:G_G = x1
# 12:112450408:G:C_C = x2                                                
# 12:112489069:G:T_T = x1
# 12:112419116:C:T_T = x3
# 12:112450416:A:G_G = x2              
# 12:112489047:C:A_A = x1                                                      
# 12:112489047:C:T_T = x2
# 12:112472968:C:T_T = x3
# 12:112489086:A:G_G = x2
# 12:112472989:G:T_T = x2
# 12:112489078:G:A_A = x2
# 12:112450359:G:C_C = x1
# 12:112472954:A:G_G = x1


## RIT1
# 1:155904738:G:C_C = x1
# 1:155904475:A:G_G = x1
# 1:155904739:C:A_A = x1
# 1:155904489:G:A_A = x2
# 1:155904494:A:T_T = x1
# 1:155904495:A:G_G = x1


## RRAS2
# 11:14294851:C:T_T = x1

## SMAD4
# 18:51078294:C:T_T = x4


## SOS1
# 2:39022779:A:G_G = x3
# 2:39022772:C:G_G = x1
# 2:39013523:A:G_G = x5
# 2:39051202:A:G_G = x2
# 2:39022786:T:G_G = x1
# 2:39022772:C:A_A = x1
# 2:39022774:T:C_C = x2
# 2:39023118:A:G_G = x1
# 2:39022773:C:T_T = x1
# 2:39023106:C:T_T = x1


## THRB
# 3:24127685:G:A_A = x1
# 3:24122892:C:T_T = x2
# 3:24122913:G:C_C = x2
# 3:24127684:C:T_T = x2
# 3:24123122:C:T_T = x1
# 3:24143512:G:A_A = x1
