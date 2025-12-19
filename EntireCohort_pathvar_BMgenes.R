# 24/06/2024

# EntireCohort_pathvar_BMgenes.R

# this script is for the Bone matrix associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)


# BMgene
# IFITM5
# FN1
# ACAN
# ADAMTS2
# B3GALT6
# B3GAT3
# B4GALT7
# CHST3
# COL10A1
# COL11A1
# COL11A2
# COL1A1
# COL1A2
# COL27A1
# COL2A1
# COL9A3
# CRTAP
# FAM20C
# FKBP10
# KDELR2
# MATN3
# MBTPS2
# MMP13
# P3H1
# SERPINH1
# SLC10A7
# SLC26A2
# SLC35D1
# CCN6
# SERPINF1









# BMgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')


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

















# BMgene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_BMgene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_BMgene_, "pathvar_BMgene_.csv", row.names=FALSE)
# system('dx upload pathvar_BMgene_.csv --path Emily-folder/results_2/pathvar_BMgene_.csv')


# BMgene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ BMgene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_BMgene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_incidence_1, "pathvar_BMgene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_incidence_1.csv --path Emily-folder/results_2/pathvar_BMgene_incidence_1.csv')


# GOIs
LR_2 <- glm(formula = cancer_ROOC_dummy ~ FN1+ACAN+ADAMTS2+B3GALT6+B3GAT3+B4GALT7+CHST3+COL10A1+COL11A1+COL11A2+COL1A1+COL27A1+COL2A1+CRTAP+FAM20C+FKBP10+KDELR2+MATN3+MMP13+P3H1+SERPINH1+SLC10A7+SLC26A2+SLC35D1+CCN6+SERPINF1, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_BMgene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_incidence_1, "pathvar_BMgene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_incidence_2.csv --path Emily-folder/results_2/pathvar_BMgene_incidence_2.csv')







### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMgene", 'FN1', 'ACAN', 'ADAMTS2', 
                                     'B3GALT6', 'B3GAT3', 'B4GALT7', 'CHST3', 
                                     'COL10A1', 'COL11A1', 'COL11A2', 'COL1A1', 
                                     'COL27A1', 'COL2A1', 
                                     'CRTAP', 'FAM20C', 'FKBP10', 'KDELR2', 'MATN3', 
                                     'MMP13', 'P3H1', 'SERPINH1', 'SLC10A7', 
                                     'SLC26A2', 'SLC35D1', 'CCN6', 'SERPINF1')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     33

with(data_setA, table(SERPINF1))
# looking good


# BMgene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ BMgene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_BMgene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_incidence_3, "pathvar_BMgene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_incidence_3.csv --path Emily-folder/results_2/pathvar_BMgene_incidence_3.csv')


# age + ethnicity + smoking + sex + BMgene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_BMgene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_incidence_4, "pathvar_BMgene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_incidence_4.csv --path Emily-folder/results_2/pathvar_BMgene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ FN1+ACAN+ADAMTS2+B3GALT6+B3GAT3+B4GALT7+CHST3+COL10A1+COL11A1+COL11A2+COL1A1+COL27A1+COL2A1+CRTAP+FAM20C+FKBP10+KDELR2+MATN3+MMP13+P3H1+SERPINH1+SLC10A7+SLC26A2+SLC35D1+CCN6+SERPINF1, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_BMgene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_incidence_5, "pathvar_BMgene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_incidence_5.csv --path Emily-folder/results_2/pathvar_BMgene_incidence_5.csv')


# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FN1+ACAN+ADAMTS2+B3GALT6+B3GAT3+B4GALT7+CHST3+COL10A1+COL11A1+COL11A2+COL1A1+COL27A1+COL2A1+CRTAP+FAM20C+FKBP10+KDELR2+MATN3+MMP13+P3H1+SERPINH1+SLC10A7+SLC26A2+SLC35D1+CCN6+SERPINF1, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_BMgene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_incidence_6, "pathvar_BMgene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_incidence_6.csv --path Emily-folder/results_2/pathvar_BMgene_incidence_6.csv')




### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "BMgene", 'FN1', 'ACAN', 'ADAMTS2', 
                                     'B3GALT6', 'B3GAT3', 'B4GALT7', 'CHST3', 'COL10A1', 
                                     'COL11A1', 'COL11A2', 'COL1A1', 'COL27A1', 
                                     'COL2A1', 'CRTAP', 'FAM20C', 'FKBP10', 
                                     'KDELR2', 'MATN3', 'MMP13', 'P3H1', 
                                     'SERPINH1', 'SLC10A7', 'SLC26A2', 'SLC35D1', 
                                     'CCN6', 'SERPINF1')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     38

with(data_setB, table(SERPINF1))



# BMgene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ BMgene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_BMgene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_incidence_7, "pathvar_BMgene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_incidence_7.csv --path Emily-folder/results_2/pathvar_BMgene_incidence_7.csv')


# all variables + BMgene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+BMgene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_BMgene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_incidence_8, "pathvar_BMgene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_incidence_8.csv --path Emily-folder/results_2/pathvar_BMgene_incidence_8.csv')







# GOIs alone
LR_9 <- glm(formula = cancer_ROOC_dummy ~ FN1+ACAN+ADAMTS2+B3GALT6+B3GAT3+B4GALT7+CHST3+COL10A1+COL11A1+COL11A2+COL1A1+COL27A1+COL2A1+CRTAP+FAM20C+FKBP10+KDELR2+MATN3+MMP13+P3H1+SERPINH1+SLC10A7+SLC26A2+SLC35D1+CCN6+SERPINF1, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_BMgene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_incidence_9, "pathvar_BMgene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_incidence_9.csv --path Emily-folder/results_2/pathvar_BMgene_incidence_9.csv')


# all variables + GOIs
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+FN1+ACAN+ADAMTS2+B3GALT6+B3GAT3+B4GALT7+CHST3+COL10A1+COL11A1+COL11A2+COL1A1+COL27A1+COL2A1+CRTAP+FAM20C+FKBP10+KDELR2+MATN3+MMP13+P3H1+SERPINH1+SLC10A7+SLC26A2+SLC35D1+CCN6+SERPINF1, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_BMgene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_incidence_10, "pathvar_BMgene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_incidence_10.csv --path Emily-folder/results_2/pathvar_BMgene_incidence_10.csv')













### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(BMgene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(BMgene, cancer_ROOC_dummy)))
# p-value = 0.1604

with(EntireCohort_pathvar, table(BMgene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(BMgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.7138746 0.7595160
# sample estimates:
#   odds ratio 
# 0.7362953

with(EntireCohort_pathvar, table(BMgene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(BMgene, Ever_smoked_dummy)))
# p-value = 0.07018

with(EntireCohort_pathvar, table(BMgene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(BMgene, Sex_Male)))
# p-value = 0.7395

with(EntireCohort_pathvar, table(BMgene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(BMgene, University_education)))
# p-value = 0.006784
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.006083 1.038783
# sample estimates:
#   odds ratio 
# 1.022331

with(EntireCohort_pathvar, table(BMgene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(BMgene, Annual_income))	)
# p-value = 0.5587

with(EntireCohort_pathvar, table(BMgene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(BMgene, Smoker_previous)))
# p-value = 0.007295
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9633061 0.9942094
# sample estimates:
#   odds ratio 
# 0.9786426

with(EntireCohort_pathvar, table(BMgene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(BMgene, Smoker_current))	)
# p-value = 0.3026



### data_setA
with(data_setA, table(BMgene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(BMgene, cancer_ROOC_dummy)))
# p-value = 0.1573

with(data_setA, table(BMgene, Ethnicity_White))			
fisher.test(with(data_setA, table(BMgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.7129983 0.7587488
# sample estimates:
#   odds ratio 
# 0.7354814

with(data_setA, table(BMgene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(BMgene, Ever_smoked_dummy)))
# p-value = 0.0819

with(data_setA, table(BMgene, Sex_Male))			
fisher.test(with(data_setA, table(BMgene, Sex_Male)))
# p-value = 0.8119



### data_setB
with(data_setB, table(BMgene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(BMgene, cancer_ROOC_dummy)))
# p-value = 0.05542

with(data_setB, table(BMgene, Ethnicity_White))			
fisher.test(with(data_setB, table(BMgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.7162958 0.7720305
# sample estimates:
#   odds ratio 
# 0.7435699

with(data_setB, table(BMgene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(BMgene, Ever_smoked_dummy)))
# p-value = 0.6277

with(data_setB, table(BMgene, Sex_Male))			
fisher.test(with(data_setB, table(BMgene, Sex_Male)))
# p-value = 0.5086

with(data_setB, table(BMgene, University_education))			
fisher.test(with(data_setB, table(BMgene, University_education)))
# p-value = 0.001744
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.010331 1.045937
# sample estimates:
#   odds ratio 
# 1.028023 

with(data_setB, table(BMgene, Annual_income))			
fisher.test(with(data_setB, table(BMgene, Annual_income)))
# p-value = 0.7625






### can plot the continuous variables: 

# make BMgene into factor
df_for_plot <- EntireCohort_pathvar
df_for_plot$BMgene <- as.factor(df_for_plot$BMgene)

# try a violin plot

# Age_at_recruitment  
df_for_plot %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ BMgene, data = df_for_plot)
# W = 1.6308e+10, p-value = 0.001587

# BMI  
df_for_plot %>% 
  ggplot(mapping = aes(x = BMgene, y = BMI)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = BMgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ BMgene, data = df_for_plot)
# W = 1.6113e+10, p-value = 0.1307

# Standing_height 
df_for_plot %>% 
  ggplot(mapping = aes(x = BMgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = BMgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ BMgene, data = df_for_plot)
# W = 1.6296e+10, p-value = 1.528e-08

# make this a bar plot
# Weekly_exercise
df_for_plot %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(BMgene), scales = "free_y") +
  theme_light()
# 



### data_setA

# make BMgene into factor
data_setA$BMgene <- as.factor(data_setA$BMgene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ BMgene, data = data_setA)
# W = 1.6032e+10, p-value = 0.001932



### data_setB

# make BMgene into factor
data_setB$BMgene <- as.factor(data_setB$BMgene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ BMgene, data = data_setB)
# W = 1.0718e+10, p-value = 0.0007529

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = BMgene, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = BMgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ BMgene, data = data_setB)
# W = 1.0659e+10, p-value = 0.2731

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = BMgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = BMgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ BMgene, data = data_setB)
# W = 1.0724e+10, p-value = 0.0002939


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(BMgene), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# BMgene Type ------------------------------------------------------------

# make sure the execute the code in the 'BMgenes' section

## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_BMgene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_type_1, "pathvar_BMgene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_type_1.csv --path Emily-folder/results_2/pathvar_BMgene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_BMgene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_type_2, "pathvar_BMgene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_type_2.csv --path Emily-folder/results_2/pathvar_BMgene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_BMgene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_type_3, "pathvar_BMgene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_type_3.csv --path Emily-folder/results_2/pathvar_BMgene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_BMgene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_type_4, "pathvar_BMgene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_type_4.csv --path Emily-folder/results_2/pathvar_BMgene_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_BMgene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_type_5, "pathvar_BMgene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_type_5.csv --path Emily-folder/results_2/pathvar_BMgene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_BMgene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_type_6, "pathvar_BMgene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_type_6.csv --path Emily-folder/results_2/pathvar_BMgene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_BMgene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_type_7, "pathvar_BMgene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_type_7.csv --path Emily-folder/results_2/pathvar_BMgene_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_BMgene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_type_8, "pathvar_BMgene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_type_8.csv --path Emily-folder/results_2/pathvar_BMgene_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_BMgene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_type_9, "pathvar_BMgene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_type_9.csv --path Emily-folder/results_2/pathvar_BMgene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_BMgene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_type_10, "pathvar_BMgene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_type_10.csv --path Emily-folder/results_2/pathvar_BMgene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_BMgene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_type_11, "pathvar_BMgene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_type_11.csv --path Emily-folder/results_2/pathvar_BMgene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_BMgene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_type_12, "pathvar_BMgene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_type_12.csv --path Emily-folder/results_2/pathvar_BMgene_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_BMgene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_type_13, "pathvar_BMgene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_type_13.csv --path Emily-folder/results_2/pathvar_BMgene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_BMgene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_type_14, "pathvar_BMgene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_type_14.csv --path Emily-folder/results_2/pathvar_BMgene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_BMgene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_BMgene_type_15, "pathvar_BMgene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_type_15.csv --path Emily-folder/results_2/pathvar_BMgene_type_15.csv')






# BMgene Age -------------------------------------------------------------

# make sure the execute the code in the 'BMgenes' section

### linear regression
# "BMgene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ BMgene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_BMgene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_BMgene_age_1, "pathvar_BMgene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_age_1.csv --path Emily-folder/results_2/pathvar_BMgene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ FN1+ACAN+ADAMTS2+B3GALT6+B3GAT3+B4GALT7+CHST3+COL10A1+COL11A1+COL11A2+COL1A1+COL27A1+COL2A1+CRTAP+FAM20C+FKBP10+KDELR2+MATN3+MMP13+P3H1+SERPINH1+SLC10A7+SLC26A2+SLC35D1+CCN6+SERPINF1, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_BMgene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_BMgene_age_2, "pathvar_BMgene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_age_2.csv --path Emily-folder/results_2/pathvar_BMgene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "BMgene", 'FN1', 'ACAN', 'ADAMTS2', 'B3GALT6', 
                                     'B3GAT3', 'B4GALT7', 'CHST3', 'COL10A1', 'COL11A1', 
                                     'COL11A2', 'COL1A1', 'COL27A1', 'COL2A1', 
                                     'CRTAP', 'FAM20C', 'FKBP10', 'KDELR2', 
                                     'MATN3', 'MMP13', 'P3H1', 'SERPINH1', 
                                     'SLC10A7', 'SLC26A2', 'SLC35D1', 'CCN6', 'SERPINF1')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     33


with(data_setD, table(SERPINF1))
# looks good



data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "BMgene", 'FN1', 'ACAN', 
                                     'ADAMTS2', 'B3GALT6', 'B3GAT3', 'B4GALT7', 'CHST3', 
                                     'COL10A1', 'COL11A1', 'COL11A2', 'COL1A1', 
                                     'COL27A1', 'COL2A1', 'CRTAP', 'FAM20C', 
                                     'FKBP10', 'KDELR2', 'MATN3', 'MMP13', 
                                     'P3H1', 'SERPINH1', 'SLC10A7', 'SLC26A2', 'SLC35D1', 
                                     'CCN6', 'SERPINF1')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    38

with(data_setE, table(SERPINF1))
# looks good



### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ BMgene, data = data_setD)
# summary(LM_15)
pathvar_BMgene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_BMgene_age_3, "pathvar_BMgene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_age_3.csv --path Emily-folder/results_2/pathvar_BMgene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMgene, data = data_setD)
# summary(LM_16)
pathvar_BMgene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_BMgene_age_4, "pathvar_BMgene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_age_4.csv --path Emily-folder/results_2/pathvar_BMgene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ FN1+ACAN+ADAMTS2+B3GALT6+B3GAT3+B4GALT7+CHST3+COL10A1+COL11A1+COL11A2+COL1A1+COL27A1+COL2A1+CRTAP+FAM20C+FKBP10+KDELR2+MATN3+MMP13+P3H1+SERPINH1+SLC10A7+SLC26A2+SLC35D1+CCN6+SERPINF1, data = data_setD)
# summary(LM_17)
pathvar_BMgene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_BMgene_age_5, "pathvar_BMgene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_age_5.csv --path Emily-folder/results_2/pathvar_BMgene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+FN1+ACAN+ADAMTS2+B3GALT6+B3GAT3+B4GALT7+CHST3+COL10A1+COL11A1+COL11A2+COL1A1+COL27A1+COL2A1+CRTAP+FAM20C+FKBP10+KDELR2+MATN3+MMP13+P3H1+SERPINH1+SLC10A7+SLC26A2+SLC35D1+CCN6+SERPINF1, data = data_setD)
# summary(LM_18)
pathvar_BMgene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_BMgene_age_6, "pathvar_BMgene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_age_6.csv --path Emily-folder/results_2/pathvar_BMgene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ BMgene, data = data_setE)
# summary(LM_15)
pathvar_BMgene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_BMgene_age_7, "pathvar_BMgene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_age_7.csv --path Emily-folder/results_2/pathvar_BMgene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+BMgene, data = data_setE)
# summary(LM_16)
pathvar_BMgene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_BMgene_age_8, "pathvar_BMgene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_age_8.csv --path Emily-folder/results_2/pathvar_BMgene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ FN1+ACAN+ADAMTS2+B3GALT6+B3GAT3+B4GALT7+CHST3+COL10A1+COL11A1+COL11A2+COL1A1+COL27A1+COL2A1+CRTAP+FAM20C+FKBP10+KDELR2+MATN3+MMP13+P3H1+SERPINH1+SLC10A7+SLC26A2+SLC35D1+CCN6+SERPINF1, data = data_setE)
# summary(LM_17)
pathvar_BMgene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_BMgene_age_9, "pathvar_BMgene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_age_9.csv --path Emily-folder/results_2/pathvar_BMgene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+FN1+ACAN+ADAMTS2+B3GALT6+B3GAT3+B4GALT7+CHST3+COL10A1+COL11A1+COL11A2+COL1A1+COL27A1+COL2A1+CRTAP+FAM20C+FKBP10+KDELR2+MATN3+MMP13+P3H1+SERPINH1+SLC10A7+SLC26A2+SLC35D1+CCN6+SERPINF1, data = data_setE)
# summary(LM_18)
pathvar_BMgene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_BMgene_age_10, "pathvar_BMgene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_BMgene_age_10.csv --path Emily-folder/results_2/pathvar_BMgene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(df_cancer_diagnosis, table(BMgene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(BMgene, Ethnicity_White)))
# p-value = 1.427e-10
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.6814063 0.8138085
# sample estimates:
#   odds ratio 
# 0.744308

with(df_cancer_diagnosis, table(BMgene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(BMgene, Ever_smoked_dummy)))
# p-value = 0.3859

with(df_cancer_diagnosis, table(BMgene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(BMgene, Sex_Male)))
# p-value = 0.43

with(df_cancer_diagnosis, table(BMgene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(BMgene, University_education)))
# p-value = 0.05333

with(df_cancer_diagnosis, table(BMgene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(BMgene, Annual_income)))
# p-value = 0.3685


### data_setD
with(data_setD, table(BMgene, Ethnicity_White))			
fisher.test(with(data_setD, table(BMgene, Ethnicity_White)))
# p-value = 2.224e-10
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.6828087 0.8163273
# sample estimates:
#   odds ratio 
# 0.7462133

with(data_setD, table(BMgene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(BMgene, Ever_smoked_dummy)))
# p-value = 0.3719

with(data_setD, table(BMgene, Sex_Male))			
fisher.test(with(data_setD, table(BMgene, Sex_Male)))
# p-value = 0.3746


### data_setE
with(data_setE, table(BMgene, Ethnicity_White))			
fisher.test(with(data_setE, table(BMgene, Ethnicity_White)))
# p-value = 1.483e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.7007405 0.8739161
# sample estimates:
#   odds ratio 
# 0.7818924 

with(data_setE, table(BMgene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(BMgene, Ever_smoked_dummy)))
# p-value = 0.4407

with(data_setE, table(BMgene, Sex_Male))			
fisher.test(with(data_setE, table(BMgene, Sex_Male)))
# p-value = 0.8875

with(data_setE, table(BMgene, University_education))			
fisher.test(with(data_setE, table(BMgene, University_education)))
# p-value = 0.0297
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.003911 1.080485
# sample estimates:
#   odds ratio 
# 1.041515

with(data_setE, table(BMgene, Annual_income))			
fisher.test(with(data_setE, table(BMgene, Annual_income)))
# p-value = 0.2281



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 

# make BMgene into factor
df_cancer_diagnosis$BMgene <- as.factor(df_cancer_diagnosis$BMgene)

# try a violin plot



# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ BMgene, data = df_cancer_diagnosis)
# 



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ BMgene, data = df_cancer_diagnosis)
# W = 897885488, p-value = 0.7147



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = BMgene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = BMgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ BMgene, data = df_cancer_diagnosis)
# W = 894556848, p-value = 0.2078



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = BMgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = BMgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ BMgene, data = df_cancer_diagnosis)
# W = 901700238, p-value = 0.01304



# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(BMgene), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$BMgene <- as.factor(data_setD$BMgene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ BMgene, data = data_setD)
# 



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ BMgene, data = data_setD)
# W = 882903574, p-value = 0.649




### data_setE

data_setE$BMgene <- as.factor(data_setE$BMgene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ BMgene, data = data_setE)
# 



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = BMgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ BMgene, data = data_setE)
# W = 569742750, p-value = 0.3863



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = BMgene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = BMgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ BMgene, data = data_setE)
# W = 570221686, p-value = 0.3025



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = BMgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = BMgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ BMgene, data = data_setE)
# W = 570665870, p-value = 0.2359


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(BMgene), scales = "free_y") +
  theme_light()
# 












# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_BMgenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_BMgenes.R')
