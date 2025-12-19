# 24/06/2024

# EntireCohort_pathvar_LTUgenes.R

# this script is for the Lysosomal trafficking and ubiquitination associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)


# LTUgene
# GORAB
# IER3IP1
# COMP
# UFSP2
# SPOP
# OBSL1
# DYM
# SBF1
# RNF13
# EXOC6B
# RAB33B
# TRAPPC6B
# ARCN1
# COPB2
# CTSA
# GALNS
# GLB1
# GNPTAB
# GUSB
# ARSB
# CUL7
# CCDC8
# TBC1D23
# UBR1
# STAMBP
















# LTUgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')


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

















# LTUgene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_LTUgene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_LTUgene_, "pathvar_LTUgene_.csv", row.names=FALSE)
# system('dx upload pathvar_LTUgene_.csv --path Emily-folder/results_2/pathvar_LTUgene_.csv')


# LTUgene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ LTUgene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_LTUgene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_incidence_1, "pathvar_LTUgene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_incidence_1.csv --path Emily-folder/results_2/pathvar_LTUgene_incidence_1.csv')


# GOIs
LR_2 <- glm(formula = cancer_ROOC_dummy ~ GORAB+IER3IP1+COMP+UFSP2+OBSL1+DYM+SBF1+EXOC6B+RAB33B+TRAPPC6B+ARCN1+COPB2+CTSA+GALNS+GLB1+GNPTAB+GUSB+ARSB+CUL7+CCDC8+TBC1D23+UBR1+STAMBP, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_LTUgene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_incidence_2, "pathvar_LTUgene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_incidence_2.csv --path Emily-folder/results_2/pathvar_LTUgene_incidence_2.csv')







### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "LTUgene", 'GORAB', 'IER3IP1', 'COMP', 'UFSP2', 
                                     'OBSL1', 'DYM', 'SBF1', 'EXOC6B', 
                                     'RAB33B', 'TRAPPC6B', 'ARCN1', 'COPB2', 'CTSA', 'GALNS', 
                                     'GLB1', 'GNPTAB', 'GUSB', 'ARSB', 'CUL7', 'CCDC8', 
                                     'TBC1D23', 'UBR1', 'STAMBP')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
#465881     30

with(data_setA, table(STAMBP))




# LTUgene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ LTUgene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_LTUgene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_incidence_3, "pathvar_LTUgene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_incidence_3.csv --path Emily-folder/results_2/pathvar_LTUgene_incidence_3.csv')


# age + ethnicity + smoking + sex + LTUgene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_LTUgene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_incidence_4, "pathvar_LTUgene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_incidence_4.csv --path Emily-folder/results_2/pathvar_LTUgene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ GORAB+IER3IP1+COMP+UFSP2+OBSL1+DYM+SBF1+EXOC6B+RAB33B+TRAPPC6B+ARCN1+COPB2+CTSA+GALNS+GLB1+GNPTAB+GUSB+ARSB+CUL7+CCDC8+TBC1D23+UBR1+STAMBP, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_LTUgene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_incidence_5, "pathvar_LTUgene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_incidence_5.csv --path Emily-folder/results_2/pathvar_LTUgene_incidence_5.csv')


# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GORAB+IER3IP1+COMP+UFSP2+OBSL1+DYM+SBF1+EXOC6B+RAB33B+TRAPPC6B+ARCN1+COPB2+CTSA+GALNS+GLB1+GNPTAB+GUSB+ARSB+CUL7+CCDC8+TBC1D23+UBR1+STAMBP, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_LTUgene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_incidence_6, "pathvar_LTUgene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_incidence_6.csv --path Emily-folder/results_2/pathvar_LTUgene_incidence_6.csv')




### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "LTUgene", 'GORAB', 'IER3IP1', 'COMP', 'UFSP2', 
                                     'OBSL1', 'DYM', 'SBF1', 'EXOC6B', 
                                     'RAB33B', 'TRAPPC6B', 'ARCN1', 'COPB2', 'CTSA', 'GALNS', 
                                     'GLB1', 'GNPTAB', 'GUSB', 'ARSB', 'CUL7', 'CCDC8', 
                                     'TBC1D23', 'UBR1', 'STAMBP')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     35

with(data_setB, table(STAMBP))




# LTUgene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ LTUgene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_LTUgene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_incidence_7, "pathvar_LTUgene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_incidence_7.csv --path Emily-folder/results_2/pathvar_LTUgene_incidence_7.csv')


# all variables + LTUgene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+LTUgene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_LTUgene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_incidence_8, "pathvar_LTUgene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_incidence_8.csv --path Emily-folder/results_2/pathvar_LTUgene_incidence_8.csv')


# GOIs alone
LR_9 <- glm(formula = cancer_ROOC_dummy ~ GORAB+IER3IP1+COMP+UFSP2+OBSL1+DYM+SBF1+EXOC6B+RAB33B+TRAPPC6B+ARCN1+COPB2+CTSA+GALNS+GLB1+GNPTAB+GUSB+ARSB+CUL7+CCDC8+TBC1D23+UBR1+STAMBP, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_LTUgene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_incidence_9, "pathvar_LTUgene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_incidence_9.csv --path Emily-folder/results_2/pathvar_LTUgene_incidence_9.csv')


# all variables + GOIs
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+GORAB+IER3IP1+COMP+UFSP2+OBSL1+DYM+SBF1+EXOC6B+RAB33B+TRAPPC6B+ARCN1+COPB2+CTSA+GALNS+GLB1+GNPTAB+GUSB+ARSB+CUL7+CCDC8+TBC1D23+UBR1+STAMBP, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_LTUgene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_incidence_10, "pathvar_LTUgene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_incidence_10.csv --path Emily-folder/results_2/pathvar_LTUgene_incidence_10.csv')










### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(LTUgene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(LTUgene, cancer_ROOC_dummy)))
# p-value = 0.7501

with(EntireCohort_pathvar, table(LTUgene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(LTUgene, Ethnicity_White)))
# p-value = 4.693e-14
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.8093177 0.8823656
# sample estimates:
#   odds ratio 
# 0.8449221 

with(EntireCohort_pathvar, table(LTUgene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(LTUgene, Ever_smoked_dummy)))
# p-value = 0.1557

with(EntireCohort_pathvar, table(LTUgene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(LTUgene, Sex_Male)))
# p-value = 0.0831

with(EntireCohort_pathvar, table(LTUgene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(LTUgene, University_education)))
# p-value = 0.6636

with(EntireCohort_pathvar, table(LTUgene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(LTUgene, Annual_income))	)
# p-value = 0.1406

with(EntireCohort_pathvar, table(LTUgene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(LTUgene, Smoker_previous)))
# p-value = 0.132

with(EntireCohort_pathvar, table(LTUgene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(LTUgene, Smoker_current))	)
# p-value = 0.9797



### data_setA
with(data_setA, table(LTUgene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(LTUgene, cancer_ROOC_dummy)))
# p-value = 0.7584

with(data_setA, table(LTUgene, Ethnicity_White))			
fisher.test(with(data_setA, table(LTUgene, Ethnicity_White)))
# p-value = 6.973e-14
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.8095479 0.8829345
# sample estimates:
#   odds ratio 
# 0.8453152

with(data_setA, table(LTUgene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(LTUgene, Ever_smoked_dummy)))
# p-value = 0.1596

with(data_setA, table(LTUgene, Sex_Male))			
fisher.test(with(data_setA, table(LTUgene, Sex_Male)))
# p-value = 0.07384



### data_setB
with(data_setB, table(LTUgene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(LTUgene, cancer_ROOC_dummy)))
# p-value = 0.6865

with(data_setB, table(LTUgene, Ethnicity_White))			
fisher.test(with(data_setB, table(LTUgene, Ethnicity_White)))
# p-value = 1.039e-09
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.8052955 0.8939100
# sample estimates:
#   odds ratio 
# 0.8482412

with(data_setB, table(LTUgene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(LTUgene, Ever_smoked_dummy)))
# p-value = 0.3631

with(data_setB, table(LTUgene, Sex_Male))			
fisher.test(with(data_setB, table(LTUgene, Sex_Male)))
# p-value = 0.07791

with(data_setB, table(LTUgene, University_education))			
fisher.test(with(data_setB, table(LTUgene, University_education)))
# p-value = 0.8377

with(data_setB, table(LTUgene, Annual_income))			
fisher.test(with(data_setB, table(LTUgene, Annual_income)))
# p-value = 0.2137






### can plot the continuous variables: 

# make LTUgene into factor
df_for_plots <- EntireCohort_pathvar
df_for_plots$LTUgene <- as.factor(df_for_plots$LTUgene)

# try a violin plot

# Age_at_recruitment  
df_for_plots %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ LTUgene, data = df_for_plots)
# W = 8705300080, p-value = 0.3839

# BMI  
df_for_plots %>% 
  ggplot(mapping = aes(x = LTUgene, y = BMI)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = LTUgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ LTUgene, data = df_for_plots)
# W = 8631423526, p-value = 0.4764

# Standing_height 
df_for_plots %>% 
  ggplot(mapping = aes(x = LTUgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = LTUgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ LTUgene, data = df_for_plots)
# W = 8624206826, p-value = 0.7625

# make this a bar plot
# Weekly_exercise
df_for_plots %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(LTUgene), scales = "free_y") +
  theme_light()
# 



### data_setA

# make LTUgene into factor
data_setA$LTUgene <- as.factor(data_setA$LTUgene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ LTUgene, data = data_setA)
# W = 8562106120, p-value = 0.3463



### data_setB

# make LTUgene into factor
data_setB$LTUgene <- as.factor(data_setB$LTUgene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ LTUgene, data = data_setB)
# W = 5722503313, p-value = 0.6327

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = LTUgene, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = LTUgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ LTUgene, data = data_setB)
# W = 5725520337, p-value = 0.5248

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = LTUgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = LTUgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ LTUgene, data = data_setB)
# W = 5694851430, p-value = 0.33


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(LTUgene), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# LTUgene Type ------------------------------------------------------------

# make sure the execute the code in the 'LTUgenes' section

## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_LTUgene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_type_1, "pathvar_LTUgene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_type_1.csv --path Emily-folder/results_2/pathvar_LTUgene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_LTUgene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_type_2, "pathvar_LTUgene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_type_2.csv --path Emily-folder/results_2/pathvar_LTUgene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_LTUgene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_type_3, "pathvar_LTUgene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_type_3.csv --path Emily-folder/results_2/pathvar_LTUgene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_LTUgene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_type_4, "pathvar_LTUgene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_type_4.csv --path Emily-folder/results_2/pathvar_LTUgene_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_LTUgene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_type_5, "pathvar_LTUgene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_type_5.csv --path Emily-folder/results_2/pathvar_LTUgene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_LTUgene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_type_6, "pathvar_LTUgene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_type_6.csv --path Emily-folder/results_2/pathvar_LTUgene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_LTUgene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_type_7, "pathvar_LTUgene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_type_7.csv --path Emily-folder/results_2/pathvar_LTUgene_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_LTUgene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_type_8, "pathvar_LTUgene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_type_8.csv --path Emily-folder/results_2/pathvar_LTUgene_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_LTUgene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_type_9, "pathvar_LTUgene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_type_9.csv --path Emily-folder/results_2/pathvar_LTUgene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_LTUgene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_type_10, "pathvar_LTUgene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_type_10.csv --path Emily-folder/results_2/pathvar_LTUgene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_LTUgene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_type_11, "pathvar_LTUgene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_type_11.csv --path Emily-folder/results_2/pathvar_LTUgene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_LTUgene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_type_12, "pathvar_LTUgene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_type_12.csv --path Emily-folder/results_2/pathvar_LTUgene_type_12.csv')



## carry out the endocrine type analysis for GOIs
LR_endocrine_GOIs <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GORAB+IER3IP1+COMP+UFSP2+OBSL1+DYM+SBF1+EXOC6B+RAB33B+TRAPPC6B+ARCN1+COPB2+CTSA+GALNS+GLB1+GNPTAB+GUSB+ARSB+CUL7+CCDC8+TBC1D23+UBR1+STAMBP, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_GOIs)
# exp(cbind(OR = coef(LR_endocrine_GOIs), confint(LR_endocrine_GOIs)))
pathvar_LTUgene_type_GOIs <- tidy(LR_endocrine_GOIs, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_type_GOIs, "pathvar_LTUgene_type_GOIs", row.names=FALSE)
system('dx upload pathvar_LTUgene_type_GOIs --path Emily-folder/results_2/pathvar_LTUgene_type_GOIs.csv')
##



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_LTUgene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_type_13, "pathvar_LTUgene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_type_13.csv --path Emily-folder/results_2/pathvar_LTUgene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_LTUgene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_type_14, "pathvar_LTUgene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_type_14.csv --path Emily-folder/results_2/pathvar_LTUgene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_LTUgene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_LTUgene_type_15, "pathvar_LTUgene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_type_15.csv --path Emily-folder/results_2/pathvar_LTUgene_type_15.csv')





# CUL7 Type -------------------------------------------------------------
# repeat the type analysis for just CUL7


LR_oral_CUL7 <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CUL7, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral_CUL7)
# exp(cbind(OR = coef(LR_oral_CUL7), confint(LR_oral_CUL7)))
pathvar_CUL7_type_1 <- tidy(LR_oral_CUL7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CUL7_type_1, "pathvar_CUL7_type_1.csv", row.names=FALSE)
system('dx upload pathvar_CUL7_type_1.csv --path Emily-folder/results_2/pathvar_CUL7_type_1.csv')



LR_digestive_CUL7 <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CUL7, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive_CUL7)
# exp(cbind(OR = coef(LR_digestive_CUL7), confint(LR_digestive_CUL7)))
pathvar_CUL7_type_2 <- tidy(LR_digestive_CUL7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CUL7_type_2, "pathvar_CUL7_type_2.csv", row.names=FALSE)
system('dx upload pathvar_CUL7_type_2.csv --path Emily-folder/results_2/pathvar_CUL7_type_2.csv')



LR_respiratory_CUL7 <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CUL7, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory_CUL7)
# exp(cbind(OR = coef(LR_respiratory_CUL7), confint(LR_respiratory_CUL7)))
pathvar_CUL7_type_3 <- tidy(LR_respiratory_CUL7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CUL7_type_3, "pathvar_CUL7_type_3.csv", row.names=FALSE)
system('dx upload pathvar_CUL7_type_3.csv --path Emily-folder/results_2/pathvar_CUL7_type_3.csv')



LR_bone_CUL7 <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CUL7, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone_CUL7)
# exp(cbind(OR = coef(LR_bone_CUL7), confint(LR_bone_CUL7)))
pathvar_CUL7_type_4 <- tidy(LR_bone_CUL7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CUL7_type_4, "pathvar_CUL7_type_4.csv", row.names=FALSE)
system('dx upload pathvar_CUL7_type_4.csv --path Emily-folder/results_2/pathvar_CUL7_type_4.csv')



LR_skin_CUL7 <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CUL7, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin_CUL7)
# exp(cbind(OR = coef(LR_skin_CUL7), confint(LR_skin_CUL7)))
pathvar_CUL7_type_5 <- tidy(LR_skin_CUL7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CUL7_type_5, "pathvar_CUL7_type_5.csv", row.names=FALSE)
system('dx upload pathvar_CUL7_type_5.csv --path Emily-folder/results_2/pathvar_CUL7_type_5.csv')



LR_mesothelium_CUL7 <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CUL7, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium_CUL7)
# exp(cbind(OR = coef(LR_mesothelium_CUL7), confint(LR_mesothelium_CUL7)))
pathvar_CUL7_type_6 <- tidy(LR_mesothelium_CUL7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CUL7_type_6, "pathvar_CUL7_type_6.csv", row.names=FALSE)
system('dx upload pathvar_CUL7_type_6.csv --path Emily-folder/results_2/pathvar_CUL7_type_6.csv')



LR_breast_CUL7 <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CUL7, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast_CUL7)
# exp(cbind(OR = coef(LR_breast_CUL7), confint(LR_breast_CUL7)))
pathvar_CUL7_type_7 <- tidy(LR_breast_CUL7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CUL7_type_7, "pathvar_CUL7_type_7.csv", row.names=FALSE)
system('dx upload pathvar_CUL7_type_7.csv --path Emily-folder/results_2/pathvar_CUL7_type_7.csv')



LR_femalegenital_CUL7 <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CUL7, 
                             data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital_CUL7)
# exp(cbind(OR = coef(LR_femalegenital_CUL7), confint(LR_femalegenital_CUL7)))
pathvar_CUL7_type_8 <- tidy(LR_femalegenital_CUL7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CUL7_type_8, "pathvar_CUL7_type_8.csv", row.names=FALSE)
system('dx upload pathvar_CUL7_type_8.csv --path Emily-folder/results_2/pathvar_CUL7_type_8.csv')



LR_malegenital_CUL7 <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CUL7, 
                           data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital_CUL7)
# exp(cbind(OR = coef(LR_malegenital_CUL7), confint(LR_malegenital_CUL7)))
pathvar_CUL7_type_9 <- tidy(LR_malegenital_CUL7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CUL7_type_9, "pathvar_CUL7_type_9.csv", row.names=FALSE)
system('dx upload pathvar_CUL7_type_9.csv --path Emily-folder/results_2/pathvar_CUL7_type_9.csv')



LR_urinary_CUL7 <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CUL7, 
                       data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary_CUL7)
# exp(cbind(OR = coef(LR_urinary_CUL7), confint(LR_urinary_CUL7)))
pathvar_CUL7_type_10 <- tidy(LR_urinary_CUL7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CUL7_type_10, "pathvar_CUL7_type_10.csv", row.names=FALSE)
system('dx upload pathvar_CUL7_type_10.csv --path Emily-folder/results_2/pathvar_CUL7_type_10.csv')



LR_cns_CUL7 <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CUL7, 
                   data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns_CUL7)
# exp(cbind(OR = coef(LR_cns_CUL7), confint(LR_cns_CUL7)))
pathvar_CUL7_type_11 <- tidy(LR_cns_CUL7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CUL7_type_11, "pathvar_CUL7_type_11.csv", row.names=FALSE)
system('dx upload pathvar_CUL7_type_11.csv --path Emily-folder/results_2/pathvar_CUL7_type_11.csv')



LR_endocrine_CUL7 <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CUL7, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine_CUL7)
# exp(cbind(OR = coef(LR_endocrine_CUL7), confint(LR_endocrine_CUL7)))
pathvar_CUL7_type_12 <- tidy(LR_endocrine_CUL7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CUL7_type_12, "pathvar_CUL7_type_12.csv", row.names=FALSE)
system('dx upload pathvar_CUL7_type_12.csv --path Emily-folder/results_2/pathvar_CUL7_type_12.csv')



LR_lymphatic_CUL7 <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CUL7, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic_CUL7)
# exp(cbind(OR = coef(LR_lymphatic_CUL7), confint(LR_lymphatic_CUL7)))
pathvar_CUL7_type_13 <- tidy(LR_lymphatic_CUL7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CUL7_type_13, "pathvar_CUL7_type_13.csv", row.names=FALSE)
system('dx upload pathvar_CUL7_type_13.csv --path Emily-folder/results_2/pathvar_CUL7_type_13.csv')



LR_secondary_CUL7 <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CUL7, 
                         data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary_CUL7)
# exp(cbind(OR = coef(LR_secondary_CUL7), confint(LR_secondary_CUL7)))
pathvar_CUL7_type_14 <- tidy(LR_secondary_CUL7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CUL7_type_14, "pathvar_CUL7_type_14.csv", row.names=FALSE)
system('dx upload pathvar_CUL7_type_14.csv --path Emily-folder/results_2/pathvar_CUL7_type_14.csv')



LR_other_CUL7 <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CUL7, 
                     data = EntireCohort_pathvar, family = binomial)
# summary(LR_other_CUL7)
# exp(cbind(OR = coef(LR_other_CUL7), confint(LR_other_CUL7)))
pathvar_CUL7_type_15 <- tidy(LR_other_CUL7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_CUL7_type_15, "pathvar_CUL7_type_15.csv", row.names=FALSE)
system('dx upload pathvar_CUL7_type_15.csv --path Emily-folder/results_2/pathvar_CUL7_type_15.csv')














# LTUgene Age -------------------------------------------------------------

# make sure the execute the code in the 'LTUgenes' section

### linear regression
# "LTUgene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ LTUgene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_LTUgene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_LTUgene_age_1, "pathvar_LTUgene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_age_1.csv --path Emily-folder/results_2/pathvar_LTUgene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ GORAB+IER3IP1+COMP+UFSP2+OBSL1+DYM+SBF1+EXOC6B+RAB33B+TRAPPC6B+ARCN1+COPB2+CTSA+GALNS+GLB1+GNPTAB+GUSB+ARSB+CUL7+CCDC8+TBC1D23+UBR1+STAMBP, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_LTUgene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_LTUgene_age_2, "pathvar_LTUgene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_age_2.csv --path Emily-folder/results_2/pathvar_LTUgene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "LTUgene", 'GORAB', 'IER3IP1', 'COMP', 'UFSP2', 
                                     'OBSL1', 'DYM', 'SBF1', 'EXOC6B', 'RAB33B', 
                                     'TRAPPC6B', 'ARCN1', 'COPB2', 'CTSA', 'GALNS', 'GLB1', 
                                     'GNPTAB', 'GUSB', 'ARSB', 'CUL7', 'CCDC8', 'TBC1D23', 'UBR1', 
                                     'STAMBP')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     30

with(data_setD, table(STAMBP))
# IER3IP1 = 44
# ARCN1 = 98




data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "LTUgene", 'GORAB', 'IER3IP1', 'COMP', 
                                     'UFSP2', 'OBSL1', 'DYM', 'SBF1', 'EXOC6B', 
                                     'RAB33B', 'TRAPPC6B', 'ARCN1', 'COPB2', 'CTSA', 'GALNS', 'GLB1', 
                                     'GNPTAB', 'GUSB', 'ARSB', 'CUL7', 'CCDC8', 'TBC1D23', 'UBR1', 'STAMBP')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    35

with(data_setE, table(STAMBP))
# GORAB = 94
# IER3IP1 = 35
# ARCN1 = 82
# GUSB = 98
# CCDC8 = 79




### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ LTUgene, data = data_setD)
# summary(LM_15)
pathvar_LTUgene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_LTUgene_age_3, "pathvar_LTUgene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_age_3.csv --path Emily-folder/results_2/pathvar_LTUgene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+LTUgene, data = data_setD)
# summary(LM_16)
pathvar_LTUgene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_LTUgene_age_4, "pathvar_LTUgene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_age_4.csv --path Emily-folder/results_2/pathvar_LTUgene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ GORAB+IER3IP1+COMP+UFSP2+OBSL1+DYM+SBF1+EXOC6B+RAB33B+TRAPPC6B+ARCN1+COPB2+CTSA+GALNS+GLB1+GNPTAB+GUSB+ARSB+CUL7+CCDC8+TBC1D23+UBR1+STAMBP, data = data_setD)
# summary(LM_17)
pathvar_LTUgene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_LTUgene_age_5, "pathvar_LTUgene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_age_5.csv --path Emily-folder/results_2/pathvar_LTUgene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GORAB+IER3IP1+COMP+UFSP2+OBSL1+DYM+SBF1+EXOC6B+RAB33B+TRAPPC6B+ARCN1+COPB2+CTSA+GALNS+GLB1+GNPTAB+GUSB+ARSB+CUL7+CCDC8+TBC1D23+UBR1+STAMBP, data = data_setD)
# summary(LM_18)
pathvar_LTUgene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_LTUgene_age_6, "pathvar_LTUgene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_age_6.csv --path Emily-folder/results_2/pathvar_LTUgene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ LTUgene, data = data_setE)
# summary(LM_15)
pathvar_LTUgene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_LTUgene_age_7, "pathvar_LTUgene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_age_7.csv --path Emily-folder/results_2/pathvar_LTUgene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+LTUgene, data = data_setE)
# summary(LM_16)
pathvar_LTUgene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_LTUgene_age_8, "pathvar_LTUgene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_age_8.csv --path Emily-folder/results_2/pathvar_LTUgene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ GORAB+IER3IP1+COMP+UFSP2+OBSL1+DYM+SBF1+EXOC6B+RAB33B+TRAPPC6B+ARCN1+COPB2+CTSA+GALNS+GLB1+GNPTAB+GUSB+ARSB+CUL7+CCDC8+TBC1D23+UBR1+STAMBP, data = data_setE)
# summary(LM_17)
pathvar_LTUgene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_LTUgene_age_9, "pathvar_LTUgene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_age_9.csv --path Emily-folder/results_2/pathvar_LTUgene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+GORAB+IER3IP1+COMP+UFSP2+OBSL1+DYM+SBF1+EXOC6B+RAB33B+TRAPPC6B+ARCN1+COPB2+CTSA+GALNS+GLB1+GNPTAB+GUSB+ARSB+CUL7+CCDC8+TBC1D23+UBR1+STAMBP, data = data_setE)
# summary(LM_18)
pathvar_LTUgene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_LTUgene_age_10, "pathvar_LTUgene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_LTUgene_age_10.csv --path Emily-folder/results_2/pathvar_LTUgene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(df_cancer_diagnosis, table(LTUgene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(LTUgene, Ethnicity_White)))
# p-value = 0.002066
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.7305312 0.9332364
# sample estimates:
#   odds ratio 
# 0.8245127

with(df_cancer_diagnosis, table(LTUgene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(LTUgene, Ever_smoked_dummy)))
# p-value = 0.1612

with(df_cancer_diagnosis, table(LTUgene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(LTUgene, Sex_Male)))
# p-value = 0.1868

with(df_cancer_diagnosis, table(LTUgene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(LTUgene, University_education)))
# p-value = 0.6903

with(df_cancer_diagnosis, table(LTUgene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(LTUgene, Annual_income)))
# p-value = 0.3967


### data_setD
with(data_setD, table(LTUgene, Ethnicity_White))			
fisher.test(with(data_setD, table(LTUgene, Ethnicity_White)))
# p-value = 0.003403
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.7357753 0.9417753
# sample estimates:
#   odds ratio 
# 0.8312341

with(data_setD, table(LTUgene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(LTUgene, Ever_smoked_dummy)))
# p-value = 0.1508

with(data_setD, table(LTUgene, Sex_Male))			
fisher.test(with(data_setD, table(LTUgene, Sex_Male)))
# p-value = 0.2191


### data_setE
with(data_setE, table(LTUgene, Ethnicity_White))			
fisher.test(with(data_setE, table(LTUgene, Ethnicity_White)))
# p-value = 0.01892
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.7191508 0.9719289
# sample estimates:
#   odds ratio 
# 0.8342167

with(data_setE, table(LTUgene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(LTUgene, Ever_smoked_dummy)))
# p-value = 0.4046

with(data_setE, table(LTUgene, Sex_Male))			
fisher.test(with(data_setE, table(LTUgene, Sex_Male)))
# p-value = 0.2602

with(data_setE, table(LTUgene, University_education))			
fisher.test(with(data_setE, table(LTUgene, University_education)))
# p-value = 0.2861

with(data_setE, table(LTUgene, Annual_income))			
fisher.test(with(data_setE, table(LTUgene, Annual_income)))
# p-value = 0.2433



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 

# make LTUgene into factor
df_cancer_diagnosis$LTUgene <- as.factor(df_cancer_diagnosis$LTUgene)

# try a violin plot



# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ LTUgene, data = df_cancer_diagnosis)
# 



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ LTUgene, data = df_cancer_diagnosis)
# W = 488551416, p-value = 0.2066



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = LTUgene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = LTUgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ LTUgene, data = df_cancer_diagnosis)
# W = 481592274, p-value = 0.7569



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = LTUgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = LTUgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ LTUgene, data = df_cancer_diagnosis)
# W = 484241932, p-value = 0.4638



# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(LTUgene), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$LTUgene <- as.factor(data_setD$LTUgene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ LTUgene, data = data_setD)
# 



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ LTUgene, data = data_setD)
# W = 479853180, p-value = 0.2661




### data_setE

data_setE$LTUgene <- as.factor(data_setE$LTUgene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ LTUgene, data = data_setE)
# 



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = LTUgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ LTUgene, data = data_setE)
# W = 311809958, p-value = 0.1694



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = LTUgene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = LTUgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ LTUgene, data = data_setE)
# W = 309972250, p-value = 0.6089



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = LTUgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = LTUgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ LTUgene, data = data_setE)
# W = 309493560, p-value = 0.7738


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(LTUgene), scales = "free_y") +
  theme_light()
# 










# CUL7 extras -----------------------------------------------------------



# look at the demographics of the CUL7 participants

### Fisher test on these datasets

### EntireCohort_pathvar
with(EntireCohort_pathvar, table(CUL7, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(CUL7, cancer_ROOC_dummy)))
#  p-value = 0.4622

with(EntireCohort_pathvar, table(CUL7, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(CUL7, Ethnicity_White)))
#  p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5434708 0.6806967
# sample estimates:
#   odds ratio 
# 0.6074215 

with(EntireCohort_pathvar, table(CUL7, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(CUL7, Ever_smoked_dummy)))
#  p-value = 0.2783

with(EntireCohort_pathvar, table(CUL7, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(CUL7, Sex_Male)))
#  p-value = 0.1729

with(EntireCohort_pathvar, table(CUL7, University_education))			
fisher.test(with(EntireCohort_pathvar, table(CUL7, University_education)))
#  p-value = 0.4928

with(EntireCohort_pathvar, table(CUL7, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(CUL7, Annual_income))	)
# p-value = 0.5885

with(EntireCohort_pathvar, table(CUL7, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(CUL7, Smoker_previous)))
#  p-value = 0.158

with(EntireCohort_pathvar, table(CUL7, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(CUL7, Smoker_current))	)
#  p-value = 0.525



### data_setD
with(data_setD, table(CUL7, Ethnicity_White))			
fisher.test(with(data_setD, table(CUL7, Ethnicity_White)))
#  p-value = 0.0001022
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.4002619 0.7345558
# sample estimates:
#   odds ratio 
# 0.5364788 

with(data_setD, table(CUL7, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(CUL7, Ever_smoked_dummy)))
#  p-value = 0.6671

with(data_setD, table(CUL7, Sex_Male))			
fisher.test(with(data_setD, table(CUL7, Sex_Male)))
#  p-value = 0.7008



### data_setE 

with(data_setE, table(CUL7, Ethnicity_White))			
fisher.test(with(data_setE, table(CUL7, Ethnicity_White)))
#  p-value = 0.0004198
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.3604818 0.7440173
# sample estimates:
#   odds ratio 
# 0.509954 

with(data_setE, table(CUL7, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(CUL7, Ever_smoked_dummy)))
#  p-value = 1

with(data_setE, table(CUL7, Sex_Male))			
fisher.test(with(data_setE, table(CUL7, Sex_Male)))
#  p-value = 0.3735

with(data_setE, table(CUL7, University_education))			
fisher.test(with(data_setE, table(CUL7, University_education)))
#  p-value = 0.1044

with(data_setE, table(CUL7, Annual_income))			
fisher.test(with(data_setE, table(CUL7, Annual_income)))
# p-value = 0.6336




### can plot the continuous variables: 

# make CUL7 into factor
df_to_plot <- EntireCohort_pathvar
df_to_plot$CUL7 <- as.factor(df_to_plot$CUL7)

# try a violin plot

# Age_at_recruitment  
df_to_plot %>% 
  ggplot(mapping = aes(x = CUL7, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = CUL7, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ CUL7, data = df_to_plot)
# W = 965256588, p-value = 0.8771

# BMI  
df_to_plot %>% 
  ggplot(mapping = aes(x = CUL7, y = BMI)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = CUL7, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ CUL7, data = df_to_plot)
#  W = 952012802, p-value = 0.5812

# Standing_height 
df_to_plot %>% 
  ggplot(mapping = aes(x = CUL7, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = CUL7, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ CUL7, data = df_to_plot)
#  W = 950944995, p-value = 0.3662

# make this a bar plot
# Weekly_exercise
df_to_plot %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(CUL7), scales = "free_y") +
  theme_light()
# 



### data_setD

# make CUL7 into factor
data_setD$CUL7 <- as.factor(data_setD$CUL7)


# jpeg(file="CUL7_Age_v_A.jpeg", width = 300, height = 300)
# 
# dev.off()
# system('dx upload CUL7_Age_v_A.jpeg --path Emily-folder/pathvar/plots/CUL7_Age_v_A.jpeg')



# Age_at_recruitment 

data_setD %>% 
  ggplot(mapping = aes(x = CUL7, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setD %>% 
  ggplot(mapping = aes(x = CUL7, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ CUL7, data = data_setD)
#  W = 53943782, p-value = 0.9691



### data_setE

# make CUL7 into factor
data_setE$CUL7 <- as.factor(data_setE$CUL7)

# Age_at_recruitment

data_setE %>% 
  ggplot(mapping = aes(x = CUL7, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setE %>% 
  ggplot(mapping = aes(x = CUL7, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ CUL7, data = data_setE)
#  W = 34560340, p-value = 0.5987



# BMI  
data_setE %>% 
  ggplot(mapping = aes(x = CUL7, y = BMI)) +
  geom_violin() +
  theme_light()
data_setE %>% 
  ggplot(mapping = aes(x = CUL7, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ CUL7, data = data_setE)
#   = 34428494, p-value = 0.4781



# Standing_height 

data_setE %>% 
  ggplot(mapping = aes(x = CUL7, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setE %>% 
  ggplot(mapping = aes(x = CUL7, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ CUL7, data = data_setE)
#  W = 34265071, p-value = 0.3485


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(CUL7), scales = "free_y") +
  theme_light()
# similar














# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_LTUgenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_LTUgenes.R')
