# 24/06/2024

# EntireCohort_pathvar_RNAPgenes.R

# this script is for the RNA processing associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)


# RNAPgene
# CRIPT
# WDR73
# TOE1
# THOC6
# EXOSC3
# PUF60
# PPIL1
# RNU4ATAC
# CDC40
# CLP1
# AARS1
# CTU2
# GON7
# LAGE3
# OSGEP
# POP1
# SEPSECS
# THG1L
# TP53RK
# TPRKB
# TRMT1
# TRMT10A
# TSEN15
# TSEN2
# TSEN54
# WDR4
# YARS1
# YRDC
# RRP7A














# RNAPgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')


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

















# RNAPgene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_RNAPgene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_RNAPgene_, "pathvar_RNAPgene_.csv", row.names=FALSE)
# system('dx upload pathvar_RNAPgene_.csv --path Emily-folder/results_2/pathvar_RNAPgene_.csv')


# RNAPgene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ RNAPgene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_RNAPgene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_incidence_1, "pathvar_RNAPgene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_incidence_1.csv --path Emily-folder/results_2/pathvar_RNAPgene_incidence_1.csv')


# GOIs
LR_2 <- glm(formula = cancer_ROOC_dummy ~ CRIPT+WDR73+TOE1+THOC6+EXOSC3+PUF60+PPIL1+CDC40+CLP1+AARS1+CTU2+GON7+LAGE3+OSGEP+POP1+SEPSECS+THG1L+TP53RK+TPRKB+TRMT1+TRMT10A+TSEN15+TSEN2+TSEN54+WDR4+YRDC+RRP7A, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_RNAPgene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_incidence_2, "pathvar_RNAPgene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_incidence_2.csv --path Emily-folder/results_2/pathvar_RNAPgene_incidence_2.csv')







### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "RNAPgene", 'CRIPT', 'WDR73', 'TOE1', 'THOC6', 
                                     'EXOSC3', 'PUF60', 'PPIL1', 'CDC40', 'CLP1', 
                                     'AARS1', 'CTU2', 'GON7', 'LAGE3', 'OSGEP', 'POP1', 
                                     'SEPSECS', 'THG1L', 'TP53RK', 'TPRKB', 'TRMT1', 'TRMT10A', 
                                     'TSEN15', 'TSEN2', 'TSEN54', 'WDR4', 'YRDC', 
                                     'RRP7A')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     34

with(data_setA, table(LAGE3))

# RNAPgene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ RNAPgene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_RNAPgene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_incidence_3, "pathvar_RNAPgene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_incidence_3.csv --path Emily-folder/results_2/pathvar_RNAPgene_incidence_3.csv')


# age + ethnicity + smoking + sex + RNAPgene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_RNAPgene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_incidence_4, "pathvar_RNAPgene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_incidence_4.csv --path Emily-folder/results_2/pathvar_RNAPgene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ CRIPT+WDR73+TOE1+THOC6+EXOSC3+PUF60+PPIL1+CDC40+CLP1+AARS1+CTU2+GON7+LAGE3+OSGEP+POP1+SEPSECS+THG1L+TP53RK+TPRKB+TRMT1+TRMT10A+TSEN15+TSEN2+TSEN54+WDR4+YRDC+RRP7A, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_RNAPgene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_incidence_5, "pathvar_RNAPgene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_incidence_5.csv --path Emily-folder/results_2/pathvar_RNAPgene_incidence_5.csv')


# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CRIPT+WDR73+TOE1+THOC6+EXOSC3+PUF60+PPIL1+CDC40+CLP1+AARS1+CTU2+GON7+LAGE3+OSGEP+POP1+SEPSECS+THG1L+TP53RK+TPRKB+TRMT1+TRMT10A+TSEN15+TSEN2+TSEN54+WDR4+YRDC+RRP7A, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_RNAPgene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_incidence_6, "pathvar_RNAPgene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_incidence_6.csv --path Emily-folder/results_2/pathvar_RNAPgene_incidence_6.csv')




### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "RNAPgene", 'CRIPT', 'WDR73', 'TOE1', 'THOC6', 
                                     'EXOSC3', 'PUF60', 'PPIL1', 'CDC40', 'CLP1', 
                                     'AARS1', 'CTU2', 'GON7', 'OSGEP', 'POP1', 
                                     'SEPSECS', 'THG1L', 'TP53RK', 'TPRKB', 'TRMT1', 'TRMT10A', 
                                     'TSEN15', 'TSEN2', 'TSEN54', 'WDR4', 'YRDC', 
                                     'RRP7A')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     39

with(data_setB, table(RRP7A))
# PUF60 = 83
# LAGE3 = 1 (removed)
# TPRKB = 71


# RNAPgene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ RNAPgene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_RNAPgene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_incidence_7, "pathvar_RNAPgene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_incidence_7.csv --path Emily-folder/results_2/pathvar_RNAPgene_incidence_7.csv')


# all variables + RNAPgene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+RNAPgene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_RNAPgene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_incidence_8, "pathvar_RNAPgene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_incidence_8.csv --path Emily-folder/results_2/pathvar_RNAPgene_incidence_8.csv')


# GOIs alone
LR_9 <- glm(formula = cancer_ROOC_dummy ~ CRIPT+WDR73+TOE1+THOC6+EXOSC3+PUF60+PPIL1+CDC40+CLP1+AARS1+CTU2+GON7+OSGEP+POP1+SEPSECS+THG1L+TP53RK+TPRKB+TRMT1+TRMT10A+TSEN15+TSEN2+TSEN54+WDR4+YRDC+RRP7A, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_RNAPgene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_incidence_9, "pathvar_RNAPgene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_incidence_9.csv --path Emily-folder/results_2/pathvar_RNAPgene_incidence_9.csv')


# all variables + GOIs
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+CRIPT+WDR73+TOE1+THOC6+EXOSC3+PUF60+PPIL1+CDC40+CLP1+AARS1+CTU2+GON7+OSGEP+POP1+SEPSECS+THG1L+TP53RK+TPRKB+TRMT1+TRMT10A+TSEN15+TSEN2+TSEN54+WDR4+YRDC+RRP7A, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_RNAPgene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_incidence_10, "pathvar_RNAPgene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_incidence_10.csv --path Emily-folder/results_2/pathvar_RNAPgene_incidence_10.csv')














### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(RNAPgene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(RNAPgene, cancer_ROOC_dummy)))
# p-value = 0.8149

with(EntireCohort_pathvar, table(RNAPgene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(RNAPgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.6793092 0.7504393
# sample estimates:
#   odds ratio 
# 0.7138361

with(EntireCohort_pathvar, table(RNAPgene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(RNAPgene, Ever_smoked_dummy)))
# p-value = 0.2621

with(EntireCohort_pathvar, table(RNAPgene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(RNAPgene, Sex_Male)))
# p-value = 0.04323
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.000754 1.052701
# sample estimates:
#   odds ratio 
# 1.026403

with(EntireCohort_pathvar, table(RNAPgene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(RNAPgene, University_education)))
# p-value = 0.6013

with(EntireCohort_pathvar, table(RNAPgene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(RNAPgene, Annual_income))	)
# p-value = 0.06781

with(EntireCohort_pathvar, table(RNAPgene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(RNAPgene, Smoker_previous)))
# p-value = 0.2266

with(EntireCohort_pathvar, table(RNAPgene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(RNAPgene, Smoker_current))	)
# p-value = 0.3306



### data_setA
with(data_setA, table(RNAPgene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(RNAPgene, cancer_ROOC_dummy)))
# p-value = 0.8615

with(data_setA, table(RNAPgene, Ethnicity_White))			
fisher.test(with(data_setA, table(RNAPgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.6769670 0.7480529
# sample estimates:
#   odds ratio 
# 0.7114639

with(data_setA, table(RNAPgene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(RNAPgene, Ever_smoked_dummy)))
# p-value = 0.2558

with(data_setA, table(RNAPgene, Sex_Male))			
fisher.test(with(data_setA, table(RNAPgene, Sex_Male)))
# p-value = 0.04109
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.000991 1.053172
# sample estimates:
#   odds ratio 
# 1.026748



### data_setB
with(data_setB, table(RNAPgene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(RNAPgene, cancer_ROOC_dummy)))
# p-value = 0.8724

with(data_setB, table(RNAPgene, Ethnicity_White))			
fisher.test(with(data_setB, table(RNAPgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.6601943 0.7438800
# sample estimates:
#   odds ratio 
# 0.7005404

with(data_setB, table(RNAPgene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(RNAPgene, Ever_smoked_dummy)))
# p-value = 0.2769

with(data_setB, table(RNAPgene, Sex_Male))			
fisher.test(with(data_setB, table(RNAPgene, Sex_Male)))
# p-value = 0.0161
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.006301 1.064442
# sample estimates:
#   odds ratio 
# 1.03498 

with(data_setB, table(RNAPgene, University_education))			
fisher.test(with(data_setB, table(RNAPgene, University_education)))
# p-value = 0.6024

with(data_setB, table(RNAPgene, Annual_income))			
fisher.test(with(data_setB, table(RNAPgene, Annual_income)))
# p-value = 0.08214






### can plot the continuous variables: 

# make RNAPgene into factor
df_for_plot <- EntireCohort_pathvar
df_for_plot$RNAPgene <- as.factor(df_for_plot$RNAPgene)

# try a violin plot

# Age_at_recruitment  
df_for_plot %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ RNAPgene, data = df_for_plot)
# W = 5718427246, p-value = 0.6262

# BMI  
df_for_plot %>% 
  ggplot(mapping = aes(x = RNAPgene, y = BMI)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = RNAPgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ RNAPgene, data = df_for_plot)
# W = 5627243151, p-value = 0.1029

# Standing_height 
df_for_plot %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ RNAPgene, data = df_for_plot)
# W = 5639213827, p-value = 0.1029

# make this a bar plot
# Weekly_exercise
df_for_plot %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(RNAPgene), scales = "free_y") +
  theme_light()
# 



### data_setA

# make RNAPgene into factor
data_setA$RNAPgene <- as.factor(data_setA$RNAPgene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ RNAPgene, data = data_setA)
# W = 5625626388, p-value = 0.6268



### data_setB

# make RNAPgene into factor
data_setB$RNAPgene <- as.factor(data_setB$RNAPgene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ RNAPgene, data = data_setB)
# W = 3753621151, p-value = 0.4829

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = RNAPgene, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = RNAPgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ RNAPgene, data = data_setB)
# W = 3729712346, p-value = 0.3957

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ RNAPgene, data = data_setB)
# W = 3710293816, p-value = 0.03489


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(RNAPgene), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# RNAPgene Type ------------------------------------------------------------

# make sure the execute the code in the 'RNAPgenes' section

## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_RNAPgene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_type_1, "pathvar_RNAPgene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_type_1.csv --path Emily-folder/results_2/pathvar_RNAPgene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_RNAPgene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_type_2, "pathvar_RNAPgene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_type_2.csv --path Emily-folder/results_2/pathvar_RNAPgene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_RNAPgene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_type_3, "pathvar_RNAPgene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_type_3.csv --path Emily-folder/results_2/pathvar_RNAPgene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_RNAPgene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_type_4, "pathvar_RNAPgene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_type_4.csv --path Emily-folder/results_2/pathvar_RNAPgene_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_RNAPgene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_type_5, "pathvar_RNAPgene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_type_5.csv --path Emily-folder/results_2/pathvar_RNAPgene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_RNAPgene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_type_6, "pathvar_RNAPgene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_type_6.csv --path Emily-folder/results_2/pathvar_RNAPgene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_RNAPgene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_type_7, "pathvar_RNAPgene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_type_7.csv --path Emily-folder/results_2/pathvar_RNAPgene_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_RNAPgene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_type_8, "pathvar_RNAPgene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_type_8.csv --path Emily-folder/results_2/pathvar_RNAPgene_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_RNAPgene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_type_9, "pathvar_RNAPgene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_type_9.csv --path Emily-folder/results_2/pathvar_RNAPgene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_RNAPgene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_type_10, "pathvar_RNAPgene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_type_10.csv --path Emily-folder/results_2/pathvar_RNAPgene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_RNAPgene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_type_11, "pathvar_RNAPgene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_type_11.csv --path Emily-folder/results_2/pathvar_RNAPgene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_RNAPgene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_type_12, "pathvar_RNAPgene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_type_12.csv --path Emily-folder/results_2/pathvar_RNAPgene_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_RNAPgene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_type_13, "pathvar_RNAPgene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_type_13.csv --path Emily-folder/results_2/pathvar_RNAPgene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_RNAPgene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_type_14, "pathvar_RNAPgene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_type_14.csv --path Emily-folder/results_2/pathvar_RNAPgene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_RNAPgene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_RNAPgene_type_15, "pathvar_RNAPgene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_type_15.csv --path Emily-folder/results_2/pathvar_RNAPgene_type_15.csv')




# TSEN15 type -------------------------------------------------------------


LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TSEN15, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_TSEN15_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TSEN15_type_1, "pathvar_TSEN15_type_1.csv", row.names=FALSE)
system('dx upload pathvar_TSEN15_type_1.csv --path Emily-folder/results_2/pathvar_TSEN15_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TSEN15, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_TSEN15_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TSEN15_type_2, "pathvar_TSEN15_type_2.csv", row.names=FALSE)
system('dx upload pathvar_TSEN15_type_2.csv --path Emily-folder/results_2/pathvar_TSEN15_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TSEN15, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_TSEN15_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TSEN15_type_3, "pathvar_TSEN15_type_3.csv", row.names=FALSE)
system('dx upload pathvar_TSEN15_type_3.csv --path Emily-folder/results_2/pathvar_TSEN15_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TSEN15, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_TSEN15_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TSEN15_type_4, "pathvar_TSEN15_type_4.csv", row.names=FALSE)
system('dx upload pathvar_TSEN15_type_4.csv --path Emily-folder/results_2/pathvar_TSEN15_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TSEN15, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_TSEN15_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TSEN15_type_5, "pathvar_TSEN15_type_5.csv", row.names=FALSE)
system('dx upload pathvar_TSEN15_type_5.csv --path Emily-folder/results_2/pathvar_TSEN15_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TSEN15, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_TSEN15_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TSEN15_type_6, "pathvar_TSEN15_type_6.csv", row.names=FALSE)
system('dx upload pathvar_TSEN15_type_6.csv --path Emily-folder/results_2/pathvar_TSEN15_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TSEN15, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_TSEN15_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TSEN15_type_7, "pathvar_TSEN15_type_7.csv", row.names=FALSE)
system('dx upload pathvar_TSEN15_type_7.csv --path Emily-folder/results_2/pathvar_TSEN15_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TSEN15, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_TSEN15_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TSEN15_type_8, "pathvar_TSEN15_type_8.csv", row.names=FALSE)
system('dx upload pathvar_TSEN15_type_8.csv --path Emily-folder/results_2/pathvar_TSEN15_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TSEN15, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_TSEN15_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TSEN15_type_9, "pathvar_TSEN15_type_9.csv", row.names=FALSE)
system('dx upload pathvar_TSEN15_type_9.csv --path Emily-folder/results_2/pathvar_TSEN15_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TSEN15, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_TSEN15_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TSEN15_type_10, "pathvar_TSEN15_type_10.csv", row.names=FALSE)
system('dx upload pathvar_TSEN15_type_10.csv --path Emily-folder/results_2/pathvar_TSEN15_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TSEN15, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_TSEN15_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TSEN15_type_11, "pathvar_TSEN15_type_11.csv", row.names=FALSE)
system('dx upload pathvar_TSEN15_type_11.csv --path Emily-folder/results_2/pathvar_TSEN15_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TSEN15, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_TSEN15_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TSEN15_type_12, "pathvar_TSEN15_type_12.csv", row.names=FALSE)
system('dx upload pathvar_TSEN15_type_12.csv --path Emily-folder/results_2/pathvar_TSEN15_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TSEN15, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_TSEN15_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TSEN15_type_13, "pathvar_TSEN15_type_13.csv", row.names=FALSE)
system('dx upload pathvar_TSEN15_type_13.csv --path Emily-folder/results_2/pathvar_TSEN15_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TSEN15, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_TSEN15_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TSEN15_type_14, "pathvar_TSEN15_type_14.csv", row.names=FALSE)
system('dx upload pathvar_TSEN15_type_14.csv --path Emily-folder/results_2/pathvar_TSEN15_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TSEN15, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_TSEN15_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TSEN15_type_15, "pathvar_TSEN15_type_15.csv", row.names=FALSE)
system('dx upload pathvar_TSEN15_type_15.csv --path Emily-folder/results_2/pathvar_TSEN15_type_15.csv')





# RNAPgene Age -------------------------------------------------------------

# make sure the execute the code in the 'RNAPgenes' section

### linear regression
# "RNAPgene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ RNAPgene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_RNAPgene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_RNAPgene_age_1, "pathvar_RNAPgene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_age_1.csv --path Emily-folder/results_2/pathvar_RNAPgene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ CRIPT+WDR73+TOE1+THOC6+EXOSC3+PUF60+PPIL1+CDC40+CLP1+AARS1+CTU2+GON7+LAGE3+OSGEP+POP1+SEPSECS+THG1L+TP53RK+TPRKB+TRMT1+TRMT10A+TSEN15+TSEN2+TSEN54+WDR4+YRDC+RRP7A, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_RNAPgene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_RNAPgene_age_2, "pathvar_RNAPgene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_age_2.csv --path Emily-folder/results_2/pathvar_RNAPgene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "RNAPgene", 'CRIPT', 'WDR73', 'TOE1', 'THOC6', 'EXOSC3', 
                                     'PUF60', 'PPIL1', 'CDC40', 'CLP1', 'AARS1', 
                                     'CTU2', 'GON7', 'OSGEP', 'POP1', 'SEPSECS', 
                                     'THG1L', 'TP53RK', 'TPRKB', 'TRMT1', 'TRMT10A', 'TSEN15', 'TSEN2', 
                                     'TSEN54', 'WDR4', 'YRDC', 'RRP7A')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     34

with(data_setD, table(RRP7A))
# CRIPT = 79
# PUF60 = 28
# PPIL1 = 34
# GON7 = 31
# LAGE3 = 1 (removed this)
# TP53RK = 91
# TPRKB = 20
# TSEN15 = 40
# WDR4 = 72
# YRDC = 92




data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "RNAPgene", 'CRIPT', 'WDR73', 'TOE1', 
                                     'THOC6', 'EXOSC3', 'PUF60', 'PPIL1', 'CDC40', 'CLP1', 
                                     'AARS1', 'CTU2', 'GON7', 'OSGEP', 'POP1', 'SEPSECS', 
                                     'THG1L', 'TP53RK', 'TPRKB', 'TRMT1', 'TRMT10A', 'TSEN15', 'TSEN2', 'TSEN54', 
                                     'WDR4', 'YRDC', 'RRP7A')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    38

with(data_setE, table(RRP7A))
# CRIPT = 64
# PUF60 = 22
# PPIL1 = 26
# CDC40 = 85
# GON7 = 24
# TP53RK = 69
# TPRKB = 20
# TSEN15 = 37
# WDR4 = 52
# YRDC = 73





### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ RNAPgene, data = data_setD)
# summary(LM_15)
pathvar_RNAPgene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_RNAPgene_age_3, "pathvar_RNAPgene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_age_3.csv --path Emily-folder/results_2/pathvar_RNAPgene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RNAPgene, data = data_setD)
# summary(LM_16)
pathvar_RNAPgene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_RNAPgene_age_4, "pathvar_RNAPgene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_age_4.csv --path Emily-folder/results_2/pathvar_RNAPgene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ CRIPT+WDR73+TOE1+THOC6+EXOSC3+PUF60+PPIL1+CDC40+CLP1+AARS1+CTU2+GON7+OSGEP+POP1+SEPSECS+THG1L+TP53RK+TPRKB+TRMT1+TRMT10A+TSEN15+TSEN2+TSEN54+WDR4+YRDC+RRP7A, data = data_setD)
# summary(LM_17)
pathvar_RNAPgene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_RNAPgene_age_5, "pathvar_RNAPgene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_age_5.csv --path Emily-folder/results_2/pathvar_RNAPgene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CRIPT+WDR73+TOE1+THOC6+EXOSC3+PUF60+PPIL1+CDC40+CLP1+AARS1+CTU2+GON7+OSGEP+POP1+SEPSECS+THG1L+TP53RK+TPRKB+TRMT1+TRMT10A+TSEN15+TSEN2+TSEN54+WDR4+YRDC+RRP7A, data = data_setD)
# summary(LM_18)
pathvar_RNAPgene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_RNAPgene_age_6, "pathvar_RNAPgene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_age_6.csv --path Emily-folder/results_2/pathvar_RNAPgene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ RNAPgene, data = data_setE)
# summary(LM_15)
pathvar_RNAPgene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_RNAPgene_age_7, "pathvar_RNAPgene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_age_7.csv --path Emily-folder/results_2/pathvar_RNAPgene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+RNAPgene, data = data_setE)
# summary(LM_16)
pathvar_RNAPgene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_RNAPgene_age_8, "pathvar_RNAPgene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_age_8.csv --path Emily-folder/results_2/pathvar_RNAPgene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ CRIPT+WDR73+TOE1+THOC6+EXOSC3+PUF60+PPIL1+CDC40+CLP1+AARS1+CTU2+GON7+OSGEP+POP1+SEPSECS+THG1L+TP53RK+TPRKB+TRMT1+TRMT10A+TSEN15+TSEN2+TSEN54+WDR4+YRDC+RRP7A, data = data_setE)
# summary(LM_17)
pathvar_RNAPgene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_RNAPgene_age_9, "pathvar_RNAPgene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_age_9.csv --path Emily-folder/results_2/pathvar_RNAPgene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+CRIPT+WDR73+TOE1+THOC6+EXOSC3+PUF60+PPIL1+CDC40+CLP1+AARS1+CTU2+GON7+OSGEP+POP1+SEPSECS+THG1L+TP53RK+TPRKB+TRMT1+TRMT10A+TSEN15+TSEN2+TSEN54+WDR4+YRDC+RRP7A, data = data_setE)
# summary(LM_18)
pathvar_RNAPgene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_RNAPgene_age_10, "pathvar_RNAPgene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_RNAPgene_age_10.csv --path Emily-folder/results_2/pathvar_RNAPgene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(df_cancer_diagnosis, table(RNAPgene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(RNAPgene, Ethnicity_White)))
# p-value = 1.608e-06
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.6103626 0.8095933
# sample estimates:
#   odds ratio 
# 0.7015499

with(df_cancer_diagnosis, table(RNAPgene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(RNAPgene, Ever_smoked_dummy)))
# p-value = 0.5452

with(df_cancer_diagnosis, table(RNAPgene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(RNAPgene, Sex_Male)))
# p-value = 0.3133

with(df_cancer_diagnosis, table(RNAPgene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(RNAPgene, University_education)))
# p-value = 0.4165

with(df_cancer_diagnosis, table(RNAPgene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(RNAPgene, Annual_income)))
# p-value = 0.9858


### data_setD
with(data_setD, table(RNAPgene, Ethnicity_White))			
fisher.test(with(data_setD, table(RNAPgene, Ethnicity_White)))
# p-value = 1.155e-06
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.6063530 0.8048711
# sample estimates:
#   odds ratio 
# 0.6972096

with(data_setD, table(RNAPgene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(RNAPgene, Ever_smoked_dummy)))
# p-value = 0.5726

with(data_setD, table(RNAPgene, Sex_Male))			
fisher.test(with(data_setD, table(RNAPgene, Sex_Male)))
# p-value = 0.3375


### data_setE
with(data_setE, table(RNAPgene, Ethnicity_White))			
fisher.test(with(data_setE, table(RNAPgene, Ethnicity_White)))
# p-value = 0.0001117
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5953574 0.8428789
# sample estimates:
#   odds ratio 
# 0.7061922

with(data_setE, table(RNAPgene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(RNAPgene, Ever_smoked_dummy)))
# p-value = 0.9509

with(data_setE, table(RNAPgene, Sex_Male))			
fisher.test(with(data_setE, table(RNAPgene, Sex_Male)))
# p-value = 0.1867

with(data_setE, table(RNAPgene, University_education))			
fisher.test(with(data_setE, table(RNAPgene, University_education)))
# p-value = 0.2198

with(data_setE, table(RNAPgene, Annual_income))			
fisher.test(with(data_setE, table(RNAPgene, Annual_income)))
# p-value = 0.9712



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 

# make RNAPgene into factor
df_cancer_diagnosis$RNAPgene <- as.factor(df_cancer_diagnosis$RNAPgene)

# try a violin plot



# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ RNAPgene, data = df_cancer_diagnosis)
# 



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ RNAPgene, data = df_cancer_diagnosis)
# W = 318825638, p-value = 0.4847



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = RNAPgene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = RNAPgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ RNAPgene, data = df_cancer_diagnosis)
# W = 316410188, p-value = 0.5757



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ RNAPgene, data = df_cancer_diagnosis)
# W = 314139572, p-value = 0.5189



# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(RNAPgene), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$RNAPgene <- as.factor(data_setD$RNAPgene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ RNAPgene, data = data_setD)
#



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ RNAPgene, data = data_setD)
# W = 313786886, p-value = 0.4805




### data_setE

data_setE$RNAPgene <- as.factor(data_setE$RNAPgene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ RNAPgene, data = data_setE)
# 



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ RNAPgene, data = data_setE)
# W = 203971274, p-value = 0.1127



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = RNAPgene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = RNAPgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ RNAPgene, data = data_setE)
# W = 202968596, p-value = 0.3161



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = RNAPgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ RNAPgene, data = data_setE)
# W = 199285864, p-value = 0.256


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(RNAPgene), scales = "free_y") +
  theme_light()
# 






# TSEN15 extras -----------------------------------------------------------



## fisher test for the linear regression data (TSEN15)

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 115361


with(df_cancer_diagnosis, table(TSEN15, Ethnicity_White))	
# Ethnicity_White
# TSEN15      0      1
# 0   3090 111775
# 1      0     40
fisher.test(with(df_cancer_diagnosis, table(TSEN15, Ethnicity_White)))
# p-value = 0.6284

with(df_cancer_diagnosis, table(TSEN15, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(TSEN15, Ever_smoked_dummy)))
# p = = 0.6265

with(df_cancer_diagnosis, table(TSEN15, Sex_Male))
# Sex_Male
# TSEN15     0     1
# 0 62219 53102
# 1    12    28
fisher.test(with(df_cancer_diagnosis, table(TSEN15, Sex_Male)))
# p-value = 0.00371
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.346180 5.903397
# sample estimates:
#   odds ratio 
# 2.733917 

with(df_cancer_diagnosis, table(TSEN15, University_education))			
fisher.test(with(df_cancer_diagnosis, table(TSEN15, University_education)))
# p-value = 0.3054

with(df_cancer_diagnosis, table(TSEN15, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(TSEN15, Annual_income)))
# p-value = 0.4167


### data_setD
# TSEN15 n=40
with(data_setD, table(TSEN15, Ethnicity_White))			
fisher.test(with(data_setD, table(TSEN15, Ethnicity_White)))
# p-value = 0.6283

with(data_setD, table(TSEN15, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(TSEN15, Ever_smoked_dummy)))
# p-value = 0.6264

with(data_setD, table(TSEN15, Sex_Male))	
# Sex_Male
# TSEN15     0     1
# 0 61758 52584
# 1    12    28
fisher.test(with(data_setD, table(TSEN15, Sex_Male)))
# p-value = 0.002432
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.349369 5.917162
# sample estimates:
#   odds ratio 
# 2.740392


### data_setE
with(data_setE, table(TSEN15, Ethnicity_White))			
fisher.test(with(data_setE, table(TSEN15, Ethnicity_White)))
# p-value = 1

with(data_setE, table(TSEN15, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(TSEN15, Ever_smoked_dummy)))
# p-value = 0.4951

with(data_setE, table(TSEN15, Sex_Male))	
# Sex_Male
# TSEN15     0     1
# 0 47614 44462
# 1    12    25
fisher.test(with(data_setE, table(TSEN15, Sex_Male)))
# p-value = 0.02102
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.079953 4.875124
# sample estimates:
#   odds ratio 
# 2.231017 

with(data_setE, table(TSEN15, University_education))			
fisher.test(with(data_setE, table(TSEN15, University_education)))
# p-value = 0.6015

with(data_setE, table(TSEN15, Annual_income))			
fisher.test(with(data_setE, table(TSEN15, Annual_income)))
# p-value = 0.4243



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 

# make TSEN15 into factor
df_cancer_diagnosis$TSEN15 <- as.factor(df_cancer_diagnosis$TSEN15)

# try a violin plot



# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TSEN15, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TSEN15, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ TSEN15, data = df_cancer_diagnosis)
# W = 2639538, p-value = 0.1137



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TSEN15, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TSEN15, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ TSEN15, data = df_cancer_diagnosis)
# W = 2299879, p-value = 0.9752



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TSEN15, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TSEN15, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ TSEN15, data = df_cancer_diagnosis)
# W = 2074778, p-value = 0.2874



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TSEN15, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TSEN15, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ TSEN15, data = df_cancer_diagnosis)
# W = 1636078, p-value = 0.00155
# plots saved in screenshots on laptop


# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(TSEN15), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$TSEN15 <- as.factor(data_setD$TSEN15)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = TSEN15, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = TSEN15, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ TSEN15, data = data_setD)
#W = 2615446, p-value = 0.1155



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = TSEN15, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = TSEN15, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ TSEN15, data = data_setD)
# W = 2277806, p-value = 0.9654




### data_setE

data_setE$TSEN15 <- as.factor(data_setE$TSEN15)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = TSEN15, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = TSEN15, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ TSEN15, data = data_setE)
# W = 1976238, p-value = 0.09158



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = TSEN15, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = TSEN15, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ TSEN15, data = data_setE)
# W = 1658534, p-value = 0.7812



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = TSEN15, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = TSEN15, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ TSEN15, data = data_setE)
# W = 1581552, p-value = 0.4511



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = TSEN15, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = TSEN15, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ TSEN15, data = data_setE)
# W = 1261656, p-value = 0.006278
# plots saved as screenshot on laptop

# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(TSEN15), scales = "free_y") +
  theme_light()
# 






## looking at proportions of variants in the TSEN15 group

# all vars
TSEN15_df <- EntireCohort_pathvar[EntireCohort_pathvar$TSEN15 == 1,]
ALL <- TSEN15_df$VAR
# 162

# all setA vars
data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     'TSEN15','VAR')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 324450      8

TSEN15_df_A <- data_setA[data_setA$TSEN15 == 1,]
setA <- TSEN15_df_A$VAR
# 162

# all setB vars
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     'TSEN15', 'VAR')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 264922     13

TSEN15_df_B <- data_setB[data_setB$TSEN15 == 1,]
setB <- TSEN15_df_B$VAR
# 130

# all setD vars
data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     'TSEN15', 'VAR')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 79270     8

TSEN15_df_D <- data_setD[data_setD$TSEN15 == 1,]
setD <- TSEN15_df_D$VAR
# 40

# all setE vars
data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", 'TSEN15','VAR')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 63717    13

TSEN15_df_E <- data_setE[data_setE$TSEN15 == 1,]
setE <- TSEN15_df_E$VAR
# 37

# save these as csv files to work with in R on my own computer
write.csv(ALL, "TSEN15_VAR_all.csv", row.names=FALSE)
system('dx upload TSEN15_VAR_all.csv --path Emily-folder/results_2/TSEN15_VAR_all.csv')

write.csv(setA, "TSEN15_VAR_setA.csv", row.names=FALSE)
system('dx upload TSEN15_VAR_setA.csv --path Emily-folder/results_2/TSEN15_VAR_setA.csv')

write.csv(setB, "TSEN15_VAR_setB.csv", row.names=FALSE)
system('dx upload TSEN15_VAR_setB.csv --path Emily-folder/results_2/TSEN15_VAR_setB.csv')

write.csv(setD, "TSEN15_VAR_setD.csv", row.names=FALSE)
system('dx upload TSEN15_VAR_setD.csv --path Emily-folder/results_2/TSEN15_VAR_setD.csv')

write.csv(setE, "TSEN15_VAR_setE.csv", row.names=FALSE)
system('dx upload TSEN15_VAR_setE.csv --path Emily-folder/results_2/TSEN15_VAR_setE.csv')






# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_RNAPgenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_RNAPgenes.R')
