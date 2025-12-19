# 09/07/2024

# EntireCohort_pathvar_Bonegenes.R

# this script is for the Bone associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)


# Bonegene
# DYNC2H1
# SGMS2
# TMEM38B
# TRPV4
# ARSL
# HSPG2
# RSPRY1
# PIEZO2
# CANT1
# NOTCH2
# ACP5
# DDR2
# INPPL1
# LIFR
# NPPC
# NPR2
# PRKG2
# PPP3CA
# TNFRSF11B











# Bonegenes ---------------------------------------------------------------

# load in the dataframe
# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')


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
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE1 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2', 'TNFRSF11B')] <- 1
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE2 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2', 'TNFRSF11B')] <- 1
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE3 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2', 'TNFRSF11B')] <- 1
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE4 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2', 'TNFRSF11B')] <- 1
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE5 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2', 'TNFRSF11B')] <- 1
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE6 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2', 'TNFRSF11B')] <- 1
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE7 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2', 'TNFRSF11B')] <- 1
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE8 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2', 'TNFRSF11B')] <- 1
EntireCohort_pathvar$Bonegene[EntireCohort_pathvar$GENE9 %in% c('DYNC2H1', 'SGMS2', 'TMEM38B', 'TRPV4', 'ARSL', 'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2', 'TNFRSF11B')] <- 1


EntireCohort_pathvar$Bonegene[is.na(EntireCohort_pathvar$Bonegene)] <- 0
EntireCohort_pathvar$Bonegene <- as.numeric(EntireCohort_pathvar$Bonegene)
unique(EntireCohort_pathvar$Bonegene)
# 0 1

table(EntireCohort_pathvar$Bonegene)
#      0      1 
# 443607  26192 




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


### 'PPP3CA'
# no variants present





# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)

















# Bonegene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_Bonegene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_Bonegene_, "pathvar_Bonegene_.csv", row.names=FALSE)
# system('dx upload pathvar_Bonegene_.csv --path Emily-folder/results_2/pathvar_Bonegene_.csv')


# Bonegene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ Bonegene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_Bonegene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_incidence_1, "pathvar_Bonegene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_incidence_1.csv --path Emily-folder/results_2/pathvar_Bonegene_incidence_1.csv')


# GOIs
# excluding SGMS2 (only 2 participants)
LR_2 <- glm(formula = cancer_ROOC_dummy ~ DYNC2H1+TMEM38B+TRPV4+ARSL+HSPG2+RSPRY1+CANT1+ACP5+DDR2+INPPL1+LIFR+NPPC+NPR2+PRKG2+TNFRSF11B, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_Bonegene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_incidence_2, "pathvar_Bonegene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_incidence_2.csv --path Emily-folder/results_2/pathvar_Bonegene_incidence_2.csv')







### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "Bonegene", 'DYNC2H1', 'TMEM38B', 'TRPV4', 
                                     'ARSL', 'HSPG2', 'RSPRY1', 'CANT1',
                                     'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 
                                     'NPR2', 'PRKG2', 'TNFRSF11B')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     21

with(data_setA, table(Bonegene))
# Bonegene
# 0      1 
# 439943  25938
# TRPV4 = 4
# ARSL = 13
# NPPC = 45



# Bonegene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ Bonegene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_Bonegene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_incidence_3, "pathvar_Bonegene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_incidence_3.csv --path Emily-folder/results_2/pathvar_Bonegene_incidence_3.csv')


# age + ethnicity + smoking + sex + Bonegene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_Bonegene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_incidence_4, "pathvar_Bonegene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_incidence_4.csv --path Emily-folder/results_2/pathvar_Bonegene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ DYNC2H1+TMEM38B+TRPV4+ARSL+HSPG2+RSPRY1+CANT1+ACP5+DDR2+INPPL1+LIFR+NPPC+NPR2+PRKG2+TNFRSF11B, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_Bonegene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_incidence_5, "pathvar_Bonegene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_incidence_5.csv --path Emily-folder/results_2/pathvar_Bonegene_incidence_5.csv')


# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DYNC2H1+TMEM38B+TRPV4+ARSL+HSPG2+RSPRY1+CANT1+ACP5+DDR2+INPPL1+LIFR+NPPC+NPR2+PRKG2+TNFRSF11B, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_Bonegene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_incidence_6, "pathvar_Bonegene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_incidence_6.csv --path Emily-folder/results_2/pathvar_Bonegene_incidence_6.csv')




### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "Bonegene", 'DYNC2H1', 'TMEM38B', 'TRPV4', 'ARSL', 
                                     'HSPG2', 'RSPRY1', 'CANT1', 'ACP5', 
                                     'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2', 'TNFRSF11B')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     26

with(data_setB, table(Bonegene))
# Bonegene
# 0      1 
# 359894  21265
# TRPV4 = 3
# ARSL = 10
# NPPC = 38



# Bonegene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ Bonegene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_Bonegene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_incidence_7, "pathvar_Bonegene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_incidence_7.csv --path Emily-folder/results_2/pathvar_Bonegene_incidence_7.csv')


# all variables + Bonegene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+Bonegene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_Bonegene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_incidence_8, "pathvar_Bonegene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_incidence_8.csv --path Emily-folder/results_2/pathvar_Bonegene_incidence_8.csv')


# GOIs alone
LR_9 <- glm(formula = cancer_ROOC_dummy ~ DYNC2H1+TMEM38B+TRPV4+ARSL+HSPG2+RSPRY1+CANT1+ACP5+DDR2+INPPL1+LIFR+NPPC+NPR2+PRKG2+TNFRSF11B, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_Bonegene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_incidence_9, "pathvar_Bonegene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_incidence_9.csv --path Emily-folder/results_2/pathvar_Bonegene_incidence_9.csv')


# all variables + GOIs
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+DYNC2H1+TMEM38B+TRPV4+ARSL+HSPG2+RSPRY1+CANT1+ACP5+DDR2+INPPL1+LIFR+NPPC+NPR2+PRKG2+TNFRSF11B, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_Bonegene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_incidence_10, "pathvar_Bonegene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_incidence_10.csv --path Emily-folder/results_2/pathvar_Bonegene_incidence_10.csv')










### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(Bonegene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(Bonegene, cancer_ROOC_dummy)))
# p-value = 0.1505

with(EntireCohort_pathvar, table(Bonegene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(Bonegene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.6167039 0.6783719
# sample estimates:
#   odds ratio 
# 0.6466718

with(EntireCohort_pathvar, table(Bonegene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(Bonegene, Ever_smoked_dummy)))
# p-value = 0.4125

with(EntireCohort_pathvar, table(Bonegene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(Bonegene, Sex_Male)))
# p-value = 0.5444

with(EntireCohort_pathvar, table(Bonegene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(Bonegene, University_education)))
# p-value = 0.5259

with(EntireCohort_pathvar, table(Bonegene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(Bonegene, Annual_income))	)
# p-value = 0.6371

with(EntireCohort_pathvar, table(Bonegene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(Bonegene, Smoker_previous)))
# p-value = 0.2633

with(EntireCohort_pathvar, table(Bonegene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(Bonegene, Smoker_current))	)
# p-value = 0.1125



### data_setA
with(data_setA, table(Bonegene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(Bonegene, cancer_ROOC_dummy)))
# p-value = 0.1267

with(data_setA, table(Bonegene, Ethnicity_White))			
fisher.test(with(data_setA, table(Bonegene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.6165617 0.6784784
# sample estimates:
#   odds ratio 
# 0.6466468 

with(data_setA, table(Bonegene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(Bonegene, Ever_smoked_dummy)))
# p-value = 0.5021

with(data_setA, table(Bonegene, Sex_Male))			
fisher.test(with(data_setA, table(Bonegene, Sex_Male)))
# p-value = 0.4926



### data_setB
with(data_setB, table(Bonegene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(Bonegene, cancer_ROOC_dummy)))
# p-value = 0.2961

with(data_setB, table(Bonegene, Ethnicity_White))			
fisher.test(with(data_setB, table(Bonegene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.6095345 0.6835308
# sample estimates:
#   odds ratio 
# 0.6452959 

with(data_setB, table(Bonegene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(Bonegene, Ever_smoked_dummy)))
# p-value = 0.5784

with(data_setB, table(Bonegene, Sex_Male))			
fisher.test(with(data_setB, table(Bonegene, Sex_Male)))
# p-value = 0.5067

with(data_setB, table(Bonegene, University_education))			
fisher.test(with(data_setB, table(Bonegene, University_education)))
# p-value = 0.8656

with(data_setB, table(Bonegene, Annual_income))			
fisher.test(with(data_setB, table(Bonegene, Annual_income)))
# p-value = 0.6893






### can plot the continuous variables: 

# make Bonegene into factor
df_for_plots <- EntireCohort_pathvar
df_for_plots$Bonegene <- as.factor(df_for_plots$Bonegene)

# try a violin plot

# Age_at_recruitment  
df_for_plots %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ Bonegene, data = df_for_plots)
# W = 5.89e+09, p-value = 0.0001586

# BMI  
df_for_plots %>% 
  ggplot(mapping = aes(x = Bonegene, y = BMI)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = Bonegene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ Bonegene, data = df_for_plots)
# W = 5714525715, p-value = 0.03286

# Standing_height 
df_for_plots %>% 
  ggplot(mapping = aes(x = Bonegene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = Bonegene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ Bonegene, data = df_for_plots)
# W = 5899561895, p-value = 2.539e-09

# make this a bar plot
# Weekly_exercise
df_for_plots %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(Bonegene), scales = "free_y") +
  theme_light()
# 



### data_setA

# make Bonegene into factor
data_setA$Bonegene <- as.factor(data_setA$Bonegene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ Bonegene, data = data_setA)
# W = 5785618086, p-value = 0.0001428



### data_setB

# make Bonegene into factor
data_setB$Bonegene <- as.factor(data_setB$Bonegene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ Bonegene, data = data_setB)
# W = 3.872e+09, p-value = 0.003518

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = Bonegene, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = Bonegene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ Bonegene, data = data_setB)
# W = 3802209250, p-value = 0.1181

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = Bonegene, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = Bonegene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ Bonegene, data = data_setB)
# W = 3914301322, p-value = 1.811e-08


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(Bonegene), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# Bonegene Type ------------------------------------------------------------

# make sure the execute the code in the 'Bonegenes' section


## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_Bonegene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_type_1, "pathvar_Bonegene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_type_1.csv --path Emily-folder/results_2/pathvar_Bonegene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_Bonegene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_type_2, "pathvar_Bonegene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_type_2.csv --path Emily-folder/results_2/pathvar_Bonegene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_Bonegene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_type_3, "pathvar_Bonegene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_type_3.csv --path Emily-folder/results_2/pathvar_Bonegene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_Bonegene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_type_4, "pathvar_Bonegene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_type_4.csv --path Emily-folder/results_2/pathvar_Bonegene_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_Bonegene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_type_5, "pathvar_Bonegene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_type_5.csv --path Emily-folder/results_2/pathvar_Bonegene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_Bonegene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_type_6, "pathvar_Bonegene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_type_6.csv --path Emily-folder/results_2/pathvar_Bonegene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_Bonegene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_type_7, "pathvar_Bonegene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_type_7.csv --path Emily-folder/results_2/pathvar_Bonegene_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_Bonegene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_type_8, "pathvar_Bonegene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_type_8.csv --path Emily-folder/results_2/pathvar_Bonegene_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_Bonegene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_type_9, "pathvar_Bonegene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_type_9.csv --path Emily-folder/results_2/pathvar_Bonegene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_Bonegene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_type_10, "pathvar_Bonegene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_type_10.csv --path Emily-folder/results_2/pathvar_Bonegene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_Bonegene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_type_11, "pathvar_Bonegene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_type_11.csv --path Emily-folder/results_2/pathvar_Bonegene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_Bonegene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_type_12, "pathvar_Bonegene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_type_12.csv --path Emily-folder/results_2/pathvar_Bonegene_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_Bonegene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_type_13, "pathvar_Bonegene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_type_13.csv --path Emily-folder/results_2/pathvar_Bonegene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_Bonegene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_type_14, "pathvar_Bonegene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_type_14.csv --path Emily-folder/results_2/pathvar_Bonegene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_Bonegene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Bonegene_type_15, "pathvar_Bonegene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_type_15.csv --path Emily-folder/results_2/pathvar_Bonegene_type_15.csv')






# Bonegene Age -------------------------------------------------------------

# make sure the execute the code in the 'Bonegenes' section

### linear regression
# "Bonegene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ Bonegene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_Bonegene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_Bonegene_age_1, "pathvar_Bonegene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_age_1.csv --path Emily-folder/results_2/pathvar_Bonegene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ DYNC2H1+TMEM38B+TRPV4+ARSL+HSPG2+RSPRY1+CANT1+ACP5+DDR2+INPPL1+LIFR+NPPC+NPR2+PRKG2+TNFRSF11B, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_Bonegene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_Bonegene_age_2, "pathvar_Bonegene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_age_2.csv --path Emily-folder/results_2/pathvar_Bonegene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "Bonegene", 'DYNC2H1', 'TMEM38B',
                                     'HSPG2', 'RSPRY1', 'CANT1', 'ACP5',
                                     'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2', 'TNFRSF11B')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     19

with(data_setD, table(Bonegene))
# Bonegene
# 0      1 
# 108145   6237
# TMEM38B = 72
# TRPV4 = 1 (removed)
# ARSL = 1 (removed)
# RSPRY1 = 99
# ACP5 = 88
# NPPC = 8
# PRKG2 = 76




data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "Bonegene", 'DYNC2H1', 'TMEM38B', 
                                     'HSPG2', 'RSPRY1', 'CANT1', 
                                     'ACP5', 'DDR2', 'INPPL1', 'LIFR', 'NPPC', 'NPR2', 'PRKG2', 'TNFRSF11B')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    24

with(data_setE, table(Bonegene))
# Bonegene
# 0     1 
# 87067  5046
# TMEM38B = 56
# RSPRY1 = 72
# ACP5 = 73
# NPPC = 7
# PRKG2 = 60




### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ Bonegene, data = data_setD)
# summary(LM_15)
pathvar_Bonegene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_Bonegene_age_3, "pathvar_Bonegene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_age_3.csv --path Emily-folder/results_2/pathvar_Bonegene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Bonegene, data = data_setD)
# summary(LM_16)
pathvar_Bonegene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_Bonegene_age_4, "pathvar_Bonegene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_age_4.csv --path Emily-folder/results_2/pathvar_Bonegene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ DYNC2H1+TMEM38B+HSPG2+RSPRY1+CANT1+ACP5+DDR2+INPPL1+LIFR+NPPC+NPR2+PRKG2+TNFRSF11B, data = data_setD)
# summary(LM_17)
pathvar_Bonegene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_Bonegene_age_5, "pathvar_Bonegene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_age_5.csv --path Emily-folder/results_2/pathvar_Bonegene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DYNC2H1+TMEM38B+HSPG2+RSPRY1+CANT1+ACP5+DDR2+INPPL1+LIFR+NPPC+NPR2+PRKG2+TNFRSF11B, data = data_setD)
# summary(LM_18)
pathvar_Bonegene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_Bonegene_age_6, "pathvar_Bonegene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_age_6.csv --path Emily-folder/results_2/pathvar_Bonegene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ Bonegene, data = data_setE)
# summary(LM_15)
pathvar_Bonegene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_Bonegene_age_7, "pathvar_Bonegene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_age_7.csv --path Emily-folder/results_2/pathvar_Bonegene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+Bonegene, data = data_setE)
# summary(LM_16)
pathvar_Bonegene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_Bonegene_age_8, "pathvar_Bonegene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_age_8.csv --path Emily-folder/results_2/pathvar_Bonegene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ DYNC2H1+TMEM38B+HSPG2+RSPRY1+CANT1+ACP5+DDR2+INPPL1+LIFR+NPPC+NPR2+PRKG2+TNFRSF11B, data = data_setE)
# summary(LM_17)
pathvar_Bonegene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_Bonegene_age_9, "pathvar_Bonegene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_age_9.csv --path Emily-folder/results_2/pathvar_Bonegene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+DYNC2H1+TMEM38B+HSPG2+RSPRY1+CANT1+ACP5+DDR2+INPPL1+LIFR+NPPC+NPR2+PRKG2+TNFRSF11B, data = data_setE)
# summary(LM_18)
pathvar_Bonegene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_Bonegene_age_10, "pathvar_Bonegene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_Bonegene_age_10.csv --path Emily-folder/results_2/pathvar_Bonegene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(df_cancer_diagnosis, table(Bonegene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(Bonegene, Ethnicity_White)))
# p-value = 9.225e-13
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5339513 0.6937821
# sample estimates:
#   odds ratio 
# 0.607625 

with(df_cancer_diagnosis, table(Bonegene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(Bonegene, Ever_smoked_dummy)))
# p-value = 0.4503

with(df_cancer_diagnosis, table(Bonegene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(Bonegene, Sex_Male)))
# p-value = 0.7453

with(df_cancer_diagnosis, table(Bonegene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(Bonegene, University_education)))
# p-value = 0.5226

with(df_cancer_diagnosis, table(Bonegene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(Bonegene, Annual_income)))
# p-value = 0.3577


### data_setD
with(data_setD, table(Bonegene, Ethnicity_White))			
fisher.test(with(data_setD, table(Bonegene, Ethnicity_White)))
# p-value = 1.691e-12
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5346067 0.6956743
# sample estimates:
#   odds ratio 
# 0.6088049 

with(data_setD, table(Bonegene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(Bonegene, Ever_smoked_dummy)))
# p-value = 0.4909

with(data_setD, table(Bonegene, Sex_Male))			
fisher.test(with(data_setD, table(Bonegene, Sex_Male)))
# p-value = 0.8141


### data_setE
with(data_setE, table(Bonegene, Ethnicity_White))			
fisher.test(with(data_setE, table(Bonegene, Ethnicity_White)))
# p-value = 5.62e-10
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5075524 0.6967158
# sample estimates:
#   odds ratio 
# 0.5931896

with(data_setE, table(Bonegene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(Bonegene, Ever_smoked_dummy)))
# p-value = 0.7634

with(data_setE, table(Bonegene, Sex_Male))			
fisher.test(with(data_setE, table(Bonegene, Sex_Male)))
# p-value = 0.8734

with(data_setE, table(Bonegene, University_education))			
fisher.test(with(data_setE, table(Bonegene, University_education)))
# p-value = 0.3567

with(data_setE, table(Bonegene, Annual_income))			
fisher.test(with(data_setE, table(Bonegene, Annual_income)))
# p-value = 0.3883



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 

# make Bonegene into factor
df_cancer_diagnosis$Bonegene <- as.factor(df_cancer_diagnosis$Bonegene)

# try a violin plot



# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ Bonegene, data = df_cancer_diagnosis)
# W = 348547754, p-value = 0.06188



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ Bonegene, data = df_cancer_diagnosis)
# WW = 345127051, p-value = 0.5911



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Bonegene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Bonegene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ Bonegene, data = df_cancer_diagnosis)
# W = 337670329, p-value = 0.1563



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Bonegene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Bonegene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ Bonegene, data = df_cancer_diagnosis)
# W = 349190432, p-value = 0.006912



# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(Bonegene), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$Bonegene <- as.factor(data_setD$Bonegene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ Bonegene, data = data_setD)
# W = 342565478, p-value = 0.03606



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ Bonegene, data = data_setD)
# W = 339034221, p-value = 0.4812




### data_setE

data_setE$Bonegene <- as.factor(data_setE$Bonegene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ Bonegene, data = data_setE)
# W = 222754367, p-value = 0.09305



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = Bonegene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ Bonegene, data = data_setE)
# W = 221248270, p-value = 0.3896



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = Bonegene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = Bonegene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ Bonegene, data = data_setE)
#  = 216813923, p-value = 0.1199



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = Bonegene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = Bonegene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ Bonegene, data = data_setE)
# W = 224613871, p-value = 0.007075


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(Bonegene), scales = "free_y") +
  theme_light()
# 












# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_Bonegenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_Bonegenes.R')
