# 24/06/24

# EntireCohort_pathvar_TGFbBMPgenes.R

# this script is for the TGFb and BMP associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)


# TGFbBMPgene
# SLC39A13
# ADAMTSL2
# SOX9
# DDRGK1
# BGN
# LEMD3
# BMP2
# BMPR1B
# GDF5
# BMP1
# FBN1
# LTBP3









# TGFbBMPgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')


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

















# TGFbBMPgene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_TGFbBMPgene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_TGFbBMPgene_, "pathvar_TGFbBMPgene_.csv", row.names=FALSE)
# system('dx upload pathvar_TGFbBMPgene_.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_.csv')


# TGFbBMPgene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ TGFbBMPgene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_TGFbBMPgene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_incidence_1, "pathvar_TGFbBMPgene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_incidence_1.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_incidence_1.csv')


# GOIs
LR_2 <- glm(formula = cancer_ROOC_dummy ~ SLC39A13+ADAMTSL2+SOX9+DDRGK1+LEMD3+BMP2+BMPR1B+GDF5+BMP1+LTBP3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_TGFbBMPgene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_incidence_2, "pathvar_TGFbBMPgene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_incidence_2.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_incidence_2.csv')







### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "TGFbBMPgene", 'SLC39A13', 'ADAMTSL2', 'SOX9', 'DDRGK1', 
                                     'LEMD3', 'BMP2', 'BMPR1B', 'GDF5', 'BMP1', 'LTBP3')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     17

with(data_setA, table(LTBP3))
# all good


# TGFbBMPgene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ TGFbBMPgene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_TGFbBMPgene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_incidence_3, "pathvar_TGFbBMPgene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_incidence_3.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_incidence_3.csv')


# age + ethnicity + smoking + sex + TGFbBMPgene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_TGFbBMPgene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_incidence_4, "pathvar_TGFbBMPgene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_incidence_4.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ SLC39A13+ADAMTSL2+SOX9+DDRGK1+LEMD3+BMP2+BMPR1B+GDF5+BMP1+LTBP3, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_TGFbBMPgene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_incidence_5, "pathvar_TGFbBMPgene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_incidence_5.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_incidence_5.csv')


# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SLC39A13+ADAMTSL2+SOX9+DDRGK1+LEMD3+BMP2+BMPR1B+GDF5+BMP1+LTBP3, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_TGFbBMPgene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_incidence_6, "pathvar_TGFbBMPgene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_incidence_6.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_incidence_6.csv')




### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "TGFbBMPgene", 'SLC39A13', 'ADAMTSL2', 'SOX9', 'DDRGK1', 
                                     'LEMD3', 'BMP2', 'BMPR1B', 'GDF5', 'BMP1', 'LTBP3')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     22

with(data_setB, table(LTBP3))


# TGFbBMPgene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ TGFbBMPgene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_TGFbBMPgene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_incidence_7, "pathvar_TGFbBMPgene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_incidence_7.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_incidence_7.csv')


# all variables + TGFbBMPgene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+TGFbBMPgene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_TGFbBMPgene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_incidence_8, "pathvar_TGFbBMPgene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_incidence_8.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_incidence_8.csv')


# GOIs alone
LR_9 <- glm(formula = cancer_ROOC_dummy ~ SLC39A13+ADAMTSL2+SOX9+DDRGK1+LEMD3+BMP2+BMPR1B+GDF5+BMP1+LTBP3, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_TGFbBMPgene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_incidence_9, "pathvar_TGFbBMPgene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_incidence_9.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_incidence_9.csv')


# all variables + GOIs
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+SLC39A13+ADAMTSL2+SOX9+DDRGK1+LEMD3+BMP2+BMPR1B+GDF5+BMP1+LTBP3, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_TGFbBMPgene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_incidence_10, "pathvar_TGFbBMPgene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_incidence_10.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_incidence_10.csv')










### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(TGFbBMPgene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(TGFbBMPgene, cancer_ROOC_dummy)))
# p-value = 0.2423

with(EntireCohort_pathvar, table(TGFbBMPgene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(TGFbBMPgene, Ethnicity_White)))
# p-value = 0.1623

with(EntireCohort_pathvar, table(TGFbBMPgene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(TGFbBMPgene, Ever_smoked_dummy)))
# p-value = 0.7328

with(EntireCohort_pathvar, table(TGFbBMPgene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(TGFbBMPgene, Sex_Male)))
# p-value = 0.2338

with(EntireCohort_pathvar, table(TGFbBMPgene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(TGFbBMPgene, University_education)))
# p-value = 0.01012
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.011387 1.088652
# sample estimates:
#   odds ratio 
# 1.049361

with(EntireCohort_pathvar, table(TGFbBMPgene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(TGFbBMPgene, Annual_income))	)
# p-value = 0.6551

with(EntireCohort_pathvar, table(TGFbBMPgene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(TGFbBMPgene, Smoker_previous)))
# p-value = 0.5979

with(EntireCohort_pathvar, table(TGFbBMPgene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(TGFbBMPgene, Smoker_current))	)
# p-value = 0.7635



### data_setA
with(data_setA, table(TGFbBMPgene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(TGFbBMPgene, cancer_ROOC_dummy)))
# p-value = 0.2445

with(data_setA, table(TGFbBMPgene, Ethnicity_White))			
fisher.test(with(data_setA, table(TGFbBMPgene, Ethnicity_White)))
# p-value = 0.1332

with(data_setA, table(TGFbBMPgene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(TGFbBMPgene, Ever_smoked_dummy)))
# p-value = 0.7256

with(data_setA, table(TGFbBMPgene, Sex_Male))			
fisher.test(with(data_setA, table(TGFbBMPgene, Sex_Male)))
# p-value = 0.2607



### data_setB
with(data_setB, table(TGFbBMPgene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(TGFbBMPgene, cancer_ROOC_dummy)))
# p-value = 0.2236

with(data_setB, table(TGFbBMPgene, Ethnicity_White))			
fisher.test(with(data_setB, table(TGFbBMPgene, Ethnicity_White)))
# p-value = 0.5838

with(data_setB, table(TGFbBMPgene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(TGFbBMPgene, Ever_smoked_dummy)))
# p-value = 0.8258

with(data_setB, table(TGFbBMPgene, Sex_Male))			
fisher.test(with(data_setB, table(TGFbBMPgene, Sex_Male)))
# p-value = 0.2246

with(data_setB, table(TGFbBMPgene, University_education))			
fisher.test(with(data_setB, table(TGFbBMPgene, University_education)))
# p-value = 0.0124
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.010950 1.094991
# sample estimates:
#   odds ratio 
# 1.052179

with(data_setB, table(TGFbBMPgene, Annual_income))			
fisher.test(with(data_setB, table(TGFbBMPgene, Annual_income)))
# p-value = 0.8679






### can plot the continuous variables: 

# make TGFbBMPgene into factor
df_for_plot <- EntireCohort_pathvar
df_for_plot$TGFbBMPgene <- as.factor(df_for_plot$TGFbBMPgene)

# try a violin plot

# Age_at_recruitment  
df_for_plot %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ TGFbBMPgene, data = df_for_plot)
# W = 3065224244, p-value = 0.1236

# BMI  
df_for_plot %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = BMI)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ TGFbBMPgene, data = df_for_plot)
# W = 3008419671, p-value = 0.6079

# Standing_height 
df_for_plot %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ TGFbBMPgene, data = df_for_plot)
# W = 3066619265, p-value = 0.00376

# make this a bar plot
# Weekly_exercise
df_for_plot %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(TGFbBMPgene), scales = "free_y") +
  theme_light()
# 



### data_setA

# make TGFbBMPgene into factor
data_setA$TGFbBMPgene <- as.factor(data_setA$TGFbBMPgene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ TGFbBMPgene, data = data_setA)
# W = 3015866370, p-value = 0.1151



### data_setB

# make TGFbBMPgene into factor
data_setB$TGFbBMPgene <- as.factor(data_setB$TGFbBMPgene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ TGFbBMPgene, data = data_setB)
# W = 2012809339, p-value = 0.04281

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ TGFbBMPgene, data = data_setB)
# W = 1.984e+09, p-value = 0.5926

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ TGFbBMPgene, data = data_setB)
# W = 2020468201, p-value = 0.006798


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(TGFbBMPgene), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# TGFbBMPgene Type ------------------------------------------------------------

# make sure the execute the code in the 'TGFbBMPgenes' section

## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_TGFbBMPgene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_type_1, "pathvar_TGFbBMPgene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_type_1.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_TGFbBMPgene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_type_2, "pathvar_TGFbBMPgene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_type_2.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_TGFbBMPgene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_type_3, "pathvar_TGFbBMPgene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_type_3.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_TGFbBMPgene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_type_4, "pathvar_TGFbBMPgene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_type_4.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_TGFbBMPgene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_type_5, "pathvar_TGFbBMPgene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_type_5.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_TGFbBMPgene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_type_6, "pathvar_TGFbBMPgene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_type_6.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_TGFbBMPgene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_type_7, "pathvar_TGFbBMPgene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_type_7.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_TGFbBMPgene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_type_8, "pathvar_TGFbBMPgene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_type_8.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_TGFbBMPgene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_type_9, "pathvar_TGFbBMPgene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_type_9.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_TGFbBMPgene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_type_10, "pathvar_TGFbBMPgene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_type_10.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_TGFbBMPgene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_type_11, "pathvar_TGFbBMPgene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_type_11.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_TGFbBMPgene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_type_12, "pathvar_TGFbBMPgene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_type_12.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_TGFbBMPgene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_type_13, "pathvar_TGFbBMPgene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_type_13.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_TGFbBMPgene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_type_14, "pathvar_TGFbBMPgene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_type_14.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_TGFbBMPgene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_TGFbBMPgene_type_15, "pathvar_TGFbBMPgene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_type_15.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_type_15.csv')









# TGFbBMPgene Age -------------------------------------------------------------

# make sure the execute the code in the 'TGFbBMPgenes' section

### linear regression
# "TGFbBMPgene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ TGFbBMPgene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_TGFbBMPgene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_TGFbBMPgene_age_1, "pathvar_TGFbBMPgene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_age_1.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ SLC39A13+ADAMTSL2+SOX9+DDRGK1+LEMD3+BMP2+BMPR1B+GDF5+BMP1+LTBP3, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_TGFbBMPgene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_TGFbBMPgene_age_2, "pathvar_TGFbBMPgene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_age_2.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "TGFbBMPgene", 'SLC39A13', 'ADAMTSL2', 'SOX9', 'DDRGK1', 
                                     'LEMD3', 'BMP2', 'BMPR1B', 'GDF5', 'BMP1', 'LTBP3')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     17

with(data_setD, table(LTBP3))
# all good


data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "TGFbBMPgene", 'SLC39A13', 'ADAMTSL2', 'SOX9', 
                                     'DDRGK1', 'LEMD3', 'BMP2', 'BMPR1B', 'GDF5', 'BMP1', 'LTBP3')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    22

with(data_setE, table(LTBP3))
# all good


### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ TGFbBMPgene, data = data_setD)
# summary(LM_15)
pathvar_TGFbBMPgene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_TGFbBMPgene_age_3, "pathvar_TGFbBMPgene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_age_3.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+TGFbBMPgene, data = data_setD)
# summary(LM_16)
pathvar_TGFbBMPgene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_TGFbBMPgene_age_4, "pathvar_TGFbBMPgene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_age_4.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ SLC39A13+ADAMTSL2+SOX9+DDRGK1+LEMD3+BMP2+BMPR1B+GDF5+BMP1+LTBP3, data = data_setD)
# summary(LM_17)
pathvar_TGFbBMPgene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_TGFbBMPgene_age_5, "pathvar_TGFbBMPgene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_age_5.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SLC39A13+ADAMTSL2+SOX9+DDRGK1+LEMD3+BMP2+BMPR1B+GDF5+BMP1+LTBP3, data = data_setD)
# summary(LM_18)
pathvar_TGFbBMPgene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_TGFbBMPgene_age_6, "pathvar_TGFbBMPgene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_age_6.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ TGFbBMPgene, data = data_setE)
# summary(LM_15)
pathvar_TGFbBMPgene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_TGFbBMPgene_age_7, "pathvar_TGFbBMPgene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_age_7.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+TGFbBMPgene, data = data_setE)
# summary(LM_16)
pathvar_TGFbBMPgene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_TGFbBMPgene_age_8, "pathvar_TGFbBMPgene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_age_8.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ SLC39A13+ADAMTSL2+SOX9+DDRGK1+LEMD3+BMP2+BMPR1B+GDF5+BMP1+LTBP3, data = data_setE)
# summary(LM_17)
pathvar_TGFbBMPgene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_TGFbBMPgene_age_9, "pathvar_TGFbBMPgene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_age_9.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+SLC39A13+ADAMTSL2+SOX9+DDRGK1+LEMD3+BMP2+BMPR1B+GDF5+BMP1+LTBP3, data = data_setE)
# summary(LM_18)
pathvar_TGFbBMPgene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_TGFbBMPgene_age_10, "pathvar_TGFbBMPgene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_TGFbBMPgene_age_10.csv --path Emily-folder/results_2/pathvar_TGFbBMPgene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(df_cancer_diagnosis, table(TGFbBMPgene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(TGFbBMPgene, Ethnicity_White)))
# p-value = 0.1747

with(df_cancer_diagnosis, table(TGFbBMPgene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(TGFbBMPgene, Ever_smoked_dummy)))
# p-value = 0.9697

with(df_cancer_diagnosis, table(TGFbBMPgene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(TGFbBMPgene, Sex_Male)))
# p-value = 0.647

with(df_cancer_diagnosis, table(TGFbBMPgene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(TGFbBMPgene, University_education)))
# p-value = 0.6312

with(df_cancer_diagnosis, table(TGFbBMPgene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(TGFbBMPgene, Annual_income)))
# p-value = 0.1683


### data_setD
with(data_setD, table(TGFbBMPgene, Ethnicity_White))			
fisher.test(with(data_setD, table(TGFbBMPgene, Ethnicity_White)))
# p-value = 0.173

with(data_setD, table(TGFbBMPgene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(TGFbBMPgene, Ever_smoked_dummy)))
# p-value = 1

with(data_setD, table(TGFbBMPgene, Sex_Male))			
fisher.test(with(data_setD, table(TGFbBMPgene, Sex_Male)))
# p-value = 0.6857


### data_setE
with(data_setE, table(TGFbBMPgene, Ethnicity_White))			
fisher.test(with(data_setE, table(TGFbBMPgene, Ethnicity_White)))
# p-value = 0.7286

with(data_setE, table(TGFbBMPgene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(TGFbBMPgene, Ever_smoked_dummy)))
# p-value = 0.6388

with(data_setE, table(TGFbBMPgene, Sex_Male))			
fisher.test(with(data_setE, table(TGFbBMPgene, Sex_Male)))
# p-value = 0.6366

with(data_setE, table(TGFbBMPgene, University_education))			
fisher.test(with(data_setE, table(TGFbBMPgene, University_education)))
# p-value = 0.4995

with(data_setE, table(TGFbBMPgene, Annual_income))			
fisher.test(with(data_setE, table(TGFbBMPgene, Annual_income)))
# p-value = 0.202



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 

# make TGFbBMPgene into factor
df_cancer_diagnosis$TGFbBMPgene <- as.factor(df_cancer_diagnosis$TGFbBMPgene)

# try a violin plot



# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ TGFbBMPgene, data = df_cancer_diagnosis)
# 



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ TGFbBMPgene, data = df_cancer_diagnosis)
# W = 169984098, p-value = 0.04178



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ TGFbBMPgene, data = df_cancer_diagnosis)
# W = 164881403, p-value = 0.9436



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ TGFbBMPgene, data = df_cancer_diagnosis)
# W = 168787996, p-value = 0.04596



# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(TGFbBMPgene), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$TGFbBMPgene <- as.factor(data_setD$TGFbBMPgene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ TGFbBMPgene, data = data_setD)
# 



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ TGFbBMPgene, data = data_setD)
# W = 167187212, p-value = 0.04236




### data_setE

data_setE$TGFbBMPgene <- as.factor(data_setE$TGFbBMPgene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ TGFbBMPgene, data = data_setE)
# 



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ TGFbBMPgene, data = data_setE)
# W = 107643266, p-value = 0.03334



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ TGFbBMPgene, data = data_setE)
# W = 104309249, p-value = 0.5797



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = TGFbBMPgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ TGFbBMPgene, data = data_setE)
# W = 107807104, p-value = 0.02393


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(TGFbBMPgene), scales = "free_y") +
  theme_light()
# 










# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_TGFbBMPgenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_TGFbBMPgenes.R')
