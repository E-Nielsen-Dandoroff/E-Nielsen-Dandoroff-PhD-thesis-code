# 09/07/2024

# EntireCohort_pathvar_Braingenes.R

# this script is for the Brain associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)


# Braingene
# NUP133
# DYNC1I2
# CASK
# SLC1A4
# MFSD2A
# PPFIBP1
# PCLO
# PCDH12
# OCLN
# FRMD4A
# NUP37








# Braingenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')


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

















# Braingene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_Braingene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_Braingene_, "pathvar_Braingene_.csv", row.names=FALSE)
# system('dx upload pathvar_Braingene_.csv --path Emily-folder/results_2/pathvar_Braingene_.csv')


# Braingene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ Braingene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_Braingene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_incidence_1, "pathvar_Braingene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_incidence_1.csv --path Emily-folder/results_2/pathvar_Braingene_incidence_1.csv')


# GOIs
LR_2 <- glm(formula = cancer_ROOC_dummy ~ NUP133+DYNC1I2+CASK+SLC1A4+MFSD2A+PPFIBP1+PCLO+PCDH12+OCLN+FRMD4A+NUP37, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_Braingene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_incidence_2, "pathvar_Braingene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_incidence_2.csv --path Emily-folder/results_2/pathvar_Braingene_incidence_2.csv')







### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "Braingene", 'NUP133', 'DYNC1I2', 'CASK', 'SLC1A4', 
                                     'MFSD2A', 'PPFIBP1', 'PCLO', 'PCDH12', 'OCLN', 'FRMD4A', 
                                     'NUP37')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     18

with(data_setA, table(NUP37))



# Braingene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ Braingene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_Braingene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_incidence_3, "pathvar_Braingene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_incidence_3.csv --path Emily-folder/results_2/pathvar_Braingene_incidence_3.csv')


# age + ethnicity + smoking + sex + Braingene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_Braingene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_incidence_4, "pathvar_Braingene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_incidence_4.csv --path Emily-folder/results_2/pathvar_Braingene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ NUP133+DYNC1I2+CASK+SLC1A4+MFSD2A+PPFIBP1+PCLO+PCDH12+OCLN+FRMD4A+NUP37, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_Braingene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_incidence_5, "pathvar_Braingene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_incidence_5.csv --path Emily-folder/results_2/pathvar_Braingene_incidence_5.csv')


# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+NUP133+DYNC1I2+CASK+SLC1A4+MFSD2A+PPFIBP1+PCLO+PCDH12+OCLN+FRMD4A+NUP37, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_Braingene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_incidence_6, "pathvar_Braingene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_incidence_6.csv --path Emily-folder/results_2/pathvar_Braingene_incidence_6.csv')




### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "Braingene", 'NUP133', 'DYNC1I2', 'CASK', 'SLC1A4', 
                                     'MFSD2A', 'PPFIBP1', 'PCLO', 'PCDH12', 'OCLN', 'FRMD4A', 
                                     'NUP37')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     23

with(data_setB, table(NUP37))


# Braingene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ Braingene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_Braingene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_incidence_7, "pathvar_Braingene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_incidence_7.csv --path Emily-folder/results_2/pathvar_Braingene_incidence_7.csv')


# all variables + Braingene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+Braingene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_Braingene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_incidence_8, "pathvar_Braingene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_incidence_8.csv --path Emily-folder/results_2/pathvar_Braingene_incidence_8.csv')


# GOIs alone
LR_9 <- glm(formula = cancer_ROOC_dummy ~ NUP133+DYNC1I2+CASK+SLC1A4+MFSD2A+PPFIBP1+PCLO+PCDH12+OCLN+FRMD4A+NUP37, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_Braingene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_incidence_9, "pathvar_Braingene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_incidence_9.csv --path Emily-folder/results_2/pathvar_Braingene_incidence_9.csv')


# all variables + GOIs
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+NUP133+DYNC1I2+CASK+SLC1A4+MFSD2A+PPFIBP1+PCLO+PCDH12+OCLN+FRMD4A+NUP37, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_Braingene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_incidence_10, "pathvar_Braingene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_incidence_10.csv --path Emily-folder/results_2/pathvar_Braingene_incidence_10.csv')










### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(Braingene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(Braingene, cancer_ROOC_dummy)))
# p-value = 0.5185

with(EntireCohort_pathvar, table(Braingene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(Braingene, Ethnicity_White)))
# p-value = 0.001053
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.039910 1.172778
# sample estimates:
#   odds ratio 
# 1.103946 

with(EntireCohort_pathvar, table(Braingene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(Braingene, Ever_smoked_dummy)))
# p-value = 0.9893

with(EntireCohort_pathvar, table(Braingene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(Braingene, Sex_Male)))
# p-value = 0.4956

with(EntireCohort_pathvar, table(Braingene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(Braingene, University_education)))
# p-value = 0.0194
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9411610 0.9947234
# sample estimates:
#   odds ratio 
# 0.9675927

with(EntireCohort_pathvar, table(Braingene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(Braingene, Annual_income))	)
# p-value = 0.9035

with(EntireCohort_pathvar, table(Braingene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(Braingene, Smoker_previous)))
# p-value = 0.9397

with(EntireCohort_pathvar, table(Braingene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(Braingene, Smoker_current))	)
# p-value = 0.3118



### data_setA
with(data_setA, table(Braingene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(Braingene, cancer_ROOC_dummy)))
# p-value = 0.5626

with(data_setA, table(Braingene, Ethnicity_White))			
fisher.test(with(data_setA, table(Braingene, Ethnicity_White)))
# p-value = 0.001119
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.039585 1.172935
# sample estimates:
#   odds ratio 
# 1.103845

with(data_setA, table(Braingene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(Braingene, Ever_smoked_dummy)))
# p-value = 0.9947

with(data_setA, table(Braingene, Sex_Male))			
fisher.test(with(data_setA, table(Braingene, Sex_Male)))
# p-value = 0.498



### data_setB
with(data_setB, table(Braingene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(Braingene, cancer_ROOC_dummy)))
# p-value = 0.2602

with(data_setB, table(Braingene, Ethnicity_White))			
fisher.test(with(data_setB, table(Braingene, Ethnicity_White)))
# p-value = 0.04989
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.000303 1.155203
# sample estimates:
#   odds ratio 
# 1.074377

with(data_setB, table(Braingene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(Braingene, Ever_smoked_dummy)))
# p-value = 0.4491

with(data_setB, table(Braingene, Sex_Male))			
fisher.test(with(data_setB, table(Braingene, Sex_Male)))
# p-value = 0.2956

with(data_setB, table(Braingene, University_education))			
fisher.test(with(data_setB, table(Braingene, University_education)))
# p-value = 0.002723
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9271054 0.9843801
# sample estimates:
#   odds ratio 
# 0.9553413

with(data_setB, table(Braingene, Annual_income))			
fisher.test(with(data_setB, table(Braingene, Annual_income)))
# p-value = 0.9869






### can plot the continuous variables: 

# make Braingene into factor
df_for_plots <- EntireCohort_pathvar
df_for_plots$Braingene <- as.factor(df_for_plots$Braingene)

# try a violin plot

# Age_at_recruitment  
df_for_plots %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ Braingene, data = df_for_plots)
# W = 5508065682, p-value = 0.8776

# BMI  
df_for_plots %>% 
  ggplot(mapping = aes(x = Braingene, y = BMI)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = Braingene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ Braingene, data = df_for_plots)
# W = 5402790072, p-value = 0.002577

# Standing_height 
df_for_plots %>% 
  ggplot(mapping = aes(x = Braingene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = Braingene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ Braingene, data = df_for_plots)
# W = 5466482080, p-value = 0.5825

# make this a bar plot
# Weekly_exercise
df_for_plots %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(Braingene), scales = "free_y") +
  theme_light()
# 



### data_setA

# make Braingene into factor
data_setA$Braingene <- as.factor(data_setA$Braingene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ Braingene, data = data_setA)
# W = 5415834062, p-value = 0.8134



### data_setB

# make Braingene into factor
data_setB$Braingene <- as.factor(data_setB$Braingene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ Braingene, data = data_setB)
# W = 3620385624, p-value = 0.8786

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = Braingene, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = Braingene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ Braingene, data = data_setB)
# W = 3565420591, p-value = 0.0005149

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = Braingene, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = Braingene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ Braingene, data = data_setB)
# W = 3602156062, p-value = 0.2936


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(Braingene), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# Braingene Type ------------------------------------------------------------

# make sure the execute the code in the 'Braingenes' section

## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_Braingene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_type_1, "pathvar_Braingene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_type_1.csv --path Emily-folder/results_2/pathvar_Braingene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_Braingene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_type_2, "pathvar_Braingene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_type_2.csv --path Emily-folder/results_2/pathvar_Braingene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_Braingene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_type_3, "pathvar_Braingene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_type_3.csv --path Emily-folder/results_2/pathvar_Braingene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_Braingene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_type_4, "pathvar_Braingene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_type_4.csv --path Emily-folder/results_2/pathvar_Braingene_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_Braingene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_type_5, "pathvar_Braingene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_type_5.csv --path Emily-folder/results_2/pathvar_Braingene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_Braingene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_type_6, "pathvar_Braingene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_type_6.csv --path Emily-folder/results_2/pathvar_Braingene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_Braingene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_type_7, "pathvar_Braingene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_type_7.csv --path Emily-folder/results_2/pathvar_Braingene_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_Braingene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_type_8, "pathvar_Braingene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_type_8.csv --path Emily-folder/results_2/pathvar_Braingene_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_Braingene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_type_9, "pathvar_Braingene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_type_9.csv --path Emily-folder/results_2/pathvar_Braingene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_Braingene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_type_10, "pathvar_Braingene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_type_10.csv --path Emily-folder/results_2/pathvar_Braingene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_Braingene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_type_11, "pathvar_Braingene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_type_11.csv --path Emily-folder/results_2/pathvar_Braingene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_Braingene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_type_12, "pathvar_Braingene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_type_12.csv --path Emily-folder/results_2/pathvar_Braingene_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_Braingene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_type_13, "pathvar_Braingene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_type_13.csv --path Emily-folder/results_2/pathvar_Braingene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_Braingene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_type_14, "pathvar_Braingene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_type_14.csv --path Emily-folder/results_2/pathvar_Braingene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_Braingene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Braingene_type_15, "pathvar_Braingene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_type_15.csv --path Emily-folder/results_2/pathvar_Braingene_type_15.csv')








# Braingene Age -------------------------------------------------------------

# make sure the execute the code in the 'Braingenes' section

### linear regression
# "Braingene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ Braingene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_Braingene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_Braingene_age_1, "pathvar_Braingene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_age_1.csv --path Emily-folder/results_2/pathvar_Braingene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ NUP133+DYNC1I2+CASK+SLC1A4+MFSD2A+PPFIBP1+PCLO+PCDH12+OCLN+FRMD4A+NUP37, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_Braingene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_Braingene_age_2, "pathvar_Braingene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_age_2.csv --path Emily-folder/results_2/pathvar_Braingene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "Braingene", 'NUP133', 'DYNC1I2', 'CASK', 'SLC1A4', 'MFSD2A', 
                                     'PPFIBP1', 'PCLO', 'PCDH12', 'OCLN', 'FRMD4A', 'NUP37')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     18

with(data_setD, table(NUP37))
# CASK = 75
# NUP37 = 82



data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "Braingene", 'NUP133', 'DYNC1I2', 'CASK', 
                                     'SLC1A4', 'MFSD2A', 'PPFIBP1', 'PCLO', 'PCDH12', 'OCLN', 'FRMD4A', 
                                     'NUP37')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    23

with(data_setE, table(NUP37))
# CASK = 46
# NUP37 = 70




### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ Braingene, data = data_setD)
# summary(LM_15)
pathvar_Braingene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_Braingene_age_3, "pathvar_Braingene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_age_3.csv --path Emily-folder/results_2/pathvar_Braingene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Braingene, data = data_setD)
# summary(LM_16)
pathvar_Braingene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_Braingene_age_4, "pathvar_Braingene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_age_4.csv --path Emily-folder/results_2/pathvar_Braingene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ NUP133+DYNC1I2+CASK+SLC1A4+MFSD2A+PPFIBP1+PCLO+PCDH12+OCLN+FRMD4A+NUP37, data = data_setD)
# summary(LM_17)
pathvar_Braingene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_Braingene_age_5, "pathvar_Braingene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_age_5.csv --path Emily-folder/results_2/pathvar_Braingene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+NUP133+DYNC1I2+CASK+SLC1A4+MFSD2A+PPFIBP1+PCLO+PCDH12+OCLN+FRMD4A+NUP37, data = data_setD)
# summary(LM_18)
pathvar_Braingene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_Braingene_age_6, "pathvar_Braingene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_age_6.csv --path Emily-folder/results_2/pathvar_Braingene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ Braingene, data = data_setE)
# summary(LM_15)
pathvar_Braingene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_Braingene_age_7, "pathvar_Braingene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_age_7.csv --path Emily-folder/results_2/pathvar_Braingene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+Braingene, data = data_setE)
# summary(LM_16)
pathvar_Braingene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_Braingene_age_8, "pathvar_Braingene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_age_8.csv --path Emily-folder/results_2/pathvar_Braingene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ NUP133+DYNC1I2+CASK+SLC1A4+MFSD2A+PPFIBP1+PCLO+PCDH12+OCLN+FRMD4A+NUP37, data = data_setE)
# summary(LM_17)
pathvar_Braingene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_Braingene_age_9, "pathvar_Braingene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_age_9.csv --path Emily-folder/results_2/pathvar_Braingene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+NUP133+DYNC1I2+CASK+SLC1A4+MFSD2A+PPFIBP1+PCLO+PCDH12+OCLN+FRMD4A+NUP37, data = data_setE)
# summary(LM_18)
pathvar_Braingene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_Braingene_age_10, "pathvar_Braingene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_Braingene_age_10.csv --path Emily-folder/results_2/pathvar_Braingene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(df_cancer_diagnosis, table(Braingene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(Braingene, Ethnicity_White)))
# p-value = 0.7384

with(df_cancer_diagnosis, table(Braingene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(Braingene, Ever_smoked_dummy)))
# p-value = 0.5373

with(df_cancer_diagnosis, table(Braingene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(Braingene, Sex_Male)))
# p-value = 0.3504

with(df_cancer_diagnosis, table(Braingene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(Braingene, University_education)))
# p-value = 0.584

with(df_cancer_diagnosis, table(Braingene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(Braingene, Annual_income)))
# p-value = 0.1622


### data_setD
with(data_setD, table(Braingene, Ethnicity_White))			
fisher.test(with(data_setD, table(Braingene, Ethnicity_White)))
# p-value = 0.8339

with(data_setD, table(Braingene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(Braingene, Ever_smoked_dummy)))
# p-value = 0.5461

with(data_setD, table(Braingene, Sex_Male))			
fisher.test(with(data_setD, table(Braingene, Sex_Male)))
# p-value = 0.3147


### data_setE
with(data_setE, table(Braingene, Ethnicity_White))			
fisher.test(with(data_setE, table(Braingene, Ethnicity_White)))
# p-value = 0.6082

with(data_setE, table(Braingene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(Braingene, Ever_smoked_dummy)))
# p-value = 0.8376

with(data_setE, table(Braingene, Sex_Male))			
fisher.test(with(data_setE, table(Braingene, Sex_Male)))
# p-value = 0.2077

with(data_setE, table(Braingene, University_education))			
fisher.test(with(data_setE, table(Braingene, University_education)))
# p-value = 0.6407

with(data_setE, table(Braingene, Annual_income))			
fisher.test(with(data_setE, table(Braingene, Annual_income)))
# p-value = 0.1656



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 

# make Braingene into factor
df_cancer_diagnosis$Braingene <- as.factor(df_cancer_diagnosis$Braingene)

# try a violin plot



# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ Braingene, data = df_cancer_diagnosis)
# 



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ Braingene, data = df_cancer_diagnosis)
# W = 306737236, p-value = 0.4279



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Braingene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Braingene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ Braingene, data = df_cancer_diagnosis)
# W = 300751700, p-value = 0.4887



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Braingene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Braingene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ Braingene, data = df_cancer_diagnosis)
# W = 300106805, p-value = 0.1976



# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(Braingene), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$Braingene <- as.factor(data_setD$Braingene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ Braingene, data = data_setD)
# 



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ Braingene, data = data_setD)
# W = 301460618, p-value = 0.5265




### data_setE

data_setE$Braingene <- as.factor(data_setE$Braingene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ Braingene, data = data_setE)
# 



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = Braingene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ Braingene, data = data_setE)
# W = 192346626, p-value = 0.9809



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = Braingene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = Braingene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ Braingene, data = data_setE)
# W = 190956256, p-value = 0.4227



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = Braingene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = Braingene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ Braingene, data = data_setE)
# W = 189951887, p-value = 0.1619


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(Braingene), scales = "free_y") +
  theme_light()
# 














# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_Braingenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_Braingenes.R')
