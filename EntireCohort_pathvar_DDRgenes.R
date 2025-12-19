# 24/06/2024

# EntireCohort_pathvar_DDRgenes.R

# this script is for the DNA damage repair associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)


# DDRgene
# PNKP
# RAD50
# TONSL
# TOP3A
# RBBP8
# TRAIP
# NHEJ1
# NSMCE2
# ERCC8
# ERCC6












# DDRgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')


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

















# DDRgene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_DDRgene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_DDRgene_, "pathvar_DDRgene_.csv", row.names=FALSE)
# system('dx upload pathvar_DDRgene_.csv --path Emily-folder/results_2/pathvar_DDRgene_.csv')


# DDRgene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ DDRgene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_DDRgene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_incidence_1, "pathvar_DDRgene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_incidence_1.csv --path Emily-folder/results_2/pathvar_DDRgene_incidence_1.csv')


# GOIs
LR_2 <- glm(formula = cancer_ROOC_dummy ~ PNKP+RAD50+TONSL+TOP3A+RBBP8+TRAIP+NHEJ1+NSMCE2+ERCC8+ERCC6, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_DDRgene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_incidence_2, "pathvar_DDRgene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_incidence_2.csv --path Emily-folder/results_2/pathvar_DDRgene_incidence_2.csv')







### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "DDRgene", 'PNKP', 'RAD50', 'TONSL', 'TOP3A', 
                                     'RBBP8', 'TRAIP', 'NHEJ1', 'NSMCE2', 'ERCC8', 'ERCC6')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     17

with(data_setA, table(ERCC6))



# DDRgene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ DDRgene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_DDRgene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_incidence_3, "pathvar_DDRgene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_incidence_3.csv --path Emily-folder/results_2/pathvar_DDRgene_incidence_3.csv')


# age + ethnicity + smoking + sex + DDRgene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_DDRgene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_incidence_4, "pathvar_DDRgene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_incidence_4.csv --path Emily-folder/results_2/pathvar_DDRgene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ PNKP+RAD50+TONSL+TOP3A+RBBP8+TRAIP+NHEJ1+NSMCE2+ERCC8+ERCC6, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_DDRgene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_incidence_5, "pathvar_DDRgene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_incidence_5.csv --path Emily-folder/results_2/pathvar_DDRgene_incidence_5.csv')


# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+PNKP+RAD50+TONSL+TOP3A+RBBP8+TRAIP+NHEJ1+NSMCE2+ERCC8+ERCC6, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_DDRgene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_incidence_6, "pathvar_DDRgene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_incidence_6.csv --path Emily-folder/results_2/pathvar_DDRgene_incidence_6.csv')




### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "DDRgene", 'PNKP', 'RAD50', 'TONSL', 'TOP3A', 
                                     'RBBP8', 'TRAIP', 'NHEJ1', 'NSMCE2', 'ERCC8', 'ERCC6')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     22

with(data_setB, table(ERCC6))




# DDRgene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ DDRgene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_DDRgene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_incidence_7, "pathvar_DDRgene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_incidence_7.csv --path Emily-folder/results_2/pathvar_DDRgene_incidence_7.csv')


# all variables + DDRgene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+DDRgene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_DDRgene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_incidence_8, "pathvar_DDRgene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_incidence_8.csv --path Emily-folder/results_2/pathvar_DDRgene_incidence_8.csv')


# GOIs alone
LR_9 <- glm(formula = cancer_ROOC_dummy ~ PNKP+RAD50+TONSL+TOP3A+RBBP8+TRAIP+NHEJ1+NSMCE2+ERCC8+ERCC6, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_DDRgene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_incidence_9, "pathvar_DDRgene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_incidence_9.csv --path Emily-folder/results_2/pathvar_DDRgene_incidence_9.csv')


# all variables + GOIs
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+PNKP+RAD50+TONSL+TOP3A+RBBP8+TRAIP+NHEJ1+NSMCE2+ERCC8+ERCC6, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_DDRgene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_incidence_10, "pathvar_DDRgene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_incidence_10.csv --path Emily-folder/results_2/pathvar_DDRgene_incidence_10.csv')










### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(DDRgene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(DDRgene, cancer_ROOC_dummy)))
# p-value = 0.6277

with(EntireCohort_pathvar, table(DDRgene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(DDRgene, Ethnicity_White)))
# p-value = 2.097e-15
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.180037 1.322094
# sample estimates:
#   odds ratio 
# 1.248577

with(EntireCohort_pathvar, table(DDRgene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(DDRgene, Ever_smoked_dummy)))
# p-value = 0.838

with(EntireCohort_pathvar, table(DDRgene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(DDRgene, Sex_Male)))
# p-value = 0.194

with(EntireCohort_pathvar, table(DDRgene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(DDRgene, University_education)))
# p-value = 0.2383

with(EntireCohort_pathvar, table(DDRgene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(DDRgene, Annual_income))	)
# p-value = 0.7878

with(EntireCohort_pathvar, table(DDRgene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(DDRgene, Smoker_previous)))
# p-value = 0.3926

with(EntireCohort_pathvar, table(DDRgene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(DDRgene, Smoker_current))	)
# p-value = 0.737



### data_setA
with(data_setA, table(DDRgene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(DDRgene, cancer_ROOC_dummy)))
# p-value = 0.6815

with(data_setA, table(DDRgene, Ethnicity_White))			
fisher.test(with(data_setA, table(DDRgene, Ethnicity_White)))
# p-value = 1.451e-15
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.182557 1.325594
# sample estimates:
#   odds ratio 
# 1.251582 

with(data_setA, table(DDRgene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(DDRgene, Ever_smoked_dummy)))
# p-value = 0.8708

with(data_setA, table(DDRgene, Sex_Male))			
fisher.test(with(data_setA, table(DDRgene, Sex_Male)))
# p-value = 0.2045



### data_setB
with(data_setB, table(DDRgene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(DDRgene, cancer_ROOC_dummy)))
# p-value = 0.6994

with(data_setB, table(DDRgene, Ethnicity_White))			
fisher.test(with(data_setB, table(DDRgene, Ethnicity_White)))
# p-value = 6.285e-13
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.191834 1.368986
# sample estimates:
#   odds ratio 
# 1.276756

with(data_setB, table(DDRgene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(DDRgene, Ever_smoked_dummy)))
# p-value = 0.3541

with(data_setB, table(DDRgene, Sex_Male))			
fisher.test(with(data_setB, table(DDRgene, Sex_Male)))
# p-value = 0.1831

with(data_setB, table(DDRgene, University_education))			
fisher.test(with(data_setB, table(DDRgene, University_education)))
# p-value = 0.221

with(data_setB, table(DDRgene, Annual_income))			
fisher.test(with(data_setB, table(DDRgene, Annual_income)))
# p-value = 0.5694






### can plot the continuous variables: 

# make DDRgene into factor
df_for_plots <- EntireCohort_pathvar
df_for_plots$DDRgene <- as.factor(df_for_plots$DDRgene)

# try a violin plot

# Age_at_recruitment  
df_for_plots %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ DDRgene, data = df_for_plots)
# W = 6806722114, p-value = 0.4045

# BMI  
df_for_plots %>% 
  ggplot(mapping = aes(x = DDRgene, y = BMI)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = DDRgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ DDRgene, data = df_for_plots)
# W = 6760835656, p-value = 0.2846

# Standing_height 
df_for_plots %>% 
  ggplot(mapping = aes(x = DDRgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = DDRgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ DDRgene, data = df_for_plots)
# W = 6750427579, p-value = 0.9896

# make this a bar plot
# Weekly_exercise
df_for_plots %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(DDRgene), scales = "free_y") +
  theme_light()
# 



### data_setA

# make DDRgene into factor
data_setA$DDRgene <- as.factor(data_setA$DDRgene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ DDRgene, data = data_setA)
# W = 6695061801, p-value = 0.3918



### data_setB

# make DDRgene into factor
data_setB$DDRgene <- as.factor(data_setB$DDRgene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ DDRgene, data = data_setB)
# W = 4.501e+09, p-value = 0.2848

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = DDRgene, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = DDRgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ DDRgene, data = data_setB)
# W = 4491204980, p-value = 0.624

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = DDRgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = DDRgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ DDRgene, data = data_setB)
# W = 4482771962, p-value = 0.9924


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(DDRgene), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# DDRgene Type ------------------------------------------------------------

# make sure the execute the code in the 'DDRgenes' section


## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_DDRgene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_type_1, "pathvar_DDRgene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_type_1.csv --path Emily-folder/results_2/pathvar_DDRgene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_DDRgene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_type_2, "pathvar_DDRgene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_type_2.csv --path Emily-folder/results_2/pathvar_DDRgene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_DDRgene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_type_3, "pathvar_DDRgene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_type_3.csv --path Emily-folder/results_2/pathvar_DDRgene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_DDRgene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_type_4, "pathvar_DDRgene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_type_4.csv --path Emily-folder/results_2/pathvar_DDRgene_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_DDRgene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_type_5, "pathvar_DDRgene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_type_5.csv --path Emily-folder/results_2/pathvar_DDRgene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_DDRgene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_type_6, "pathvar_DDRgene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_type_6.csv --path Emily-folder/results_2/pathvar_DDRgene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_DDRgene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_type_7, "pathvar_DDRgene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_type_7.csv --path Emily-folder/results_2/pathvar_DDRgene_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_DDRgene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_type_8, "pathvar_DDRgene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_type_8.csv --path Emily-folder/results_2/pathvar_DDRgene_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_DDRgene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_type_9, "pathvar_DDRgene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_type_9.csv --path Emily-folder/results_2/pathvar_DDRgene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_DDRgene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_type_10, "pathvar_DDRgene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_type_10.csv --path Emily-folder/results_2/pathvar_DDRgene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_DDRgene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_type_11, "pathvar_DDRgene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_type_11.csv --path Emily-folder/results_2/pathvar_DDRgene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_DDRgene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_type_12, "pathvar_DDRgene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_type_12.csv --path Emily-folder/results_2/pathvar_DDRgene_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_DDRgene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_type_13, "pathvar_DDRgene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_type_13.csv --path Emily-folder/results_2/pathvar_DDRgene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_DDRgene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_type_14, "pathvar_DDRgene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_type_14.csv --path Emily-folder/results_2/pathvar_DDRgene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_DDRgene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_DDRgene_type_15, "pathvar_DDRgene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_type_15.csv --path Emily-folder/results_2/pathvar_DDRgene_type_15.csv')







# DDRgene Age -------------------------------------------------------------

# make sure the execute the code in the 'DDRgenes' section

### linear regression
# "DDRgene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ DDRgene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_DDRgene_age_1 <- tidy(LM_1, conf.int = TRUE,)
write.csv(pathvar_DDRgene_age_1, "pathvar_DDRgene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_age_1.csv --path Emily-folder/results_2/pathvar_DDRgene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ PNKP+RAD50+TONSL+TOP3A+RBBP8+TRAIP+NHEJ1+NSMCE2+ERCC8+ERCC6, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_DDRgene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_DDRgene_age_2, "pathvar_DDRgene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_age_2.csv --path Emily-folder/results_2/pathvar_DDRgene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "DDRgene", 'PNKP', 'RAD50', 'TONSL', 'TOP3A', 'RBBP8', 
                                     'TRAIP', 'NHEJ1', 'NSMCE2', 'ERCC8', 'ERCC6')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     17

with(data_setD, table(ERCC6))
# NSMCE2 = 29



data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "DDRgene", 'PNKP', 'RAD50', 'TONSL', 
                                     'TOP3A', 'RBBP8', 'TRAIP', 'NHEJ1', 'NSMCE2', 'ERCC8', 'ERCC6')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    22

with(data_setE, table(ERCC6))
# NSMCE2 = 19




### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ DDRgene, data = data_setD)
# summary(LM_15)
pathvar_DDRgene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_DDRgene_age_3, "pathvar_DDRgene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_age_3.csv --path Emily-folder/results_2/pathvar_DDRgene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+DDRgene, data = data_setD)
# summary(LM_16)
pathvar_DDRgene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_DDRgene_age_4, "pathvar_DDRgene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_age_4.csv --path Emily-folder/results_2/pathvar_DDRgene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ PNKP+RAD50+TONSL+TOP3A+RBBP8+TRAIP+NHEJ1+NSMCE2+ERCC8+ERCC6, data = data_setD)
# summary(LM_17)
pathvar_DDRgene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_DDRgene_age_5, "pathvar_DDRgene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_age_5.csv --path Emily-folder/results_2/pathvar_DDRgene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+PNKP+RAD50+TONSL+TOP3A+RBBP8+TRAIP+NHEJ1+NSMCE2+ERCC8+ERCC6, data = data_setD)
# summary(LM_18)
pathvar_DDRgene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_DDRgene_age_6, "pathvar_DDRgene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_age_6.csv --path Emily-folder/results_2/pathvar_DDRgene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ DDRgene, data = data_setE)
# summary(LM_15)
pathvar_DDRgene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_DDRgene_age_7, "pathvar_DDRgene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_age_7.csv --path Emily-folder/results_2/pathvar_DDRgene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+DDRgene, data = data_setE)
# summary(LM_16)
pathvar_DDRgene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_DDRgene_age_8, "pathvar_DDRgene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_age_8.csv --path Emily-folder/results_2/pathvar_DDRgene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ PNKP+RAD50+TONSL+TOP3A+RBBP8+TRAIP+NHEJ1+NSMCE2+ERCC8+ERCC6, data = data_setE)
# summary(LM_17)
pathvar_DDRgene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_DDRgene_age_9, "pathvar_DDRgene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_age_9.csv --path Emily-folder/results_2/pathvar_DDRgene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+PNKP+RAD50+TONSL+TOP3A+RBBP8+TRAIP+NHEJ1+NSMCE2+ERCC8+ERCC6, data = data_setE)
# summary(LM_18)
pathvar_DDRgene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_DDRgene_age_10, "pathvar_DDRgene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_DDRgene_age_10.csv --path Emily-folder/results_2/pathvar_DDRgene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 110850


with(df_cancer_diagnosis, table(DDRgene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(DDRgene, Ethnicity_White)))
# p-value = 0.03288
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.013673 1.398099
# sample estimates:
#   odds ratio 
# 1.187141

with(df_cancer_diagnosis, table(DDRgene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(DDRgene, Ever_smoked_dummy)))
# p-value = 1

with(df_cancer_diagnosis, table(DDRgene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(DDRgene, Sex_Male)))
# p-value = 0.1632

with(df_cancer_diagnosis, table(DDRgene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(DDRgene, University_education)))
# p-value = 0.7203

with(df_cancer_diagnosis, table(DDRgene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(DDRgene, Annual_income)))
# p-value = 0.7087


### data_setD
with(data_setD, table(DDRgene, Ethnicity_White))			
fisher.test(with(data_setD, table(DDRgene, Ethnicity_White)))
# p-value = 0.02646
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.021422 1.412755
# sample estimates:
#   odds ratio 
# 1.197843 

with(data_setD, table(DDRgene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(DDRgene, Ever_smoked_dummy)))
# p-value = 0.9498

with(data_setD, table(DDRgene, Sex_Male))			
fisher.test(with(data_setD, table(DDRgene, Sex_Male)))
# p-value = 0.1802


### data_setE
with(data_setE, table(DDRgene, Ethnicity_White))			
fisher.test(with(data_setE, table(DDRgene, Ethnicity_White)))
# p-value = 0.06129

with(data_setE, table(DDRgene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(DDRgene, Ever_smoked_dummy)))
# p-value = 0.5094

with(data_setE, table(DDRgene, Sex_Male))			
fisher.test(with(data_setE, table(DDRgene, Sex_Male)))
# p-value = 0.2186

with(data_setE, table(DDRgene, University_education))			
fisher.test(with(data_setE, table(DDRgene, University_education)))
# p-value = 0.9771

with(data_setE, table(DDRgene, Annual_income))			
fisher.test(with(data_setE, table(DDRgene, Annual_income)))
# p-value = 0.4893



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 

# make DDRgene into factor
df_cancer_diagnosis$DDRgene <- as.factor(df_cancer_diagnosis$DDRgene)

# try a violin plot



# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ DDRgene, data = df_cancer_diagnosis)
# 



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ DDRgene, data = df_cancer_diagnosis)
# W = 381840899, p-value = 0.4232



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = DDRgene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = DDRgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ DDRgene, data = df_cancer_diagnosis)
# W = 377898258, p-value = 0.8523



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = DDRgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = DDRgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ DDRgene, data = df_cancer_diagnosis)
# W = 378177690, p-value = 0.9897



# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(DDRgene), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$DDRgene <- as.factor(data_setD$DDRgene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ DDRgene, data = data_setD)
# 



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ DDRgene, data = data_setD)
# W = 375182384, p-value = 0.4349




### data_setE

data_setE$DDRgene <- as.factor(data_setE$DDRgene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ DDRgene, data = data_setE)
# 



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = DDRgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ DDRgene, data = data_setE)
# W = 244580336, p-value = 0.2637



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = DDRgene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = DDRgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ DDRgene, data = data_setE)
# W = 241847028, p-value = 0.7421



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = DDRgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = DDRgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ DDRgene, data = data_setE)
# W = 243084012, p-value = 0.745


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(DDRgene), scales = "free_y") +
  theme_light()
# 










# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_DDRgenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_DDRgenes.R')
