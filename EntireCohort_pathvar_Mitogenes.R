# 09/07/2024

# EntireCohort_pathvar_Mitogenes.R

# this script is for the Mitochondria associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)


# Mitogene
# PYCR1
# PISD
# SLC25A19
# BCS1L
# ATAD3A
# PDHA1
# PAM16
# GPX4
# AIFM1







# Mitogenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')


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

















# Mitogene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_Mitogene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_Mitogene_, "pathvar_Mitogene_.csv", row.names=FALSE)
# system('dx upload pathvar_Mitogene_.csv --path Emily-folder/results_2/pathvar_Mitogene_.csv')


# Mitogene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ Mitogene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_Mitogene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_incidence_1, "pathvar_Mitogene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_incidence_1.csv --path Emily-folder/results_2/pathvar_Mitogene_incidence_1.csv')


# GOIs
LR_2 <- glm(formula = cancer_ROOC_dummy ~ PYCR1+PISD+SLC25A19+BCS1L+ATAD3A+PDHA1+PAM16+GPX4, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_Mitogene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_incidence_2, "pathvar_Mitogene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_incidence_2.csv --path Emily-folder/results_2/pathvar_Mitogene_incidence_2.csv')







### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "Mitogene", 'PYCR1', 'PISD', 'SLC25A19', 'BCS1L', 
                                     'ATAD3A', 'PDHA1', 'PAM16', 'GPX4')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     15

with(data_setA, table(GPX4))



# Mitogene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ Mitogene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_Mitogene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_incidence_3, "pathvar_Mitogene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_incidence_3.csv --path Emily-folder/results_2/pathvar_Mitogene_incidence_3.csv')


# age + ethnicity + smoking + sex + Mitogene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_Mitogene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_incidence_4, "pathvar_Mitogene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_incidence_4.csv --path Emily-folder/results_2/pathvar_Mitogene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ PYCR1+PISD+SLC25A19+BCS1L+ATAD3A+PDHA1+PAM16+GPX4, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_Mitogene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_incidence_5, "pathvar_Mitogene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_incidence_5.csv --path Emily-folder/results_2/pathvar_Mitogene_incidence_5.csv')


# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+PYCR1+PISD+SLC25A19+BCS1L+ATAD3A+PDHA1+PAM16+GPX4, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_Mitogene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_incidence_6, "pathvar_Mitogene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_incidence_6.csv --path Emily-folder/results_2/pathvar_Mitogene_incidence_6.csv')




### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "Mitogene", 'PYCR1', 'PISD', 'SLC25A19', 'BCS1L', 
                                     'ATAD3A', 'PDHA1', 'PAM16', 'GPX4')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     20

with(data_setB, table(GPX4))


# Mitogene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ Mitogene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_Mitogene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_incidence_7, "pathvar_Mitogene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_incidence_7.csv --path Emily-folder/results_2/pathvar_Mitogene_incidence_7.csv')


# all variables + Mitogene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+Mitogene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_Mitogene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_incidence_8, "pathvar_Mitogene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_incidence_8.csv --path Emily-folder/results_2/pathvar_Mitogene_incidence_8.csv')


# GOIs alone
LR_9 <- glm(formula = cancer_ROOC_dummy ~ PYCR1+PISD+SLC25A19+BCS1L+ATAD3A+PDHA1+PAM16+GPX4, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_Mitogene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_incidence_9, "pathvar_Mitogene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_incidence_9.csv --path Emily-folder/results_2/pathvar_Mitogene_incidence_9.csv')


# all variables + GOIs
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+PYCR1+PISD+SLC25A19+BCS1L+ATAD3A+PDHA1+PAM16+GPX4, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_Mitogene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_incidence_10, "pathvar_Mitogene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_incidence_10.csv --path Emily-folder/results_2/pathvar_Mitogene_incidence_10.csv')










### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(Mitogene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(Mitogene, cancer_ROOC_dummy)))
# p-value = 0.3829

with(EntireCohort_pathvar, table(Mitogene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(Mitogene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5704627 0.6874810
# sample estimates:
#   odds ratio 
# 0.6257258

with(EntireCohort_pathvar, table(Mitogene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(Mitogene, Ever_smoked_dummy)))
# p-value = 0.01009
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.8885641 0.9844278
# sample estimates:
#   odds ratio 
# 0.9352327

with(EntireCohort_pathvar, table(Mitogene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(Mitogene, Sex_Male)))
# p-value = 0.1178

with(EntireCohort_pathvar, table(Mitogene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(Mitogene, University_education)))
# p-value = 0.02554
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.007326 1.121197
# sample estimates:
#   odds ratio 
# 1.062852

with(EntireCohort_pathvar, table(Mitogene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(Mitogene, Annual_income))	)
# p-value = 0.4084

with(EntireCohort_pathvar, table(Mitogene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(Mitogene, Smoker_previous)))
# p-value = 0.154

with(EntireCohort_pathvar, table(Mitogene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(Mitogene, Smoker_current))	)
# p-value = 0.252



### data_setA
with(data_setA, table(Mitogene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(Mitogene, cancer_ROOC_dummy)))
# p-value = 0.4589

with(data_setA, table(Mitogene, Ethnicity_White))			
fisher.test(with(data_setA, table(Mitogene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5715314 0.6893943
# sample estimates:
#   odds ratio 
# 0.6271575

with(data_setA, table(Mitogene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(Mitogene, Ever_smoked_dummy)))
# p-value = 0.01108
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.8890307 0.9852098
# sample estimates:
#   odds ratio 
# 0.9358501

with(data_setA, table(Mitogene, Sex_Male))			
fisher.test(with(data_setA, table(Mitogene, Sex_Male)))
# p-value = 0.1317



### data_setB
with(data_setB, table(Mitogene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(Mitogene, cancer_ROOC_dummy)))
# p-value = 0.601

with(data_setB, table(Mitogene, Ethnicity_White))			
fisher.test(with(data_setB, table(Mitogene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5446081 0.6804574
# sample estimates:
#   odds ratio 
# 0.6079813

with(data_setB, table(Mitogene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(Mitogene, Ever_smoked_dummy)))
# p-value = 0.262

with(data_setB, table(Mitogene, Sex_Male))			
fisher.test(with(data_setB, table(Mitogene, Sex_Male)))
# p-value = 0.2915

with(data_setB, table(Mitogene, University_education))			
fisher.test(with(data_setB, table(Mitogene, University_education)))
# p-value = 0.06124

with(data_setB, table(Mitogene, Annual_income))			
fisher.test(with(data_setB, table(Mitogene, Annual_income)))
# p-value = 0.3576






### can plot the continuous variables: 

# make Mitogene into factor
df_for_plots <- EntireCohort_pathvar
df_for_plots$Mitogene <- as.factor(df_for_plots$Mitogene)

# try a violin plot

# Age_at_recruitment  
df_for_plots %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ Mitogene, data = df_for_plots)
# W = 1460677948, p-value = 0.04305

# BMI  
df_for_plots %>% 
  ggplot(mapping = aes(x = Mitogene, y = BMI)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = Mitogene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ Mitogene, data = df_for_plots)
# W = 1430530245, p-value = 0.8011

# Standing_height 
df_for_plots %>% 
  ggplot(mapping = aes(x = Mitogene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = Mitogene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ Mitogene, data = df_for_plots)
# W = 1430586516, p-value = 0.9988

# make this a bar plot
# Weekly_exercise
df_for_plots %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(Mitogene), scales = "free_y") +
  theme_light()
# 



### data_setA

# make Mitogene into factor
data_setA$Mitogene <- as.factor(data_setA$Mitogene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ Mitogene, data = data_setA)
# W = 1433278182, p-value = 0.05818



### data_setB

# make Mitogene into factor
data_setB$Mitogene <- as.factor(data_setB$Mitogene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ Mitogene, data = data_setB)
# W = 954535284, p-value = 0.04819

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = Mitogene, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = Mitogene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ Mitogene, data = data_setB)
# W = 941059274, p-value = 0.8182

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = Mitogene, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = Mitogene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ Mitogene, data = data_setB)
# W = 939182076, p-value = 0.9895


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(Mitogene), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# Mitogene Type ------------------------------------------------------------

# make sure the execute the code in the 'Mitogenes' section

## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_Mitogene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_type_1, "pathvar_Mitogene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_type_1.csv --path Emily-folder/results_2/pathvar_Mitogene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_Mitogene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_type_2, "pathvar_Mitogene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_type_2.csv --path Emily-folder/results_2/pathvar_Mitogene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_Mitogene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_type_3, "pathvar_Mitogene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_type_3.csv --path Emily-folder/results_2/pathvar_Mitogene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_Mitogene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_type_4, "pathvar_Mitogene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_type_4.csv --path Emily-folder/results_2/pathvar_Mitogene_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_Mitogene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_type_5, "pathvar_Mitogene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_type_5.csv --path Emily-folder/results_2/pathvar_Mitogene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_Mitogene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_type_6, "pathvar_Mitogene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_type_6.csv --path Emily-folder/results_2/pathvar_Mitogene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_Mitogene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_type_7, "pathvar_Mitogene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_type_7.csv --path Emily-folder/results_2/pathvar_Mitogene_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_Mitogene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_type_8, "pathvar_Mitogene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_type_8.csv --path Emily-folder/results_2/pathvar_Mitogene_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_Mitogene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_type_9, "pathvar_Mitogene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_type_9.csv --path Emily-folder/results_2/pathvar_Mitogene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_Mitogene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_type_10, "pathvar_Mitogene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_type_10.csv --path Emily-folder/results_2/pathvar_Mitogene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_Mitogene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_type_11, "pathvar_Mitogene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_type_11.csv --path Emily-folder/results_2/pathvar_Mitogene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_Mitogene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_type_12, "pathvar_Mitogene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_type_12.csv --path Emily-folder/results_2/pathvar_Mitogene_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_Mitogene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_type_13, "pathvar_Mitogene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_type_13.csv --path Emily-folder/results_2/pathvar_Mitogene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_Mitogene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_type_14, "pathvar_Mitogene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_type_14.csv --path Emily-folder/results_2/pathvar_Mitogene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_Mitogene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Mitogene_type_15, "pathvar_Mitogene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_type_15.csv --path Emily-folder/results_2/pathvar_Mitogene_type_15.csv')






# Mitogene Age -------------------------------------------------------------

# make sure the execute the code in the 'Mitogenes' section

### linear regression
# "Mitogene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ Mitogene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_Mitogene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_Mitogene_age_1, "pathvar_Mitogene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_age_1.csv --path Emily-folder/results_2/pathvar_Mitogene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ PYCR1+PISD+SLC25A19+BCS1L+ATAD3A+PDHA1+PAM16+GPX4, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_Mitogene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_Mitogene_age_2, "pathvar_Mitogene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_age_2.csv --path Emily-folder/results_2/pathvar_Mitogene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "Mitogene", 'PYCR1', 'PISD', 'SLC25A19', 'BCS1L', 'ATAD3A', 
                                     'PDHA1', 'PAM16', 'GPX4')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     15

with(data_setD, table(GPX4))
# SLC25A19 = 67
# PDHA1 = 54
# PAM16 = 61
# GPX4 = 54



data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "Mitogene", 'PYCR1', 'PISD', 'SLC25A19', 
                                     'BCS1L', 'ATAD3A', 'PDHA1', 'PAM16', 'GPX4')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    20

with(data_setE, table(GPX4))
# PYCR1 = 83
# SLC25A19 = 56
# PDHA1 = 46
# PAM16 = 48
# GPX4 = 46




### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ Mitogene, data = data_setD)
# summary(LM_15)
pathvar_Mitogene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_Mitogene_age_3, "pathvar_Mitogene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_age_3.csv --path Emily-folder/results_2/pathvar_Mitogene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Mitogene, data = data_setD)
# summary(LM_16)
pathvar_Mitogene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_Mitogene_age_4, "pathvar_Mitogene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_age_4.csv --path Emily-folder/results_2/pathvar_Mitogene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ PYCR1+PISD+SLC25A19+BCS1L+ATAD3A+PDHA1+PAM16+GPX4, data = data_setD)
# summary(LM_17)
pathvar_Mitogene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_Mitogene_age_5, "pathvar_Mitogene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_age_5.csv --path Emily-folder/results_2/pathvar_Mitogene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+PYCR1+PISD+SLC25A19+BCS1L+ATAD3A+PDHA1+PAM16+GPX4, data = data_setD)
# summary(LM_18)
pathvar_Mitogene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_Mitogene_age_6, "pathvar_Mitogene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_age_6.csv --path Emily-folder/results_2/pathvar_Mitogene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ Mitogene, data = data_setE)
# summary(LM_15)
pathvar_Mitogene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_Mitogene_age_7, "pathvar_Mitogene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_age_7.csv --path Emily-folder/results_2/pathvar_Mitogene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+Mitogene, data = data_setE)
# summary(LM_16)
pathvar_Mitogene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_Mitogene_age_8, "pathvar_Mitogene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_age_8.csv --path Emily-folder/results_2/pathvar_Mitogene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ PYCR1+PISD+SLC25A19+BCS1L+ATAD3A+PDHA1+PAM16+GPX4, data = data_setE)
# summary(LM_17)
pathvar_Mitogene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_Mitogene_age_9, "pathvar_Mitogene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_age_9.csv --path Emily-folder/results_2/pathvar_Mitogene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+PYCR1+PISD+SLC25A19+BCS1L+ATAD3A+PDHA1+PAM16+GPX4, data = data_setE)
# summary(LM_18)
pathvar_Mitogene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_Mitogene_age_10, "pathvar_Mitogene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_Mitogene_age_10.csv --path Emily-folder/results_2/pathvar_Mitogene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(df_cancer_diagnosis, table(Mitogene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(Mitogene, Ethnicity_White)))
# p-value = 3.326e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.4403135 0.7375975
# sample estimates:
#   odds ratio 
# 0.5656235

with(df_cancer_diagnosis, table(Mitogene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(Mitogene, Ever_smoked_dummy)))
# p-value = 0.7827

with(df_cancer_diagnosis, table(Mitogene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(Mitogene, Sex_Male)))
# p-value = 0.3373

with(df_cancer_diagnosis, table(Mitogene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(Mitogene, University_education)))
# p-value = 0.9535

with(df_cancer_diagnosis, table(Mitogene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(Mitogene, Annual_income)))
# p-value = 0.886


### data_setD
with(data_setD, table(Mitogene, Ethnicity_White))			
fisher.test(with(data_setD, table(Mitogene, Ethnicity_White)))
# p-value = 3.098e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.4370573 0.7321687
# sample estimates:
#   odds ratio 
# 0.56143

with(data_setD, table(Mitogene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(Mitogene, Ever_smoked_dummy)))
# p-value = 0.7823

with(data_setD, table(Mitogene, Sex_Male))			
fisher.test(with(data_setD, table(Mitogene, Sex_Male)))
# p-value = 0.407


### data_setE
with(data_setE, table(Mitogene, Ethnicity_White))			
fisher.test(with(data_setE, table(Mitogene, Ethnicity_White)))
# p-value = 1.133e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.3709267 0.6729782
# sample estimates:
#   odds ratio 
# 0.494563

with(data_setE, table(Mitogene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(Mitogene, Ever_smoked_dummy)))
# p-value = 0.7104

with(data_setE, table(Mitogene, Sex_Male))			
fisher.test(with(data_setE, table(Mitogene, Sex_Male)))
# p-value = 0.2834

with(data_setE, table(Mitogene, University_education))			
fisher.test(with(data_setE, table(Mitogene, University_education)))
# p-value = 1

with(data_setE, table(Mitogene, Annual_income))			
fisher.test(with(data_setE, table(Mitogene, Annual_income)))
# p-value = 0.7441



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 

# make Mitogene into factor
df_cancer_diagnosis$Mitogene <- as.factor(df_cancer_diagnosis$Mitogene)

# try a violin plot



# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ Mitogene, data = df_cancer_diagnosis)
#



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ Mitogene, data = df_cancer_diagnosis)
# W = 79810476, p-value = 0.3047



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Mitogene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Mitogene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ Mitogene, data = df_cancer_diagnosis)
# W = 77530367, p-value = 0.6875



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Mitogene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Mitogene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ Mitogene, data = df_cancer_diagnosis)
# W = 80621225, p-value = 0.04004



# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(Mitogene), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$Mitogene <- as.factor(data_setD$Mitogene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ Mitogene, data = data_setD)
# 



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ Mitogene, data = data_setD)
# W = 78536204, p-value = 0.3258




### data_setE

data_setE$Mitogene <- as.factor(data_setE$Mitogene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ Mitogene, data = data_setE)
# 



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = Mitogene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ Mitogene, data = data_setE)
# W = 50351614, p-value = 0.6042



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = Mitogene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = Mitogene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ Mitogene, data = data_setE)
# W = 49043836, p-value = 0.3141



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = Mitogene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = Mitogene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ Mitogene, data = data_setE)
# W = 51636687, p-value = 0.04372


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(Mitogene), scales = "free_y") +
  theme_light()
# 












# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_Mitogenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_Mitogenes.R')
