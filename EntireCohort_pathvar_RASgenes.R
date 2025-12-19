# 20/06/2024

# EntireCohort_pathvar_RASgenes.R

# this script is for the RAS associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)


# RASgene
# PHEX
# FGFR1
# FGFR3
# RASA2
# RALA
# SHOC2
# SPRED2
# MRAS
# FGD1









# RASgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')


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

















# RASgene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_Rasgene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_Rasgene_, "pathvar_Rasgene_.csv", row.names=FALSE)
# system('dx upload pathvar_Rasgene_.csv --path Emily-folder/results_2/pathvar_Rasgene_.csv')


# RASgene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ RASgene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_Rasgene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_incidence_1, "pathvar_Rasgene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_incidence_1.csv --path Emily-folder/results_2/pathvar_Rasgene_incidence_1.csv')


# GOIs
LR_2 <- glm(formula = cancer_ROOC_dummy ~ PHEX+FGFR3+RASA2+SHOC2+SPRED2+FGD1, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_Rasgene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_incidence_2, "pathvar_Rasgene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_incidence_2.csv --path Emily-folder/results_2/pathvar_Rasgene_incidence_2.csv')







### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "RASgene", 'PHEX', 'FGFR3', 'RASA2', 
                                     'SHOC2', 'SPRED2', 'FGD1')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     13

with(data_setA, table(FGD1))


# RASgene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ RASgene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_Rasgene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_incidence_3, "pathvar_Rasgene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_incidence_3.csv --path Emily-folder/results_2/pathvar_Rasgene_incidence_3.csv')


# age + ethnicity + smoking + sex + RASgene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_Rasgene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_incidence_4, "pathvar_Rasgene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_incidence_4.csv --path Emily-folder/results_2/pathvar_Rasgene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ PHEX+FGFR3+RASA2+SHOC2+SPRED2+FGD1, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_Rasgene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_incidence_5, "pathvar_Rasgene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_incidence_5.csv --path Emily-folder/results_2/pathvar_Rasgene_incidence_5.csv')


# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+PHEX+FGFR3+RASA2+SHOC2+SPRED2+FGD1, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_Rasgene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_incidence_6, "pathvar_Rasgene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_incidence_6.csv --path Emily-folder/results_2/pathvar_Rasgene_incidence_6.csv')




### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "RASgene", 'PHEX', 'FGFR3', 'RASA2', 
                                     'SHOC2', 'SPRED2', 'FGD1')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     18

with(data_setB, table(FGD1))
# FGFR3 = 12
# SHOC2 = 4

# RASgene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ RASgene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_Rasgene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_incidence_7, "pathvar_Rasgene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_incidence_7.csv --path Emily-folder/results_2/pathvar_Rasgene_incidence_7.csv')


# all variables + RASgene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+RASgene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_Rasgene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_incidence_8, "pathvar_Rasgene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_incidence_8.csv --path Emily-folder/results_2/pathvar_Rasgene_incidence_8.csv')



# GOIs alone
# LR_9 <- glm(formula = cancer_ROOC_dummy ~ PHEX+FGFR3+RASA2+SHOC2+SPRED2+FGD1, 
#             data = data_setB, family = binomial)

## this was giving the following error message: 
# Error in profile.glm(object, which = parm, alpha = (1 - level)/4, trace = trace) : 
#   profiling has found a better solution, so original fit had not converged

# removing just SHOC2 still gives the same error message
# removing FGFR3 no longer gives the error message

# removed SHOC2 (only 4 participants)
LR_9 <- glm(formula = cancer_ROOC_dummy ~ PHEX+RASA2+SHOC2+SPRED2+FGD1, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_Rasgene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_incidence_9, "pathvar_Rasgene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_incidence_9.csv --path Emily-folder/results_2/pathvar_Rasgene_incidence_9.csv')



# all variables + GOIs
# LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+PHEX+FGFR3+RASA2+SHOC2+SPRED2+FGD1, 
#              data = data_setB, family = binomial)

## same error message

# tried with FGFR3 removed
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+PHEX+RASA2+SHOC2+SPRED2+FGD1, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_Rasgene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_incidence_10, "pathvar_Rasgene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_incidence_10.csv --path Emily-folder/results_2/pathvar_Rasgene_incidence_10.csv')










### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(RASgene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(RASgene, cancer_ROOC_dummy)))
# p-value = 0.3665

with(EntireCohort_pathvar, table(RASgene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(RASgene, Ethnicity_White)))
# p-value = 3.714e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5579118 0.7666985
# sample estimates:
#   odds ratio 
# 0.6522771

with(EntireCohort_pathvar, table(RASgene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(RASgene, Ever_smoked_dummy)))
# p-value = 0.8445

with(EntireCohort_pathvar, table(RASgene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(RASgene, Sex_Male)))
# p-value = 0.881

with(EntireCohort_pathvar, table(RASgene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(RASgene, University_education)))
# p-value = 0.2822

with(EntireCohort_pathvar, table(RASgene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(RASgene, Annual_income))	)
# p-value = 0.7708

with(EntireCohort_pathvar, table(RASgene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(RASgene, Smoker_previous)))
# p-value = 0.3458

with(EntireCohort_pathvar, table(RASgene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(RASgene, Smoker_current))	)
# p-value = 0.01005
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.039660 1.347956
# sample estimates:
#   odds ratio 
# 1.185705



### data_setA
with(data_setA, table(RASgene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(RASgene, cancer_ROOC_dummy)))
# p-value = 0.3508

with(data_setA, table(RASgene, Ethnicity_White))			
fisher.test(with(data_setA, table(RASgene, Ethnicity_White)))
# p-value = 3.319e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5560260 0.7648145
# sample estimates:
#   odds ratio 
# 0.6503597

with(data_setA, table(RASgene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(RASgene, Ever_smoked_dummy)))
# p-value = 0.9304

with(data_setA, table(RASgene, Sex_Male))			
fisher.test(with(data_setA, table(RASgene, Sex_Male)))
# p-value = 0.9829



### data_setB
with(data_setB, table(RASgene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(RASgene, cancer_ROOC_dummy)))
# p-value = 0.632

with(data_setB, table(RASgene, Ethnicity_White))			
fisher.test(with(data_setB, table(RASgene, Ethnicity_White)))
# p-value = 0.001417
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5902233 0.8824742
# sample estimates:
#   odds ratio 
# 0.71845

with(data_setB, table(RASgene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(RASgene, Ever_smoked_dummy)))
# p-value = 0.8267

with(data_setB, table(RASgene, Sex_Male))			
fisher.test(with(data_setB, table(RASgene, Sex_Male)))
# p-value = 0.8118

with(data_setB, table(RASgene, University_education))			
fisher.test(with(data_setB, table(RASgene, University_education)))
# p-value = 0.6373

with(data_setB, table(RASgene, Annual_income))			
fisher.test(with(data_setB, table(RASgene, Annual_income)))
# p-value = 0.8503






### can plot the continuous variables: 

# make RASgene into factor
df_for_plots <- EntireCohort_pathvar
df_for_plots$RASgene <- as.factor(df_for_plots$RASgene)

# try a violin plot

# Age_at_recruitment  
df_for_plots %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ RASgene, data = df_for_plots)
# W = 523032630, p-value = 0.3946

# BMI  
df_for_plots %>% 
  ggplot(mapping = aes(x = RASgene, y = BMI)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = RASgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ RASgene, data = df_for_plots)
# W = 506527972, p-value = 0.3631

# Standing_height 
df_for_plots %>% 
  ggplot(mapping = aes(x = RASgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = RASgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ RASgene, data = df_for_plots)
# W = 527780479, p-value = 0.02443

# make this a bar plot
# Weekly_exercise
df_for_plots %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(RASgene), scales = "free_y") +
  theme_light()
# 



### data_setA

# make RASgene into factor
data_setA$RASgene <- as.factor(data_setA$RASgene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ RASgene, data = data_setA)
# W = 512712511, p-value = 0.4282



### data_setB

# make RASgene into factor
data_setB$RASgene <- as.factor(data_setB$RASgene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ RASgene, data = data_setB)
# W = 344091784, p-value = 0.101

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = RASgene, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = RASgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ RASgene, data = data_setB)
# W = 333272926, p-value = 0.4832

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = RASgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = RASgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ RASgene, data = data_setB)
# W = 345632809, p-value = 0.0485


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(RASgene), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# RASgene Type ------------------------------------------------------------

# make sure the execute the code in the 'RASgenes' section

## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_Rasgene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_type_1, "pathvar_Rasgene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_type_1.csv --path Emily-folder/results_2/pathvar_Rasgene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_Rasgene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_type_2, "pathvar_Rasgene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_type_2.csv --path Emily-folder/results_2/pathvar_Rasgene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_Rasgene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_type_3, "pathvar_Rasgene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_type_3.csv --path Emily-folder/results_2/pathvar_Rasgene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_Rasgene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_type_4, "pathvar_Rasgene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_type_4.csv --path Emily-folder/results_2/pathvar_Rasgene_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_Rasgene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_type_5, "pathvar_Rasgene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_type_5.csv --path Emily-folder/results_2/pathvar_Rasgene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_Rasgene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_type_6, "pathvar_Rasgene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_type_6.csv --path Emily-folder/results_2/pathvar_Rasgene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_Rasgene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_type_7, "pathvar_Rasgene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_type_7.csv --path Emily-folder/results_2/pathvar_Rasgene_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_Rasgene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_type_8, "pathvar_Rasgene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_type_8.csv --path Emily-folder/results_2/pathvar_Rasgene_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_Rasgene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_type_9, "pathvar_Rasgene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_type_9.csv --path Emily-folder/results_2/pathvar_Rasgene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_Rasgene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_type_10, "pathvar_Rasgene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_type_10.csv --path Emily-folder/results_2/pathvar_Rasgene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_Rasgene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_type_11, "pathvar_Rasgene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_type_11.csv --path Emily-folder/results_2/pathvar_Rasgene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_Rasgene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_type_12, "pathvar_Rasgene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_type_12.csv --path Emily-folder/results_2/pathvar_Rasgene_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_Rasgene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_type_13, "pathvar_Rasgene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_type_13.csv --path Emily-folder/results_2/pathvar_Rasgene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_Rasgene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_type_14, "pathvar_Rasgene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_type_14.csv --path Emily-folder/results_2/pathvar_Rasgene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_Rasgene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Rasgene_type_15, "pathvar_Rasgene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_type_15.csv --path Emily-folder/results_2/pathvar_Rasgene_type_15.csv')





# RASgene Age -------------------------------------------------------------

# make sure the execute the code in the 'RASgenes' section

### linear regression
# "RASgene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ RASgene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_Rasgene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_Rasgene_age_1, "pathvar_Rasgene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_age_1.csv --path Emily-folder/results_2/pathvar_Rasgene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ PHEX+FGFR3+RASA2+SHOC2+SPRED2+FGD1, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_Rasgene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_Rasgene_age_2, "pathvar_Rasgene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_age_2.csv --path Emily-folder/results_2/pathvar_Rasgene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "RASgene", 'PHEX', 'FGFR3', 'RASA2', 
                                     'SHOC2', 'SPRED2', 'FGD1')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     13

with(data_setD, table(FGD1))
# FGFR3 = 1
# SHOC2 = 2


data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "RASgene", 'PHEX', 'FGFR3', 
                                     'RASA2', 'SHOC2', 'SPRED2', 'FGD1')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    18


with(data_setE, table(FGD1))
# FGFR3 = 0
# SHOC2 = 1



### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ RASgene, data = data_setD)
# summary(LM_15)
pathvar_Rasgene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_Rasgene_age_3, "pathvar_Rasgene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_age_3.csv --path Emily-folder/results_2/pathvar_Rasgene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+RASgene, data = data_setD)
# summary(LM_16)
pathvar_Rasgene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_Rasgene_age_4, "pathvar_Rasgene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_age_4.csv --path Emily-folder/results_2/pathvar_Rasgene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ PHEX+RASA2+SPRED2+FGD1, data = data_setD)
# summary(LM_17)
pathvar_Rasgene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_Rasgene_age_5, "pathvar_Rasgene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_age_5.csv --path Emily-folder/results_2/pathvar_Rasgene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+PHEX+RASA2+SPRED2+FGD1, data = data_setD)
# summary(LM_18)
pathvar_Rasgene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_Rasgene_age_6, "pathvar_Rasgene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_age_6.csv --path Emily-folder/results_2/pathvar_Rasgene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ RASgene, data = data_setE)
# summary(LM_15)
pathvar_Rasgene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_Rasgene_age_7, "pathvar_Rasgene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_age_7.csv --path Emily-folder/results_2/pathvar_Rasgene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+RASgene, data = data_setE)
# summary(LM_16)
pathvar_Rasgene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_Rasgene_age_8, "pathvar_Rasgene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_age_8.csv --path Emily-folder/results_2/pathvar_Rasgene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ PHEX+RASA2+SPRED2+FGD1, data = data_setE)
# summary(LM_17)
pathvar_Rasgene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_Rasgene_age_9, "pathvar_Rasgene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_age_9.csv --path Emily-folder/results_2/pathvar_Rasgene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+PHEX+RASA2+SPRED2+FGD1, data = data_setE)
# summary(LM_18)
pathvar_Rasgene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_Rasgene_age_10, "pathvar_Rasgene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_Rasgene_age_10.csv --path Emily-folder/results_2/pathvar_Rasgene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(df_cancer_diagnosis, table(RASgene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(RASgene, Ethnicity_White)))
# p-value = 0.0719

with(df_cancer_diagnosis, table(RASgene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(RASgene, Ever_smoked_dummy)))
# p-value = 0.8528

with(df_cancer_diagnosis, table(RASgene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(RASgene, Sex_Male)))
# p-value = 0.244

with(df_cancer_diagnosis, table(RASgene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(RASgene, University_education)))
# p-value = 0.883

with(df_cancer_diagnosis, table(RASgene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(RASgene, Annual_income)))
# p-value = 0.5475


### data_setD
with(data_setD, table(RASgene, Ethnicity_White))			
fisher.test(with(data_setD, table(RASgene, Ethnicity_White)))
# p-value = 0.06952

with(data_setD, table(RASgene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(RASgene, Ever_smoked_dummy)))
# p-value = 0.6751

with(data_setD, table(RASgene, Sex_Male))			
fisher.test(with(data_setD, table(RASgene, Sex_Male)))
# p-value = 0.2411


### data_setE
with(data_setE, table(RASgene, Ethnicity_White))			
fisher.test(with(data_setE, table(RASgene, Ethnicity_White)))
# p-value = 0.1751

with(data_setE, table(RASgene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(RASgene, Ever_smoked_dummy)))
# p-value = 0.5335

with(data_setE, table(RASgene, Sex_Male))			
fisher.test(with(data_setE, table(RASgene, Sex_Male)))
# p-value = 0.4539

with(data_setE, table(RASgene, University_education))			
fisher.test(with(data_setE, table(RASgene, University_education)))
# p-value = 0.7911

with(data_setE, table(RASgene, Annual_income))			
fisher.test(with(data_setE, table(RASgene, Annual_income)))
# p-value = 0.3946



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 

# make RASgene into factor
df_cancer_diagnosis$RASgene <- as.factor(df_cancer_diagnosis$RASgene)

# try a violin plot



# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ RASgene, data = df_cancer_diagnosis)
#



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ RASgene, data = df_cancer_diagnosis)
# W = 28048511, p-value = 0.7419



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = RASgene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = RASgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ RASgene, data = df_cancer_diagnosis)
# W = 27743922, p-value = 0.9008



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = RASgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = RASgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ RASgene, data = df_cancer_diagnosis)
# W = 27790813, p-value = 0.8846



# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(RASgene), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$RASgene <- as.factor(data_setD$RASgene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ RASgene, data = data_setD)
# 



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ RASgene, data = data_setD)
# W = 27343225, p-value = 0.895




### data_setE

data_setE$RASgene <- as.factor(data_setE$RASgene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ RASgene, data = data_setE)
# 



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = RASgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ RASgene, data = data_setE)
# W = 18202884, p-value = 0.3707



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = RASgene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = RASgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ RASgene, data = data_setE)
# W = 17431586, p-value = 0.5394



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = RASgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = RASgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ RASgene, data = data_setE)
# W = 17964982, p-value = 0.6677


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(RASgene), scales = "free_y") +
  theme_light()
# 










# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_RASgenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_RASgenes.R')
