# 20/06/2024

# EntireCohort_pathvar_THgenes.R

# this script is for the Thyroid and parathyroid hormone associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)


# Thyroid Hormone - associated genes
# PRKAR1A
# SIK3
# TBCE
# TRIP11
# PTHLH
# PTH1R
# THRA
# THRB
# GNAS





# THgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')


# go thr all the genes to check that they're present
grep("GNAS", EntireCohort_pathvar$GENE)
# all seem to be present except for these genes
# PRKAR1A = 0
# THRA = 0



# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$THgene <- NA
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE1 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS', 'THRB')] <- 1
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE2 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS', 'THRB')] <- 1
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE3 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS', 'THRB')] <- 1
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE4 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS', 'THRB')] <- 1
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE5 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS', 'THRB')] <- 1
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE6 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS', 'THRB')] <- 1
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE7 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS', 'THRB')] <- 1
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE8 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS', 'THRB')] <- 1
EntireCohort_pathvar$THgene[EntireCohort_pathvar$GENE9 %in% c('SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 'GNAS', 'THRB')] <- 1


EntireCohort_pathvar$THgene[is.na(EntireCohort_pathvar$THgene)] <- 0
EntireCohort_pathvar$THgene <- as.numeric(EntireCohort_pathvar$THgene)
unique(EntireCohort_pathvar$THgene)
# 0 1 

table(EntireCohort_pathvar$THgene)
#      0      1 
# 462283   7516



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

















# THgene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_THgene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_THgene_, "pathvar_THgene_.csv", row.names=FALSE)
# system('dx upload pathvar_THgene_.csv --path Emily-folder/results_2/pathvar_THgene_.csv')


# THgene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ THgene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_THgene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_incidence_1, "pathvar_THgene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_THgene_incidence_1.csv --path Emily-folder/results_2/pathvar_THgene_incidence_1.csv')


# GOIs
LR_2 <- glm(formula = cancer_ROOC_dummy ~ SIK3+TBCE+TRIP11+PTHLH+PTH1R+GNAS+THRB, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_THgene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_incidence_2, "pathvar_THgene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_THgene_incidence_2.csv --path Emily-folder/results_2/pathvar_THgene_incidence_2.csv')






### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "THgene", 'SIK3', 'TBCE', 'TRIP11', 'PTHLH', 
                                     'PTH1R', 'GNAS', 'THRB')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     13

with(data_setA, table(THgene))
# THgene
# 0      1 
# 458426   7455

# THgene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ THgene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_THgene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_incidence_3, "pathvar_THgene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_THgene_incidence_3.csv --path Emily-folder/results_2/pathvar_THgene_incidence_3.csv')


# age + ethnicity + smoking + sex + THgene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_THgene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_incidence_4, "pathvar_THgene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_THgene_incidence_4.csv --path Emily-folder/results_2/pathvar_THgene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ SIK3+TBCE+TRIP11+PTHLH+PTH1R+GNAS+THRB, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_THgene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_incidence_5, "pathvar_THgene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_THgene_incidence_5.csv --path Emily-folder/results_2/pathvar_THgene_incidence_5.csv')


# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SIK3+TBCE+TRIP11+PTHLH+PTH1R+GNAS+THRB, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_THgene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_incidence_6, "pathvar_THgene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_THgene_incidence_6.csv --path Emily-folder/results_2/pathvar_THgene_incidence_6.csv')




### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "THgene", 'SIK3', 'TBCE', 'TRIP11', 'PTHLH', 
                                     'PTH1R', 'GNAS', 'THRB')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     18

with(data_setB, table(THgene))
# THgene
# 0      1 
# 375085   6074 

# THgene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ THgene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_THgene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_incidence_7, "pathvar_THgene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_THgene_incidence_7.csv --path Emily-folder/results_2/pathvar_THgene_incidence_7.csv')


# all variables + THgene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+THgene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_THgene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_incidence_8, "pathvar_THgene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_THgene_incidence_8.csv --path Emily-folder/results_2/pathvar_THgene_incidence_8.csv')


# GOIs alone
LR_9 <- glm(formula = cancer_ROOC_dummy ~ SIK3+TBCE+TRIP11+PTHLH+PTH1R+GNAS+THRB, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_THgene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_incidence_9, "pathvar_THgene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_THgene_incidence_9.csv --path Emily-folder/results_2/pathvar_THgene_incidence_9.csv')


# all variables + GOIs
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+SIK3+TBCE+TRIP11+PTHLH+PTH1R+GNAS+THRB, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_THgene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_incidence_10, "pathvar_THgene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_THgene_incidence_10.csv --path Emily-folder/results_2/pathvar_THgene_incidence_10.csv')










### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(THgene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(THgene, cancer_ROOC_dummy)))
# p-value = 0.9345

with(EntireCohort_pathvar, table(THgene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(THgene, Ethnicity_White)))
# p-value = 0.8759

with(EntireCohort_pathvar, table(THgene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(THgene, Ever_smoked_dummy)))
# p-value = 0.2638

with(EntireCohort_pathvar, table(THgene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(THgene, Sex_Male)))
# p-value = 0.7351

with(EntireCohort_pathvar, table(THgene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(THgene, University_education)))
# p-value = 0.5094

with(EntireCohort_pathvar, table(THgene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(THgene, Annual_income))	)
# p-value = 0.5256

with(EntireCohort_pathvar, table(THgene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(THgene, Smoker_previous)))
# p-value = 0.8639

with(EntireCohort_pathvar, table(THgene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(THgene, Smoker_current))	)
# p-value = 0.8645



### data_setA
with(data_setA, table(THgene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(THgene, cancer_ROOC_dummy)))
# p-value = 0.9781

with(data_setA, table(THgene, Ethnicity_White))			
fisher.test(with(data_setA, table(THgene, Ethnicity_White)))
# p-value = 0.7938

with(data_setA, table(THgene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(THgene, Ever_smoked_dummy)))
# p-value = 0.2947

with(data_setA, table(THgene, Sex_Male))			
fisher.test(with(data_setA, table(THgene, Sex_Male)))
# p-value = 0.6144



### data_setB
with(data_setB, table(THgene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(THgene, cancer_ROOC_dummy)))
# p-value = 0.8904

with(data_setB, table(THgene, Ethnicity_White))			
fisher.test(with(data_setB, table(THgene, Ethnicity_White)))
# p-value = 0.8012

with(data_setB, table(THgene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(THgene, Ever_smoked_dummy)))
# p-value = 0.1316

with(data_setB, table(THgene, Sex_Male))			
fisher.test(with(data_setB, table(THgene, Sex_Male)))
# p-value = 0.6882

with(data_setB, table(THgene, University_education))			
fisher.test(with(data_setB, table(THgene, University_education)))
# p-value = 0.8927

with(data_setB, table(THgene, Annual_income))			
fisher.test(with(data_setB, table(THgene, Annual_income)))
# p-value = 0.608






### can plot the continuous variables: 


# make THgene into factor
df_for_plot <- EntireCohort_pathvar
df_for_plot$THgene <- as.factor(df_for_plot$THgene)

# try a violin plot

# Age_at_recruitment  
df_for_plot %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ THgene, data = df_for_plot)
# W = 1732640874, p-value = 0.6919

# BMI  
df_for_plot %>% 
  ggplot(mapping = aes(x = THgene, y = BMI)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = THgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ THgene, data = df_for_plot)
# W = 1728620824, p-value = 0.5674

# Standing_height 
df_for_plot %>% 
  ggplot(mapping = aes(x = THgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_for_plot %>% 
  ggplot(mapping = aes(x = THgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ THgene, data = df_for_plot)
# W = 1722234724, p-value = 0.7337

# make this a bar plot
# Weekly_exercise
df_for_plot %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(THgene), scales = "free_y") +
  theme_light()
# 



### data_setA

# make THgene into factor
data_setA$THgene <- as.factor(data_setA$THgene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ THgene, data = data_setA)
# W = 1704530164, p-value = 0.7118



### data_setB

# make THgene into factor
data_setB$THgene <- as.factor(data_setB$THgene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ THgene, data = data_setB)
# W = 1142180449, p-value = 0.72

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = THgene, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = THgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ THgene, data = data_setB)
# W = 1139468975, p-value = 0.9685

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = THgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = THgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ THgene, data = data_setB)
# W = 1131534880, p-value = 0.3715


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(THgene), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# THgene Type ------------------------------------------------------------

# make sure the execute the code in the 'THgenes' section

## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_THgene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_type_1, "pathvar_THgene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_THgene_type_1.csv --path Emily-folder/results_2/pathvar_THgene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_THgene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_type_2, "pathvar_THgene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_THgene_type_2.csv --path Emily-folder/results_2/pathvar_THgene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_THgene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_type_3, "pathvar_THgene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_THgene_type_3.csv --path Emily-folder/results_2/pathvar_THgene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_THgene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_type_4, "pathvar_THgene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_THgene_type_4.csv --path Emily-folder/results_2/pathvar_THgene_type_4.csv')
# some warnings - check this


LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_THgene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_type_5, "pathvar_THgene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_THgene_type_5.csv --path Emily-folder/results_2/pathvar_THgene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_THgene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_type_6, "pathvar_THgene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_THgene_type_6.csv --path Emily-folder/results_2/pathvar_THgene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_THgene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_type_7, "pathvar_THgene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_THgene_type_7.csv --path Emily-folder/results_2/pathvar_THgene_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_THgene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_type_8, "pathvar_THgene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_THgene_type_8.csv --path Emily-folder/results_2/pathvar_THgene_type_8.csv')
# some warnings - check this


LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_THgene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_type_9, "pathvar_THgene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_THgene_type_9.csv --path Emily-folder/results_2/pathvar_THgene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_THgene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_type_10, "pathvar_THgene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_THgene_type_10.csv --path Emily-folder/results_2/pathvar_THgene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_THgene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_type_11, "pathvar_THgene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_THgene_type_11.csv --path Emily-folder/results_2/pathvar_THgene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_THgene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_type_12, "pathvar_THgene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_THgene_type_12.csv --path Emily-folder/results_2/pathvar_THgene_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_THgene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_type_13, "pathvar_THgene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_THgene_type_13.csv --path Emily-folder/results_2/pathvar_THgene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_THgene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_type_14, "pathvar_THgene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_THgene_type_14.csv --path Emily-folder/results_2/pathvar_THgene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_THgene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_THgene_type_15, "pathvar_THgene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_THgene_type_15.csv --path Emily-folder/results_2/pathvar_THgene_type_15.csv')













# THgene Age -------------------------------------------------------------

# make sure the execute the code in the 'THgenes' section

### linear regression
# "THgene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ THgene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_THgene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_THgene_age_1, "pathvar_THgene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_THgene_age_1.csv --path Emily-folder/results_2/pathvar_THgene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ SIK3+TBCE+TRIP11+PTHLH+PTH1R+GNAS+THRB, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_THgene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_THgene_age_2, "pathvar_THgene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_THgene_age_2.csv --path Emily-folder/results_2/pathvar_THgene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "THgene", 'SIK3', 'TBCE', 'TRIP11', 'PTHLH', 'PTH1R', 
                                     'GNAS', 'THRB')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     13

# check that there are enough participants in each gene category still
with(data_setD, table(THgene))
# THgene
# 0      1 
# 112555   1827


data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "THgene", 'SIK3', 'TBCE', 'TRIP11', 
                                     'PTHLH', 'PTH1R', 'GNAS', 'THRB')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    18

# check that there are enough participants in each gene category still
with(data_setE, table(THgene))
# THgene
# 0     1 
# 90649  1464




### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ THgene, data = data_setD)
# summary(LM_15)
pathvar_THgene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_THgene_age_3, "pathvar_THgene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_THgene_age_3.csv --path Emily-folder/results_2/pathvar_THgene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+THgene, data = data_setD)
# summary(LM_16)
pathvar_THgene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_THgene_age_4, "pathvar_THgene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_THgene_age_4.csv --path Emily-folder/results_2/pathvar_THgene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ SIK3+TBCE+TRIP11+PTHLH+PTH1R+GNAS+THRB, data = data_setD)
# summary(LM_17)
pathvar_THgene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_THgene_age_5, "pathvar_THgene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_THgene_age_5.csv --path Emily-folder/results_2/pathvar_THgene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+SIK3+TBCE+TRIP11+PTHLH+PTH1R+GNAS+THRB, data = data_setD)
# summary(LM_18)
pathvar_THgene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_THgene_age_6, "pathvar_THgene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_THgene_age_6.csv --path Emily-folder/results_2/pathvar_THgene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ THgene, data = data_setE)
# summary(LM_15)
pathvar_THgene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_THgene_age_7, "pathvar_THgene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_THgene_age_7.csv --path Emily-folder/results_2/pathvar_THgene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+THgene, data = data_setE)
# summary(LM_16)
pathvar_THgene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_THgene_age_8, "pathvar_THgene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_THgene_age_8.csv --path Emily-folder/results_2/pathvar_THgene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ SIK3+TBCE+TRIP11+PTHLH+PTH1R+GNAS+THRB, data = data_setE)
# summary(LM_17)
pathvar_THgene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_THgene_age_9, "pathvar_THgene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_THgene_age_9.csv --path Emily-folder/results_2/pathvar_THgene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+SIK3+TBCE+TRIP11+PTHLH+PTH1R+GNAS+THRB, data = data_setE)
# summary(LM_18)
pathvar_THgene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_THgene_age_10, "pathvar_THgene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_THgene_age_10.csv --path Emily-folder/results_2/pathvar_THgene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(df_cancer_diagnosis, table(THgene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(THgene, Ethnicity_White)))
# p-value = 0.2751

with(df_cancer_diagnosis, table(THgene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(THgene, Ever_smoked_dummy)))
# p-value = 0.3408

with(df_cancer_diagnosis, table(THgene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(THgene, Sex_Male)))
# p-value = 0.3343

with(df_cancer_diagnosis, table(THgene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(THgene, University_education)))
# p-value = 0.8774

with(df_cancer_diagnosis, table(THgene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(THgene, Annual_income)))
# p-value = 0.5914


### data_setD
with(data_setD, table(THgene, Ethnicity_White))			
fisher.test(with(data_setD, table(THgene, Ethnicity_White)))
# p-value = 0.3056

with(data_setD, table(THgene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(THgene, Ever_smoked_dummy)))
# p-value = 0.3527

with(data_setD, table(THgene, Sex_Male))			
fisher.test(with(data_setD, table(THgene, Sex_Male)))
# p-value = 0.2464


### data_setE
with(data_setE, table(THgene, Ethnicity_White))			
fisher.test(with(data_setE, table(THgene, Ethnicity_White)))
# p-value = 0.1079

with(data_setE, table(THgene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(THgene, Ever_smoked_dummy)))
# p-value = 0.3661

with(data_setE, table(THgene, Sex_Male))			
fisher.test(with(data_setE, table(THgene, Sex_Male)))
# p-value = 0.1965

with(data_setE, table(THgene, University_education))			
fisher.test(with(data_setE, table(THgene, University_education)))
# p-value = 0.9777

with(data_setE, table(THgene, Annual_income))			
fisher.test(with(data_setE, table(THgene, Annual_income)))
# p-value = 0.5642



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 


# make THgene into factor
df_cancer_diagnosis$THgene <- as.factor(df_cancer_diagnosis$THgene)

# try a violin plot



# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ THgene, data = df_cancer_diagnosis)
# W = 104413663, p-value = 0.8607



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ THgene, data = df_cancer_diagnosis)
# W = 106339128, p-value = 0.2367



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = THgene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = THgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ THgene, data = df_cancer_diagnosis)
# W = 103927333, p-value = 0.961



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = THgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = THgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ THgene, data = df_cancer_diagnosis)
# W = 102382031, p-value = 0.1912



# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(THgene), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$THgene <- as.factor(data_setD$THgene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ THgene, data = data_setD)
# W = 102426036, p-value = 0.779



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ THgene, data = data_setD)
# W = 104331794, p-value = 0.2793




### data_setE

data_setE$THgene <- as.factor(data_setE$THgene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ THgene, data = data_setE)
# W = 65374006, p-value = 0.331



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = THgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ THgene, data = data_setE)
# W = 66940338, p-value = 0.5616



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = THgene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = THgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ THgene, data = data_setE)
# W = 66177460, p-value = 0.8603



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = THgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = THgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ THgene, data = data_setE)
# W = 64805294, p-value = 0.1245


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(THgene), scales = "free_y") +
  theme_light()
# 












# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_THgenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_THgenes.R')
