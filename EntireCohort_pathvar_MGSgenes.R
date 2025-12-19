# 28/02/24

# EntireCohort_pathvar_MGSgenes.R


# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)




# MGSgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')


# tag by MGS gene
# could separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")

EntireCohort_pathvar$MGSgene <- NA
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE1 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE2 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE3 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE4 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE5 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE6 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE7 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE8 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1
EntireCohort_pathvar$MGSgene[EntireCohort_pathvar$GENE9 %in% c("ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "GMNN", "MCM3", "MCM5", "MCM7")] <- 1

EntireCohort_pathvar$MGSgene[is.na(EntireCohort_pathvar$MGSgene)] <- 0
EntireCohort_pathvar$MGSgene <- as.numeric(EntireCohort_pathvar$MGSgene)
unique(EntireCohort_pathvar$MGSgene)
# 0 1

table(EntireCohort_pathvar$MGSgene)
# 0      1 
# 448986  20813



# now to make a column for each of the MGS genes
### ORC1
EntireCohort_pathvar$ORC1 <- NA
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE1 %in% c("ORC1")] <- 1
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE2 %in% c("ORC1")] <- 1
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE3 %in% c("ORC1")] <- 1
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE4 %in% c("ORC1")] <- 1
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE5 %in% c("ORC1")] <- 1
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE6 %in% c("ORC1")] <- 1
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE7 %in% c("ORC1")] <- 1
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE8 %in% c("ORC1")] <- 1
EntireCohort_pathvar$ORC1[EntireCohort_pathvar$GENE9 %in% c("ORC1")] <- 1

EntireCohort_pathvar$ORC1[is.na(EntireCohort_pathvar$ORC1)] <- 0
EntireCohort_pathvar$ORC1 <- as.numeric(EntireCohort_pathvar$ORC1)
unique(EntireCohort_pathvar$ORC1)
# 0 1

table(EntireCohort_pathvar$ORC1)
# 0      1 
# 468182   1617

### ORC4
EntireCohort_pathvar$ORC4 <- NA
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE1 %in% c("ORC4")] <- 1
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE2 %in% c("ORC4")] <- 1
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE3 %in% c("ORC4")] <- 1
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE4 %in% c("ORC4")] <- 1
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE5 %in% c("ORC4")] <- 1
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE6 %in% c("ORC4")] <- 1
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE7 %in% c("ORC4")] <- 1
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE8 %in% c("ORC4")] <- 1
EntireCohort_pathvar$ORC4[EntireCohort_pathvar$GENE9 %in% c("ORC4")] <- 1

EntireCohort_pathvar$ORC4[is.na(EntireCohort_pathvar$ORC4)] <- 0
EntireCohort_pathvar$ORC4 <- as.numeric(EntireCohort_pathvar$ORC4)
unique(EntireCohort_pathvar$ORC4)
# 0 1

table(EntireCohort_pathvar$ORC4)
# 0      1 
# 466067   3732 

### ORC6
EntireCohort_pathvar$ORC6 <- NA
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE1 %in% c("ORC6")] <- 1
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE2 %in% c("ORC6")] <- 1
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE3 %in% c("ORC6")] <- 1
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE4 %in% c("ORC6")] <- 1
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE5 %in% c("ORC6")] <- 1
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE6 %in% c("ORC6")] <- 1
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE7 %in% c("ORC6")] <- 1
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE8 %in% c("ORC6")] <- 1
EntireCohort_pathvar$ORC6[EntireCohort_pathvar$GENE9 %in% c("ORC6")] <- 1

EntireCohort_pathvar$ORC6[is.na(EntireCohort_pathvar$ORC6)] <- 0
EntireCohort_pathvar$ORC6 <- as.numeric(EntireCohort_pathvar$ORC6)
unique(EntireCohort_pathvar$ORC6)
# 0 1

table(EntireCohort_pathvar$ORC6)
# 0      1 
# 468884    915

### CDT1
EntireCohort_pathvar$CDT1 <- NA
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE1 %in% c("CDT1")] <- 1
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE2 %in% c("CDT1")] <- 1
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE3 %in% c("CDT1")] <- 1
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE4 %in% c("CDT1")] <- 1
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE5 %in% c("CDT1")] <- 1
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE6 %in% c("CDT1")] <- 1
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE7 %in% c("CDT1")] <- 1
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE8 %in% c("CDT1")] <- 1
EntireCohort_pathvar$CDT1[EntireCohort_pathvar$GENE9 %in% c("CDT1")] <- 1

EntireCohort_pathvar$CDT1[is.na(EntireCohort_pathvar$CDT1)] <- 0
EntireCohort_pathvar$CDT1 <- as.numeric(EntireCohort_pathvar$CDT1)
unique(EntireCohort_pathvar$CDT1)
# 0 1

table(EntireCohort_pathvar$CDT1)
# 0      1 
# 468174   1625

### CDC6
EntireCohort_pathvar$CDC6 <- NA
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE1 %in% c("CDC6")] <- 1
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE2 %in% c("CDC6")] <- 1
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE3 %in% c("CDC6")] <- 1
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE4 %in% c("CDC6")] <- 1
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE5 %in% c("CDC6")] <- 1
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE6 %in% c("CDC6")] <- 1
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE7 %in% c("CDC6")] <- 1
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE8 %in% c("CDC6")] <- 1
EntireCohort_pathvar$CDC6[EntireCohort_pathvar$GENE9 %in% c("CDC6")] <- 1

EntireCohort_pathvar$CDC6[is.na(EntireCohort_pathvar$CDC6)] <- 0
EntireCohort_pathvar$CDC6 <- as.numeric(EntireCohort_pathvar$CDC6)
unique(EntireCohort_pathvar$CDC6)
# 0 1

table(EntireCohort_pathvar$CDC6)
# 0      1 
# 469347    452

### CDC45
EntireCohort_pathvar$CDC45 <- NA
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE1 %in% c("CDC45")] <- 1
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE2 %in% c("CDC45")] <- 1
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE3 %in% c("CDC45")] <- 1
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE4 %in% c("CDC45")] <- 1
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE5 %in% c("CDC45")] <- 1
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE6 %in% c("CDC45")] <- 1
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE7 %in% c("CDC45")] <- 1
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE8 %in% c("CDC45")] <- 1
EntireCohort_pathvar$CDC45[EntireCohort_pathvar$GENE9 %in% c("CDC45")] <- 1

EntireCohort_pathvar$CDC45[is.na(EntireCohort_pathvar$CDC45)] <- 0
EntireCohort_pathvar$CDC45 <- as.numeric(EntireCohort_pathvar$CDC45)
unique(EntireCohort_pathvar$CDC45)
# 0 1

table(EntireCohort_pathvar$CDC45)
# 0      1 
# 467756   2043 

### DONSON
EntireCohort_pathvar$DONSON <- NA
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE1 %in% c("DONSON")] <- 1
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE2 %in% c("DONSON")] <- 1
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE3 %in% c("DONSON")] <- 1
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE4 %in% c("DONSON")] <- 1
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE5 %in% c("DONSON")] <- 1
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE6 %in% c("DONSON")] <- 1
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE7 %in% c("DONSON")] <- 1
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE8 %in% c("DONSON")] <- 1
EntireCohort_pathvar$DONSON[EntireCohort_pathvar$GENE9 %in% c("DONSON")] <- 1

EntireCohort_pathvar$DONSON[is.na(EntireCohort_pathvar$DONSON)] <- 0
EntireCohort_pathvar$DONSON <- as.numeric(EntireCohort_pathvar$DONSON)
unique(EntireCohort_pathvar$DONSON)
# 0 1

table(EntireCohort_pathvar$DONSON)
# 0      1 
# 466255   3544 

### GINS2
EntireCohort_pathvar$GINS2 <- NA
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE1 %in% c("GINS2")] <- 1
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE2 %in% c("GINS2")] <- 1
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE3 %in% c("GINS2")] <- 1
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE4 %in% c("GINS2")] <- 1
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE5 %in% c("GINS2")] <- 1
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE6 %in% c("GINS2")] <- 1
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE7 %in% c("GINS2")] <- 1
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE8 %in% c("GINS2")] <- 1
EntireCohort_pathvar$GINS2[EntireCohort_pathvar$GENE9 %in% c("GINS2")] <- 1

EntireCohort_pathvar$GINS2[is.na(EntireCohort_pathvar$GINS2)] <- 0
EntireCohort_pathvar$GINS2 <- as.numeric(EntireCohort_pathvar$GINS2)
unique(EntireCohort_pathvar$GINS2)
# 0 1

table(EntireCohort_pathvar$GINS2)
# 0      1 
# 469245    554 

### GINS3
EntireCohort_pathvar$GINS3 <- NA
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE1 %in% c("GINS3")] <- 1
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE2 %in% c("GINS3")] <- 1
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE3 %in% c("GINS3")] <- 1
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE4 %in% c("GINS3")] <- 1
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE5 %in% c("GINS3")] <- 1
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE6 %in% c("GINS3")] <- 1
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE7 %in% c("GINS3")] <- 1
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE8 %in% c("GINS3")] <- 1
EntireCohort_pathvar$GINS3[EntireCohort_pathvar$GENE9 %in% c("GINS3")] <- 1

EntireCohort_pathvar$GINS3[is.na(EntireCohort_pathvar$GINS3)] <- 0
EntireCohort_pathvar$GINS3 <- as.numeric(EntireCohort_pathvar$GINS3)
unique(EntireCohort_pathvar$GINS3)
# 0 1

table(EntireCohort_pathvar$GINS3)
# 0      1 
# 469413    386 

### MCM3
EntireCohort_pathvar$MCM3 <- NA
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE1 %in% c("MCM3")] <- 1
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE2 %in% c("MCM3")] <- 1
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE3 %in% c("MCM3")] <- 1
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE4 %in% c("MCM3")] <- 1
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE5 %in% c("MCM3")] <- 1
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE6 %in% c("MCM3")] <- 1
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE7 %in% c("MCM3")] <- 1
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE8 %in% c("MCM3")] <- 1
EntireCohort_pathvar$MCM3[EntireCohort_pathvar$GENE9 %in% c("MCM3")] <- 1

EntireCohort_pathvar$MCM3[is.na(EntireCohort_pathvar$MCM3)] <- 0
EntireCohort_pathvar$MCM3 <- as.numeric(EntireCohort_pathvar$MCM3)
unique(EntireCohort_pathvar$MCM3)
# 0 1

table(EntireCohort_pathvar$MCM3)
# 0      1 
# 467331   2468

### MCM5
EntireCohort_pathvar$MCM5 <- NA
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE1 %in% c("MCM5")] <- 1
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE2 %in% c("MCM5")] <- 1
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE3 %in% c("MCM5")] <- 1
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE4 %in% c("MCM5")] <- 1
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE5 %in% c("MCM5")] <- 1
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE6 %in% c("MCM5")] <- 1
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE7 %in% c("MCM5")] <- 1
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE8 %in% c("MCM5")] <- 1
EntireCohort_pathvar$MCM5[EntireCohort_pathvar$GENE9 %in% c("MCM5")] <- 1

EntireCohort_pathvar$MCM5[is.na(EntireCohort_pathvar$MCM5)] <- 0
EntireCohort_pathvar$MCM5 <- as.numeric(EntireCohort_pathvar$MCM5)
unique(EntireCohort_pathvar$MCM5)
# 0 1

table(EntireCohort_pathvar$MCM5)
# 0      1 
# 467771   2028

### MCM7
EntireCohort_pathvar$MCM7 <- NA
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE1 %in% c("MCM7")] <- 1
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE2 %in% c("MCM7")] <- 1
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE3 %in% c("MCM7")] <- 1
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE4 %in% c("MCM7")] <- 1
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE5 %in% c("MCM7")] <- 1
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE6 %in% c("MCM7")] <- 1
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE7 %in% c("MCM7")] <- 1
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE8 %in% c("MCM7")] <- 1
EntireCohort_pathvar$MCM7[EntireCohort_pathvar$GENE9 %in% c("MCM7")] <- 1

EntireCohort_pathvar$MCM7[is.na(EntireCohort_pathvar$MCM7)] <- 0
EntireCohort_pathvar$MCM7 <- as.numeric(EntireCohort_pathvar$MCM7)
unique(EntireCohort_pathvar$MCM7)
# 0 1

table(EntireCohort_pathvar$MCM7)
# 0      1 
# 467936   1863


# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)








# MGSgenes Incidence --------------------------------------------

# Start things out by repeating the logistic regression (cancer incidence) for the MGS genes

# now that we have a dummy variable for people with MGS gene variants (MGSgene),
# we can do logistic regression analysis against cancer_ROOC_dummy

# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_MGSgene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_MGSgene_, "pathvar_MGSgene_.csv", row.names=FALSE)
# system('dx upload pathvar_MGSgene_.csv --path Emily-folder/results_2/pathvar_MGSgene_.csv')


# MGSgene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ MGSgene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_MGSgene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_incidence_1, "pathvar_MGSgene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_incidence_1.csv --path Emily-folder/results_2/pathvar_MGSgene_incidence_1.csv')


# "ORC1", "ORC4", "ORC6", "CDT1", "CDC6", "CDC45", "DONSON", "GINS2", "GINS3", "MCM3", "MCM5", "MCM7"
LR_2 <- glm(formula = cancer_ROOC_dummy ~ ORC1+ORC4+ORC6+CDT1+CDC6+CDC45+DONSON+GINS2+GINS3+MCM3+MCM5+MCM7, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_MGSgene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_incidence_2, "pathvar_MGSgene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_incidence_2.csv --path Emily-folder/results_2/pathvar_MGSgene_incidence_2.csv')


# won't save each individual gene - just here for reference within R
# "ORC1"
LR_3 <- glm(formula = cancer_ROOC_dummy ~ ORC1, data = EntireCohort_pathvar, family = binomial)
summary(LR_3)
exp(cbind(OR = coef(LR_3), confint(LR_3)))
# "ORC4"
LR_4 <- glm(formula = cancer_ROOC_dummy ~ ORC4, data = EntireCohort_pathvar, family = binomial)
summary(LR_4)
exp(cbind(OR = coef(LR_4), confint(LR_4)))
# "ORC6"
LR_5 <- glm(formula = cancer_ROOC_dummy ~ ORC6, data = EntireCohort_pathvar, family = binomial)
summary(LR_5)
exp(cbind(OR = coef(LR_5), confint(LR_5)))
# "CDT1"
LR_6 <- glm(formula = cancer_ROOC_dummy ~ CDT1, data = EntireCohort_pathvar, family = binomial)
summary(LR_6)
exp(cbind(OR = coef(LR_6), confint(LR_6)))
# "CDC6"
LR_7 <- glm(formula = cancer_ROOC_dummy ~ CDC6, data = EntireCohort_pathvar, family = binomial)
summary(LR_7)
exp(cbind(OR = coef(LR_7), confint(LR_7)))
# "CDC45"
LR_8 <- glm(formula = cancer_ROOC_dummy ~ CDC45, data = EntireCohort_pathvar, family = binomial)
summary(LR_8)
exp(cbind(OR = coef(LR_8), confint(LR_8)))
# "DONSON"
LR_9 <- glm(formula = cancer_ROOC_dummy ~ DONSON, data = EntireCohort_pathvar, family = binomial)
summary(LR_9)
exp(cbind(OR = coef(LR_9), confint(LR_9)))
# "GINS2"
LR_10 <- glm(formula = cancer_ROOC_dummy ~ GINS2, data = EntireCohort_pathvar, family = binomial)
summary(LR_10)
exp(cbind(OR = coef(LR_10), confint(LR_10)))
# "GINS3"
LR_11 <- glm(formula = cancer_ROOC_dummy ~ GINS3, data = EntireCohort_pathvar, family = binomial)
summary(LR_11)
exp(cbind(OR = coef(LR_11), confint(LR_11)))
# "MCM3"
LR_12 <- glm(formula = cancer_ROOC_dummy ~ MCM3, data = EntireCohort_pathvar, family = binomial)
summary(LR_12)
exp(cbind(OR = coef(LR_12), confint(LR_12)))
#  "MCM5"
LR_13 <- glm(formula = cancer_ROOC_dummy ~ MCM5, data = EntireCohort_pathvar, family = binomial)
summary(LR_13)
exp(cbind(OR = coef(LR_13), confint(LR_13)))
# "MCM7"
LR_14 <- glm(formula = cancer_ROOC_dummy ~ MCM7, data = EntireCohort_pathvar, family = binomial)
summary(LR_14)
exp(cbind(OR = coef(LR_14), confint(LR_14)))








### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "MGSgene", "ORC1", "ORC4", "ORC6", "CDT1", "CDC6", 
                                     "CDC45", "DONSON", "GINS2", "GINS3", "MCM3", "MCM5", "MCM7")]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     19




# MGSgene by itself
LR_15 <- glm(formula = cancer_ROOC_dummy ~ MGSgene, data = data_setA, family = binomial)
# summary(LR_15)
# exp(cbind(OR = coef(LR_15), confint(LR_15)))
pathvar_MGSgene_incidence_3 <- tidy(LR_15, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_incidence_3, "pathvar_MGSgene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_incidence_3.csv --path Emily-folder/results_2/pathvar_MGSgene_incidence_3.csv')


# age + ethnicity + smoking + sex + MGSgene
LR_16 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, 
             data = data_setA, family = binomial)
# summary(LR_16)
# exp(cbind(OR = coef(LR_16), confint(LR_16)))
pathvar_MGSgene_incidence_4 <- tidy(LR_16, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_incidence_4, "pathvar_MGSgene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_incidence_4.csv --path Emily-folder/results_2/pathvar_MGSgene_incidence_4.csv')



# GOIs by themselves
LR_17 <- glm(formula = cancer_ROOC_dummy ~ ORC1+ORC4+ORC6+CDT1+CDC6+CDC45+DONSON+GINS2+GINS3+MCM3+MCM5+MCM7, 
             data = data_setA, family = binomial)
# summary(LR_17)
# exp(cbind(OR = coef(LR_17), confint(LR_17)))
pathvar_MGSgene_incidence_5 <- tidy(LR_17, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_incidence_5, "pathvar_MGSgene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_incidence_5.csv --path Emily-folder/results_2/pathvar_MGSgene_incidence_5.csv')



# age + ethnicity + smoking + sex + GOIs 
LR_18 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ORC1+ORC4+ORC6+CDT1+CDC6+CDC45+DONSON+GINS2+GINS3+MCM3+MCM5+MCM7, 
             data = data_setA, family = binomial)
# summary(LR_18)
# exp(cbind(OR = coef(LR_18), confint(LR_18)))
pathvar_MGSgene_incidence_6 <- tidy(LR_18, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_incidence_6, "pathvar_MGSgene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_incidence_6.csv --path Emily-folder/results_2/pathvar_MGSgene_incidence_6.csv')





# wont save each individual gene
# age + ethnicity + smoking + sex + ORC6
LR_19 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment + Ethnicity_White + 
               Ever_smoked_dummy + Sex_Male + ORC6, 
             data = data_setA, family = binomial)
summary(LR_19)
exp(cbind(OR = coef(LR_19), confint(LR_19)))

# age + ethnicity + smoking + sex + CDT1
LR_20 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment + Ethnicity_White + 
               Ever_smoked_dummy + Sex_Male + CDT1, 
             data = data_setA, family = binomial)
summary(LR_20)
exp(cbind(OR = coef(LR_20), confint(LR_20)))

# age + ethnicity + smoking + sex + DONSON
LR_21 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment + Ethnicity_White + 
               Ever_smoked_dummy + Sex_Male + DONSON, 
             data = data_setA, family = binomial)
summary(LR_21)
exp(cbind(OR = coef(LR_21), confint(LR_21)))







### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "MGSgene", "ORC1", "ORC4", "ORC6", "CDT1", "CDC6", 
                                     "CDC45", "DONSON", "GINS2", "GINS3", "MCM3", "MCM5", "MCM7")]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     24




# MGSgene only
LR_22 <- glm(formula = cancer_ROOC_dummy ~ MGSgene, data = data_setB, family = binomial)
# summary(LR_22)
# exp(cbind(OR = coef(LR_22), confint(LR_22)))
pathvar_MGSgene_incidence_7 <- tidy(LR_22, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_incidence_7, "pathvar_MGSgene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_incidence_7.csv --path Emily-folder/results_2/pathvar_MGSgene_incidence_7.csv')


# all variables + MGSgene
LR_23 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+MGSgene, 
             data = data_setB, family = binomial)
# summary(LR_23)
# exp(cbind(OR = coef(LR_23), confint(LR_23)))
pathvar_MGSgene_incidence_8 <- tidy(LR_23, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_incidence_8, "pathvar_MGSgene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_incidence_8.csv --path Emily-folder/results_2/pathvar_MGSgene_incidence_8.csv')



# GOIs alone
LR_24 <- glm(formula = cancer_ROOC_dummy ~ ORC1+ORC4+ORC6+CDT1+CDC6+CDC45+DONSON+GINS2+GINS3+MCM3+MCM5+MCM7, 
             data = data_setB, family = binomial)
# summary(LR_24)
# exp(cbind(OR = coef(LR_24), confint(LR_24)))
pathvar_MGSgene_incidence_9 <- tidy(LR_24, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_incidence_9, "pathvar_MGSgene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_incidence_9.csv --path Emily-folder/results_2/pathvar_MGSgene_incidence_9.csv')


# all variables + GOIs
LR_25 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+ORC1+ORC4+ORC6+CDT1+CDC6+CDC45+DONSON+GINS2+GINS3+MCM3+MCM5+MCM7, 
             data = data_setB, family = binomial)
# summary(LR_25)
# exp(cbind(OR = coef(LR_25), confint(LR_25)))
pathvar_MGSgene_incidence_10 <- tidy(LR_25, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_incidence_10, "pathvar_MGSgene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_incidence_10.csv --path Emily-folder/results_2/pathvar_MGSgene_incidence_10.csv')












### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(MGSgene, cancer_ROOC_dummy))			
fisher.test(with(EntireCohort_pathvar, table(MGSgene, cancer_ROOC_dummy))	)
# p-value = 0.3717

with(EntireCohort_pathvar, table(MGSgene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(MGSgene, Ethnicity_White)))
# p-value < 2.2e-16
# 95 percent confidence interval:
#   0.5598989 0.6201634
# sample estimates:
#   odds ratio 
# 0.5891156

with(EntireCohort_pathvar, table(MGSgene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(MGSgene, Ever_smoked_dummy)))
# p-value = 0.1315

with(EntireCohort_pathvar, table(MGSgene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(MGSgene, Sex_Male)))
# p-value = 0.87

with(EntireCohort_pathvar, table(MGSgene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(MGSgene, University_education)))
# p-value = 0.3082

with(EntireCohort_pathvar, table(MGSgene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(MGSgene, Annual_income)))
# p-value = 0.7851

with(EntireCohort_pathvar, table(MGSgene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(MGSgene, Smoker_previous)))
# p-value = 0.3354

with(EntireCohort_pathvar, table(MGSgene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(MGSgene, Smoker_current)))
# p-value = 0.5631

### data_setA
with(data_setA, table(MGSgene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(MGSgene, cancer_ROOC_dummy)))
# p-value = 0.3142

with(data_setA, table(MGSgene, Ethnicity_White))			
fisher.test(with(data_setA, table(MGSgene, Ethnicity_White)))
# p-value < 2.2e-16
# 95 percent confidence interval:
#   0.5584351 0.6187701
# sample estimates:
#   odds ratio 
# 0.5876681 

with(data_setA, table(MGSgene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(MGSgene, Ever_smoked_dummy)))
# p-value = 0.1743

with(data_setA, table(MGSgene, Sex_Male))			
fisher.test(with(data_setA, table(MGSgene, Sex_Male)))
# p-value = 0.9373



### data_setB

with(data_setB, table(MGSgene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(MGSgene, cancer_ROOC_dummy)))
# p-value = 0.3166

with(data_setB, table(MGSgene, Ethnicity_White))			
fisher.test(with(data_setB, table(MGSgene, Ethnicity_White)))
# p-value < 2.2e-16
# 95 percent confidence interval:
#   0.5531270 0.6257177
# sample estimates:
#   odds ratio 
# 0.5881111

with(data_setB, table(MGSgene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(MGSgene, Ever_smoked_dummy)))
# p-value = 0.605

with(data_setB, table(MGSgene, Sex_Male))			
fisher.test(with(data_setB, table(MGSgene, Sex_Male)))
# p-value = 0.8619

with(data_setB, table(MGSgene, University_education))			
fisher.test(with(data_setB, table(MGSgene, University_education)))
# p-value = 0.3227

with(data_setB, table(MGSgene, Annual_income))			
fisher.test(with(data_setB, table(MGSgene, Annual_income)))
# p-value = 0.7813







### can plot the continuous variables: 

# make MGSgene into factor in  new df
df_to_plot <- EntireCohort_pathvar
df_to_plot$MGSgene <- as.factor(df_to_plot$MGSgene)

# try a violin plot

# Age_at_recruitment  
df_to_plot %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_to_plot %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ MGSgene, data = df_to_plot)
# W = 4753212598, p-value = 2.344e-05


# BMI  
df_to_plot %>% 
  ggplot(mapping = aes(x = MGSgene, y = BMI)) +
  geom_violin() +
  theme_light()

df_to_plot %>% 
  ggplot(mapping = aes(x = MGSgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ MGSgene, data = df_to_plot)
# W = 4612515150, p-value = 0.2848


# Standing_height 
df_to_plot %>% 
  ggplot(mapping = aes(x = MGSgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_to_plot %>% 
  ggplot(mapping = aes(x = MGSgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ MGSgene, data = df_to_plot)
# W = 4701056556, p-value = 0.001845


# make this a bar plot
# Weekly_exercise
df_to_plot %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(MGSgene), scales = "free_y") +
  theme_light()




### data_setA

# make MGSgene into factor
data_setA$MGSgene <- as.factor(data_setA$MGSgene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setA %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ MGSgene, data = data_setA)
# W = 4671579095, p-value = 2.082e-05




### data_setB

# make MGSgene into factor
data_setB$MGSgene <- as.factor(data_setB$MGSgene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setB %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ MGSgene, data = data_setB)
# W = 3109704648, p-value = 3.291e-05

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = MGSgene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setB %>% 
  ggplot(mapping = aes(x = MGSgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ MGSgene, data = data_setB)
# W = 3040204800, p-value = 0.3999


# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = MGSgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setB %>% 
  ggplot(mapping = aes(x = MGSgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ MGSgene, data = data_setB)
# W = 3086444038, p-value = 0.01313



# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(MGSgene), scales = "free_y") +
  theme_light()





### fisher test for GOIs of interest

# ORC6
with(EntireCohort_pathvar, table(ORC6, cancer_ROOC_dummy))			
fisher.test(with(EntireCohort_pathvar, table(ORC6, cancer_ROOC_dummy)))
# p-value = 0.01404
# 95 percent confidence interval:
#   1.036433 1.395434
# sample estimates:
#   odds ratio 
# 1.204159

with(EntireCohort_pathvar, table(ORC6, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(ORC6, Ethnicity_White)))
# p-value = 0.0002618
# 95 percent confidence interval:
#   1.329431 3.071234
# sample estimates:
#   odds ratio 
# 1.976236

with(EntireCohort_pathvar, table(ORC6, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(ORC6, Ever_smoked_dummy)))
# p-value = 0.2364

with(EntireCohort_pathvar, table(ORC6, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(ORC6, Sex_Male)))
# p-value = 0.8681

with(EntireCohort_pathvar, table(ORC6, University_education))			
fisher.test(with(EntireCohort_pathvar, table(ORC6, University_education)))
# p-value = 0.07035
# 95 percent confidence interval:
#   0.755793 1.010516
# sample estimates:
#   odds ratio 
# 0.8748736 

with(EntireCohort_pathvar, table(ORC6, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(ORC6, Annual_income)))
# p-value = 0.4851


# DONSON
with(EntireCohort_pathvar, table(DONSON, cancer_ROOC_dummy))			
fisher.test(with(EntireCohort_pathvar, table(DONSON, cancer_ROOC_dummy)))
# p-value = 0.001589
# 95 percent confidence interval:
#   0.8093713 0.9525781
# sample estimates:
#   odds ratio 
# 0.878484

with(EntireCohort_pathvar, table(DONSON, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(DONSON, Ethnicity_White)))
# p-value < 2.2e-16
# 95 percent confidence interval:
#   0.5023671 0.6361491
# sample estimates:
#   odds ratio 
# 0.5645347 

with(EntireCohort_pathvar, table(DONSON, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(DONSON, Ever_smoked_dummy)))
# p-value = 0.01371
# 95 percent confidence interval:
#   0.8589786 0.9834051
# sample estimates:
#   odds ratio 
# 0.9190196 

with(EntireCohort_pathvar, table(DONSON, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(DONSON, Sex_Male)))
# p-value = 0.879

with(EntireCohort_pathvar, table(DONSON, University_education))			
fisher.test(with(EntireCohort_pathvar, table(DONSON, University_education)))
# p-value = 0.5507

with(EntireCohort_pathvar, table(DONSON, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(DONSON, Annual_income)))
# p-value = 0.5871









# could do plots for interesting GOIs














# MGSgenes Type ---------------------------------------------------

# make sure the execute the code in the 'MGSgenes' section


## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses


# LR for each cancer category




LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_MGSgene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_type_1, "pathvar_MGSgene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_type_1.csv --path Emily-folder/results_2/pathvar_MGSgene_type_1.csv')




LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_MGSgene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_type_2, "pathvar_MGSgene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_type_2.csv --path Emily-folder/results_2/pathvar_MGSgene_type_2.csv')




LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_MGSgene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_type_3, "pathvar_MGSgene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_type_3.csv --path Emily-folder/results_2/pathvar_MGSgene_type_3.csv')




LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_MGSgene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_type_4, "pathvar_MGSgene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_type_4.csv --path Emily-folder/results_2/pathvar_MGSgene_type_4.csv')




LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_MGSgene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_type_5, "pathvar_MGSgene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_type_5.csv --path Emily-folder/results_2/pathvar_MGSgene_type_5.csv')




LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_MGSgene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_type_6, "pathvar_MGSgene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_type_6.csv --path Emily-folder/results_2/pathvar_MGSgene_type_6.csv')




LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_MGSgene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_type_7, "pathvar_MGSgene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_type_7.csv --path Emily-folder/results_2/pathvar_MGSgene_type_7.csv')




LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_MGSgene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_type_8, "pathvar_MGSgene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_type_8.csv --path Emily-folder/results_2/pathvar_MGSgene_type_8.csv')




LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_MGSgene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_type_9, "pathvar_MGSgene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_type_9.csv --path Emily-folder/results_2/pathvar_MGSgene_type_9.csv')




LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_MGSgene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_type_10, "pathvar_MGSgene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_type_10.csv --path Emily-folder/results_2/pathvar_MGSgene_type_10.csv')




LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_MGSgene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_type_11, "pathvar_MGSgene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_type_11.csv --path Emily-folder/results_2/pathvar_MGSgene_type_11.csv')




LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_MGSgene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_type_12, "pathvar_MGSgene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_type_12.csv --path Emily-folder/results_2/pathvar_MGSgene_type_12.csv')




LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_MGSgene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_type_13, "pathvar_MGSgene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_type_13.csv --path Emily-folder/results_2/pathvar_MGSgene_type_13.csv')




LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_MGSgene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_type_14, "pathvar_MGSgene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_type_14.csv --path Emily-folder/results_2/pathvar_MGSgene_type_14.csv')




LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_MGSgene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_MGSgene_type_15, "pathvar_MGSgene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_type_15.csv --path Emily-folder/results_2/pathvar_MGSgene_type_15.csv')













# MGSgenes Age  -------------------------------------------------------


# make sure the execute the code in the 'MGS genes' section



### linear regression
# "MGSgene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ MGSgene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_MGSgene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_MGSgene_age_1, "pathvar_MGSgene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_age_1.csv --path Emily-folder/results_2/pathvar_MGSgene_age_1.csv')



# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ ORC1+ORC4+ORC6+CDT1+CDC6+CDC45+DONSON+GINS2+GINS3+MCM3+MCM5+MCM7, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_MGSgene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_MGSgene_age_2, "pathvar_MGSgene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_age_2.csv --path Emily-folder/results_2/pathvar_MGSgene_age_2.csv')



# won't save the individual gene analyses as they don't differ from when they're together
LM_3 <- lm(Age_at_cancer_diagnosis_earliest ~ ORC1, data = EntireCohort_pathvar)
summary(LM_3)

LM_4 <- lm(Age_at_cancer_diagnosis_earliest ~ ORC4, data = EntireCohort_pathvar)
summary(LM_4)

LM_5 <- lm(Age_at_cancer_diagnosis_earliest ~ ORC6, data = EntireCohort_pathvar)
summary(LM_5)

LM_6 <- lm(Age_at_cancer_diagnosis_earliest ~ CDT1, data = EntireCohort_pathvar)
summary(LM_6)

LM_7 <- lm(Age_at_cancer_diagnosis_earliest ~ CDC6, data = EntireCohort_pathvar)
summary(LM_7)

LM_8 <- lm(Age_at_cancer_diagnosis_earliest ~ CDC45, data = EntireCohort_pathvar)
summary(LM_8)

LM_9 <- lm(Age_at_cancer_diagnosis_earliest ~ DONSON, data = EntireCohort_pathvar)
summary(LM_9)

LM_10 <- lm(Age_at_cancer_diagnosis_earliest ~ GINS2, data = EntireCohort_pathvar)
summary(LM_10)

LM_11 <- lm(Age_at_cancer_diagnosis_earliest ~ GINS3, data = EntireCohort_pathvar)
summary(LM_11)

LM_12 <- lm(Age_at_cancer_diagnosis_earliest ~ MCM3, data = EntireCohort_pathvar)
summary(LM_12)

LM_13 <- lm(Age_at_cancer_diagnosis_earliest ~ MCM5, data = EntireCohort_pathvar)
summary(LM_13)

LM_14 <- lm(Age_at_cancer_diagnosis_earliest ~ MCM7, data = EntireCohort_pathvar)
summary(LM_14)







## multivariate analysis



# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "MGSgene", "ORC1", "ORC4", "ORC6", "CDT1", "CDC6", 
                                     "CDC45", "DONSON", "GINS2", "GINS3", "MCM3", "MCM5", "MCM7")]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     19


data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "MGSgene", "ORC1", "ORC4", "ORC6", "CDT1", "CDC6", 
                                     "CDC45", "DONSON", "GINS2", "GINS3", "MCM3", "MCM5", "MCM7")]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    24



### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ MGSgene, data = data_setD)
# summary(LM_15)
pathvar_MGSgene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_MGSgene_age_3, "pathvar_MGSgene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_age_3.csv --path Emily-folder/results_2/pathvar_MGSgene_age_3.csv')



LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSgene, data = data_setD)
# summary(LM_16)
# confint(LM_16)
pathvar_MGSgene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_MGSgene_age_4, "pathvar_MGSgene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_age_4.csv --path Emily-folder/results_2/pathvar_MGSgene_age_4.csv')



LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ ORC1+ORC4+ORC6+CDT1+CDC6+CDC45+DONSON+GINS2+GINS3+MCM3+MCM5+MCM7, data = data_setD)
# summary(LM_17)
pathvar_MGSgene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_MGSgene_age_5, "pathvar_MGSgene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_age_5.csv --path Emily-folder/results_2/pathvar_MGSgene_age_5.csv')



LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+ORC1+ORC4+ORC6+CDT1+CDC6+CDC45+DONSON+GINS2+GINS3+MCM3+MCM5+MCM7, data = data_setD)
# summary(LM_18)
# confint(LM_18)
pathvar_MGSgene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_MGSgene_age_6, "pathvar_MGSgene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_age_6.csv --path Emily-folder/results_2/pathvar_MGSgene_age_6.csv')





### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ MGSgene, data = data_setE)
# summary(LM_15)
pathvar_MGSgene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_MGSgene_age_7, "pathvar_MGSgene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_age_7.csv --path Emily-folder/results_2/pathvar_MGSgene_age_7.csv')



LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+MGSgene, data = data_setE)
# summary(LM_16)
pathvar_MGSgene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_MGSgene_age_8, "pathvar_MGSgene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_age_8.csv --path Emily-folder/results_2/pathvar_MGSgene_age_8.csv')



LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ ORC1+ORC4+ORC6+CDT1+CDC6+CDC45+DONSON+GINS2+GINS3+MCM3+MCM5+MCM7, data = data_setE)
# summary(LM_17)
pathvar_MGSgene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_MGSgene_age_9, "pathvar_MGSgene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_age_9.csv --path Emily-folder/results_2/pathvar_MGSgene_age_9.csv')



LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+ORC1+ORC4+ORC6+CDT1+CDC6+CDC45+DONSON+GINS2+GINS3+MCM3+MCM5+MCM7, data = data_setE)
# summary(LM_18)
pathvar_MGSgene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_MGSgene_age_10, "pathvar_MGSgene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_MGSgene_age_10.csv --path Emily-folder/results_2/pathvar_MGSgene_age_10.csv')










## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar[complete.cases(EntireCohort_pathvar$Age_at_cancer_diagnosis_earliest), ]



with(df_cancer_diagnosis, table(MGSgene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(MGSgene, Ethnicity_White)))
# p-value = 2.358e-13
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.4921333 0.6558974
# sample estimates:
#   odds ratio 
# 0.5669916

with(df_cancer_diagnosis, table(MGSgene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(MGSgene, Ever_smoked_dummy)))
# p-value = 1

with(df_cancer_diagnosis, table(MGSgene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(MGSgene, Sex_Male)))
# p-value = 0.791

with(df_cancer_diagnosis, table(MGSgene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(MGSgene, University_education)))
# p-value = 0.3263

with(df_cancer_diagnosis, table(MGSgene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(MGSgene, Annual_income)))
# p-value = 0.03971
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.003198 1.169414
# sample estimates:
#   odds ratio 
# 1.083555


### data_setD
with(data_setD, table(MGSgene, Ethnicity_White))			
fisher.test(with(data_setD, table(MGSgene, Ethnicity_White)))
# p-value = 2.483e-13
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.4910620 0.6553223
# sample estimates:
#   odds ratio 
# 0.5661067

with(data_setD, table(MGSgene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(MGSgene, Ever_smoked_dummy)))
# p-value = 1

with(data_setD, table(MGSgene, Sex_Male))			
fisher.test(with(data_setD, table(MGSgene, Sex_Male)))
# p-value = 0.894


### data_setE
with(data_setE, table(MGSgene, Ethnicity_White))			
fisher.test(with(data_setE, table(MGSgene, Ethnicity_White)))
# p-value = 2.682e-08
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.4924522 0.7043897
# sample estimates:
#   odds ratio 
# 0.5870354

with(data_setE, table(MGSgene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(MGSgene, Ever_smoked_dummy)))
# p-value = 0.8101

with(data_setE, table(MGSgene, Sex_Male))			
fisher.test(with(data_setE, table(MGSgene, Sex_Male)))
# p-value = 0.5519

with(data_setE, table(MGSgene, University_education))			
fisher.test(with(data_setE, table(MGSgene, University_education)))
# p-value = 0.2697

with(data_setE, table(MGSgene, Annual_income))			
fisher.test(with(data_setE, table(MGSgene, Annual_income)))
# p-value = 0.0559








## make plots for the linear regression data


# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 
# df_cancer_diagnosis



# make MGSgene into factor
df_cancer_diagnosis$MGSgene <- as.factor(df_cancer_diagnosis$MGSgene)

# try a violin plot

# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ MGSgene, data = df_cancer_diagnosis)
# 

# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ MGSgene, data = df_cancer_diagnosis)
# W = 265631038, p-value = 0.0001928


# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = MGSgene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = MGSgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ MGSgene, data = df_cancer_diagnosis)
# W = 252502760, p-value = 0.1591


# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = MGSgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = MGSgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ MGSgene, data = df_cancer_diagnosis)
# W = 255277413, p-value = 0.6812


# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(MGSgene), scales = "free_y") +
  theme_light()




### data_setD

data_setD$MGSgene <- as.factor(data_setD$MGSgene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ MGSgene, data = data_setD)
# 


# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ MGSgene, data = data_setD)
# W = 260769456, p-value = 0.0001513





### data_setE

data_setE$MGSgene <- as.factor(data_setE$MGSgene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ MGSgene, data = data_setE)
# 



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = MGSgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ MGSgene, data = data_setE)
# W = 166814025, p-value = 0.003118


# BMI
data_setE %>% 
  ggplot(mapping = aes(x = MGSgene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = MGSgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ MGSgene, data = data_setE)
# W = 160135266, p-value = 0.1721


# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = MGSgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = MGSgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ MGSgene, data = data_setE)
# W = 161067825, p-value = 0.4456


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(MGSgene), scales = "free_y") +
  theme_light()
















# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_MGSgenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_MGSgenes.R')
