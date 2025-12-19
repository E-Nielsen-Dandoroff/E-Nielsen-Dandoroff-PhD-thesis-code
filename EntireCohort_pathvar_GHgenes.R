# 18/06/2024

# EntireCohort_pathvar_GHgenes.R

# this script is for the Growth Hormone associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)



# GHgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')


# go thr all the genes to check that they're present
grep("STAT5B", EntireCohort_pathvar$GENE)
# all seem to be present

# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$GHgene <- NA
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE1 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE2 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE3 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE4 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE5 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE6 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE7 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE8 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1
EntireCohort_pathvar$GHgene[EntireCohort_pathvar$GENE9 %in% c('BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B')] <- 1


EntireCohort_pathvar$GHgene[is.na(EntireCohort_pathvar$GHgene)] <- 0
EntireCohort_pathvar$GHgene <- as.numeric(EntireCohort_pathvar$GHgene)
unique(EntireCohort_pathvar$GHgene)
# 0 1

table(EntireCohort_pathvar$GHgene)
#      0      1 
# 456368  13431



# now to make a column for each of the GHgenes
### 'BTK'
EntireCohort_pathvar$BTK <- NA
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE1 %in% c("BTK")] <- 1
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE2 %in% c("BTK")] <- 1
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE3 %in% c("BTK")] <- 1
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE4 %in% c("BTK")] <- 1
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE5 %in% c("BTK")] <- 1
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE6 %in% c("BTK")] <- 1
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE7 %in% c("BTK")] <- 1
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE8 %in% c("BTK")] <- 1
EntireCohort_pathvar$BTK[EntireCohort_pathvar$GENE9 %in% c("BTK")] <- 1

EntireCohort_pathvar$BTK[is.na(EntireCohort_pathvar$BTK)] <- 0
EntireCohort_pathvar$BTK <- as.numeric(EntireCohort_pathvar$BTK)
unique(EntireCohort_pathvar$BTK)
# 0 1

table(EntireCohort_pathvar$BTK)
# 0      1 
# 469710     89





### 'GH1'
EntireCohort_pathvar$GH1 <- NA
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE1 %in% c("GH1")] <- 1
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE2 %in% c("GH1")] <- 1
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE3 %in% c("GH1")] <- 1
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE4 %in% c("GH1")] <- 1
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE5 %in% c("GH1")] <- 1
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE6 %in% c("GH1")] <- 1
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE7 %in% c("GH1")] <- 1
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE8 %in% c("GH1")] <- 1
EntireCohort_pathvar$GH1[EntireCohort_pathvar$GENE9 %in% c("GH1")] <- 1

EntireCohort_pathvar$GH1[is.na(EntireCohort_pathvar$GH1)] <- 0
EntireCohort_pathvar$GH1 <- as.numeric(EntireCohort_pathvar$GH1)
unique(EntireCohort_pathvar$GH1)
# 0 1

table(EntireCohort_pathvar$GH1)
# 0      1 
# 469552    247




### 'GHR'
EntireCohort_pathvar$GHR <- NA
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE1 %in% c("GHR")] <- 1
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE2 %in% c("GHR")] <- 1
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE3 %in% c("GHR")] <- 1
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE4 %in% c("GHR")] <- 1
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE5 %in% c("GHR")] <- 1
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE6 %in% c("GHR")] <- 1
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE7 %in% c("GHR")] <- 1
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE8 %in% c("GHR")] <- 1
EntireCohort_pathvar$GHR[EntireCohort_pathvar$GENE9 %in% c("GHR")] <- 1

EntireCohort_pathvar$GHR[is.na(EntireCohort_pathvar$GHR)] <- 0
EntireCohort_pathvar$GHR <- as.numeric(EntireCohort_pathvar$GHR)
unique(EntireCohort_pathvar$GHR)
# 0 1

table(EntireCohort_pathvar$GHR)
# 0      1 
# 467645   2154




### 'GHRH'
EntireCohort_pathvar$GHRH <- NA
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE1 %in% c("GHRH")] <- 1
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE2 %in% c("GHRH")] <- 1
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE3 %in% c("GHRH")] <- 1
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE4 %in% c("GHRH")] <- 1
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE5 %in% c("GHRH")] <- 1
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE6 %in% c("GHRH")] <- 1
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE7 %in% c("GHRH")] <- 1
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE8 %in% c("GHRH")] <- 1
EntireCohort_pathvar$GHRH[EntireCohort_pathvar$GENE9 %in% c("GHRH")] <- 1

EntireCohort_pathvar$GHRH[is.na(EntireCohort_pathvar$GHRH)] <- 0
EntireCohort_pathvar$GHRH <- as.numeric(EntireCohort_pathvar$GHRH)
unique(EntireCohort_pathvar$GHRH)
# 0 1

table(EntireCohort_pathvar$GHRH)
# 0      1 
# 469712     87




### 'GHRHR'
EntireCohort_pathvar$GHRHR <- NA
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE1 %in% c("GHRHR")] <- 1
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE2 %in% c("GHRHR")] <- 1
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE3 %in% c("GHRHR")] <- 1
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE4 %in% c("GHRHR")] <- 1
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE5 %in% c("GHRHR")] <- 1
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE6 %in% c("GHRHR")] <- 1
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE7 %in% c("GHRHR")] <- 1
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE8 %in% c("GHRHR")] <- 1
EntireCohort_pathvar$GHRHR[EntireCohort_pathvar$GENE9 %in% c("GHRHR")] <- 1

EntireCohort_pathvar$GHRHR[is.na(EntireCohort_pathvar$GHRHR)] <- 0
EntireCohort_pathvar$GHRHR <- as.numeric(EntireCohort_pathvar$GHRHR)
unique(EntireCohort_pathvar$GHRHR)
# 0 1

table(EntireCohort_pathvar$GHRHR)
# 0      1 
# 469138    661




### 'GHSR'
EntireCohort_pathvar$GHSR <- NA
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE1 %in% c("GHSR")] <- 1
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE2 %in% c("GHSR")] <- 1
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE3 %in% c("GHSR")] <- 1
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE4 %in% c("GHSR")] <- 1
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE5 %in% c("GHSR")] <- 1
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE6 %in% c("GHSR")] <- 1
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE7 %in% c("GHSR")] <- 1
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE8 %in% c("GHSR")] <- 1
EntireCohort_pathvar$GHSR[EntireCohort_pathvar$GENE9 %in% c("GHSR")] <- 1

EntireCohort_pathvar$GHSR[is.na(EntireCohort_pathvar$GHSR)] <- 0
EntireCohort_pathvar$GHSR <- as.numeric(EntireCohort_pathvar$GHSR)
unique(EntireCohort_pathvar$GHSR)
# 0 1

table(EntireCohort_pathvar$GHSR)
# 0      1 
# 468753   1046




### 'HESX1'
EntireCohort_pathvar$HESX1 <- NA
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE1 %in% c("HESX1")] <- 1
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE2 %in% c("HESX1")] <- 1
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE3 %in% c("HESX1")] <- 1
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE4 %in% c("HESX1")] <- 1
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE5 %in% c("HESX1")] <- 1
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE6 %in% c("HESX1")] <- 1
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE7 %in% c("HESX1")] <- 1
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE8 %in% c("HESX1")] <- 1
EntireCohort_pathvar$HESX1[EntireCohort_pathvar$GENE9 %in% c("HESX1")] <- 1

EntireCohort_pathvar$HESX1[is.na(EntireCohort_pathvar$HESX1)] <- 0
EntireCohort_pathvar$HESX1 <- as.numeric(EntireCohort_pathvar$HESX1)
unique(EntireCohort_pathvar$HESX1)
# 0 1

table(EntireCohort_pathvar$HESX1)
# 0      1 
# 469289    510




### 'IGF1'
EntireCohort_pathvar$IGF1 <- NA
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE1 %in% c("IGF1")] <- 1
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE2 %in% c("IGF1")] <- 1
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE3 %in% c("IGF1")] <- 1
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE4 %in% c("IGF1")] <- 1
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE5 %in% c("IGF1")] <- 1
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE6 %in% c("IGF1")] <- 1
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE7 %in% c("IGF1")] <- 1
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE8 %in% c("IGF1")] <- 1
EntireCohort_pathvar$IGF1[EntireCohort_pathvar$GENE9 %in% c("IGF1")] <- 1

EntireCohort_pathvar$IGF1[is.na(EntireCohort_pathvar$IGF1)] <- 0
EntireCohort_pathvar$IGF1 <- as.numeric(EntireCohort_pathvar$IGF1)
unique(EntireCohort_pathvar$IGF1)
# 0 1

table(EntireCohort_pathvar$IGF1)
# 0      1 
# 469390    409




### 'IGF1R'
EntireCohort_pathvar$IGF1R <- NA
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE1 %in% c("IGF1R")] <- 1
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE2 %in% c("IGF1R")] <- 1
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE3 %in% c("IGF1R")] <- 1
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE4 %in% c("IGF1R")] <- 1
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE5 %in% c("IGF1R")] <- 1
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE6 %in% c("IGF1R")] <- 1
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE7 %in% c("IGF1R")] <- 1
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE8 %in% c("IGF1R")] <- 1
EntireCohort_pathvar$IGF1R[EntireCohort_pathvar$GENE9 %in% c("IGF1R")] <- 1

EntireCohort_pathvar$IGF1R[is.na(EntireCohort_pathvar$IGF1R)] <- 0
EntireCohort_pathvar$IGF1R <- as.numeric(EntireCohort_pathvar$IGF1R)
unique(EntireCohort_pathvar$IGF1R)
# 0 1

table(EntireCohort_pathvar$IGF1R)
# 0      1 
# 467800   1999




### 'IGFALS'
EntireCohort_pathvar$IGFALS <- NA
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE1 %in% c("IGFALS")] <- 1
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE2 %in% c("IGFALS")] <- 1
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE3 %in% c("IGFALS")] <- 1
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE4 %in% c("IGFALS")] <- 1
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE5 %in% c("IGFALS")] <- 1
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE6 %in% c("IGFALS")] <- 1
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE7 %in% c("IGFALS")] <- 1
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE8 %in% c("IGFALS")] <- 1
EntireCohort_pathvar$IGFALS[EntireCohort_pathvar$GENE9 %in% c("IGFALS")] <- 1

EntireCohort_pathvar$IGFALS[is.na(EntireCohort_pathvar$IGFALS)] <- 0
EntireCohort_pathvar$IGFALS <- as.numeric(EntireCohort_pathvar$IGFALS)
unique(EntireCohort_pathvar$IGFALS)
# 0 1

table(EntireCohort_pathvar$IGFALS)
# 0      1 
# 468516   1283




### 'LHX3'
EntireCohort_pathvar$LHX3 <- NA
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE1 %in% c("LHX3")] <- 1
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE2 %in% c("LHX3")] <- 1
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE3 %in% c("LHX3")] <- 1
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE4 %in% c("LHX3")] <- 1
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE5 %in% c("LHX3")] <- 1
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE6 %in% c("LHX3")] <- 1
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE7 %in% c("LHX3")] <- 1
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE8 %in% c("LHX3")] <- 1
EntireCohort_pathvar$LHX3[EntireCohort_pathvar$GENE9 %in% c("LHX3")] <- 1

EntireCohort_pathvar$LHX3[is.na(EntireCohort_pathvar$LHX3)] <- 0
EntireCohort_pathvar$LHX3 <- as.numeric(EntireCohort_pathvar$LHX3)
unique(EntireCohort_pathvar$LHX3)
# 0 1

table(EntireCohort_pathvar$LHX3)
# 0      1 
# 468886    913




### 'LHX4'
EntireCohort_pathvar$LHX4 <- NA
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE1 %in% c("LHX4")] <- 1
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE2 %in% c("LHX4")] <- 1
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE3 %in% c("LHX4")] <- 1
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE4 %in% c("LHX4")] <- 1
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE5 %in% c("LHX4")] <- 1
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE6 %in% c("LHX4")] <- 1
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE7 %in% c("LHX4")] <- 1
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE8 %in% c("LHX4")] <- 1
EntireCohort_pathvar$LHX4[EntireCohort_pathvar$GENE9 %in% c("LHX4")] <- 1

EntireCohort_pathvar$LHX4[is.na(EntireCohort_pathvar$LHX4)] <- 0
EntireCohort_pathvar$LHX4 <- as.numeric(EntireCohort_pathvar$LHX4)
unique(EntireCohort_pathvar$LHX4)
# 0 1

table(EntireCohort_pathvar$LHX4)
# 0      1 
# 469232    567




### 'OTX2'
EntireCohort_pathvar$OTX2 <- NA
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE1 %in% c("OTX2")] <- 1
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE2 %in% c("OTX2")] <- 1
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE3 %in% c("OTX2")] <- 1
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE4 %in% c("OTX2")] <- 1
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE5 %in% c("OTX2")] <- 1
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE6 %in% c("OTX2")] <- 1
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE7 %in% c("OTX2")] <- 1
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE8 %in% c("OTX2")] <- 1
EntireCohort_pathvar$OTX2[EntireCohort_pathvar$GENE9 %in% c("OTX2")] <- 1

EntireCohort_pathvar$OTX2[is.na(EntireCohort_pathvar$OTX2)] <- 0
EntireCohort_pathvar$OTX2 <- as.numeric(EntireCohort_pathvar$OTX2)
unique(EntireCohort_pathvar$OTX2)
# 0 1

table(EntireCohort_pathvar$OTX2)
# 0      1 
# 469074    725




### 'PNPLA6'
EntireCohort_pathvar$PNPLA6 <- NA
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE1 %in% c("PNPLA6")] <- 1
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE2 %in% c("PNPLA6")] <- 1
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE3 %in% c("PNPLA6")] <- 1
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE4 %in% c("PNPLA6")] <- 1
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE5 %in% c("PNPLA6")] <- 1
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE6 %in% c("PNPLA6")] <- 1
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE7 %in% c("PNPLA6")] <- 1
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE8 %in% c("PNPLA6")] <- 1
EntireCohort_pathvar$PNPLA6[EntireCohort_pathvar$GENE9 %in% c("PNPLA6")] <- 1

EntireCohort_pathvar$PNPLA6[is.na(EntireCohort_pathvar$PNPLA6)] <- 0
EntireCohort_pathvar$PNPLA6 <- as.numeric(EntireCohort_pathvar$PNPLA6)
unique(EntireCohort_pathvar$PNPLA6)
# 0 1

table(EntireCohort_pathvar$PNPLA6)
# 0      1 
# 468124   1675




### 'POU1F1'
EntireCohort_pathvar$POU1F1 <- NA
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE1 %in% c("POU1F1")] <- 1
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE2 %in% c("POU1F1")] <- 1
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE3 %in% c("POU1F1")] <- 1
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE4 %in% c("POU1F1")] <- 1
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE5 %in% c("POU1F1")] <- 1
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE6 %in% c("POU1F1")] <- 1
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE7 %in% c("POU1F1")] <- 1
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE8 %in% c("POU1F1")] <- 1
EntireCohort_pathvar$POU1F1[EntireCohort_pathvar$GENE9 %in% c("POU1F1")] <- 1

EntireCohort_pathvar$POU1F1[is.na(EntireCohort_pathvar$POU1F1)] <- 0
EntireCohort_pathvar$POU1F1 <- as.numeric(EntireCohort_pathvar$POU1F1)
unique(EntireCohort_pathvar$POU1F1)
# 0 1

table(EntireCohort_pathvar$POU1F1)
# 0      1 
# 469210    589




### 'PROP1'
EntireCohort_pathvar$PROP1 <- NA
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE1 %in% c("PROP1")] <- 1
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE2 %in% c("PROP1")] <- 1
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE3 %in% c("PROP1")] <- 1
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE4 %in% c("PROP1")] <- 1
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE5 %in% c("PROP1")] <- 1
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE6 %in% c("PROP1")] <- 1
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE7 %in% c("PROP1")] <- 1
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE8 %in% c("PROP1")] <- 1
EntireCohort_pathvar$PROP1[EntireCohort_pathvar$GENE9 %in% c("PROP1")] <- 1

EntireCohort_pathvar$PROP1[is.na(EntireCohort_pathvar$PROP1)] <- 0
EntireCohort_pathvar$PROP1 <- as.numeric(EntireCohort_pathvar$PROP1)
unique(EntireCohort_pathvar$PROP1)
# 0 1

table(EntireCohort_pathvar$PROP1)
# 0      1 
# 469703     96




### 'STAT5B'
EntireCohort_pathvar$STAT5B <- NA
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE1 %in% c("STAT5B")] <- 1
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE2 %in% c("STAT5B")] <- 1
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE3 %in% c("STAT5B")] <- 1
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE4 %in% c("STAT5B")] <- 1
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE5 %in% c("STAT5B")] <- 1
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE6 %in% c("STAT5B")] <- 1
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE7 %in% c("STAT5B")] <- 1
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE8 %in% c("STAT5B")] <- 1
EntireCohort_pathvar$STAT5B[EntireCohort_pathvar$GENE9 %in% c("STAT5B")] <- 1

EntireCohort_pathvar$STAT5B[is.na(EntireCohort_pathvar$STAT5B)] <- 0
EntireCohort_pathvar$STAT5B <- as.numeric(EntireCohort_pathvar$STAT5B)
unique(EntireCohort_pathvar$STAT5B)
# 0 1

table(EntireCohort_pathvar$STAT5B)
# 0      1 
# 469249    550





# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)

















# GHgene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_GHgene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_GHgene_, "pathvar_GHgene_.csv", row.names=FALSE)
# system('dx upload pathvar_GHgene_.csv --path Emily-folder/results_2/pathvar_GHgene_.csv')


# GHgene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ GHgene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_GHgene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_incidence_1, "pathvar_GHgene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_incidence_1.csv --path Emily-folder/results_2/pathvar_GHgene_incidence_1.csv')


# 'BTK','GH1','GHR','GHRH','GHRHR','GHSR','HESX1','IGF1','IGF1R','IGFALS','LHX3','LHX4','OTX2','PNPLA6','POU1F1','PROP1','STAT5B'
LR_2 <- glm(formula = cancer_ROOC_dummy ~ BTK+GH1+GHR+GHRH+GHRHR+GHSR+HESX1+IGF1+IGF1R+IGFALS+LHX3+LHX4+OTX2+PNPLA6+POU1F1+PROP1+STAT5B, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_GHgene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_incidence_2, "pathvar_GHgene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_incidence_2.csv --path Emily-folder/results_2/pathvar_GHgene_incidence_2.csv')






### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "GHgene", 'BTK', 'GH1', 'GHR', 'GHRH', 'GHRHR', 
                                     'GHSR', 'HESX1', 'IGF1', 'IGF1R', 'IGFALS', 'LHX3', 
                                     'LHX4', 'OTX2', 'PNPLA6', 'POU1F1', 'PROP1', 'STAT5B')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     24

with(data_setA, table(STAT5B))

# GHgene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ GHgene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_GHgene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_incidence_3, "pathvar_GHgene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_incidence_3.csv --path Emily-folder/results_2/pathvar_GHgene_incidence_3.csv')


# age + ethnicity + smoking + sex + GHgene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_GHgene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_incidence_4, "pathvar_GHgene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_incidence_4.csv --path Emily-folder/results_2/pathvar_GHgene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ BTK+GH1+GHR+GHRH+GHRHR+GHSR+HESX1+IGF1+IGF1R+IGFALS+LHX3+LHX4+OTX2+PNPLA6+POU1F1+PROP1+STAT5B, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_GHgene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_incidence_5, "pathvar_GHgene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_incidence_5.csv --path Emily-folder/results_2/pathvar_GHgene_incidence_5.csv')


# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BTK+GH1+GHR+GHRH+GHRHR+GHSR+HESX1+IGF1+IGF1R+IGFALS+LHX3+LHX4+OTX2+PNPLA6+POU1F1+PROP1+STAT5B, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_GHgene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_incidence_6, "pathvar_GHgene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_incidence_6.csv --path Emily-folder/results_2/pathvar_GHgene_incidence_6.csv')




### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "GHgene", 'BTK', 'GH1', 'GHR', 'GHRH', 'GHRHR', 
                                     'GHSR', 'HESX1', 'IGF1', 'IGF1R', 'IGFALS', 'LHX3', 
                                     'LHX4', 'OTX2', 'PNPLA6', 'POU1F1', 'PROP1', 'STAT5B')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     29

with(data_setB, table(STAT5B))


# GHgene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ GHgene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_GHgene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_incidence_7, "pathvar_GHgene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_incidence_7.csv --path Emily-folder/results_2/pathvar_GHgene_incidence_7.csv')


# all variables + GHgene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+GHgene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_GHgene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_incidence_8, "pathvar_GHgene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_incidence_8.csv --path Emily-folder/results_2/pathvar_GHgene_incidence_8.csv')


# GOIs alone
LR_9 <- glm(formula = cancer_ROOC_dummy ~ BTK+GH1+GHR+GHRH+GHRHR+GHSR+HESX1+IGF1+IGF1R+IGFALS+LHX3+LHX4+OTX2+PNPLA6+POU1F1+PROP1+STAT5B, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_GHgene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_incidence_9, "pathvar_GHgene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_incidence_9.csv --path Emily-folder/results_2/pathvar_GHgene_incidence_9.csv')


# all variables + GOIs
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+BTK+GH1+GHR+GHRH+GHRHR+GHSR+HESX1+IGF1+IGF1R+IGFALS+LHX3+LHX4+OTX2+PNPLA6+POU1F1+PROP1+STAT5B, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_GHgene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_incidence_10, "pathvar_GHgene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_incidence_10.csv --path Emily-folder/results_2/pathvar_GHgene_incidence_10.csv')










### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(GHgene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(GHgene, cancer_ROOC_dummy)))
# p-value = 0.4833

with(EntireCohort_pathvar, table(GHgene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(GHgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.6138881 0.6995354
# sample estimates:
#   odds ratio 
# 0.6550611

with(EntireCohort_pathvar, table(GHgene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(GHgene, Ever_smoked_dummy)))
# p-value = 0.05888

with(EntireCohort_pathvar, table(GHgene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(GHgene, Sex_Male)))
# p-value = 0.568

with(EntireCohort_pathvar, table(GHgene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(GHgene, University_education)))
# p-value = 0.5229

with(EntireCohort_pathvar, table(GHgene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(GHgene, Annual_income))	)
# p-value = 0.4339

with(EntireCohort_pathvar, table(GHgene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(GHgene, Smoker_previous)))
# p-value = 0.05281

with(EntireCohort_pathvar, table(GHgene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(GHgene, Smoker_current))	)
# p-value = 0.07426



### data_setA
with(data_setA, table(GHgene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(GHgene, cancer_ROOC_dummy)))
# p-value = 0.6266

with(data_setA, table(GHgene, Ethnicity_White))			
fisher.test(with(data_setA, table(GHgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.6115765 0.6971068
# sample estimates:
#   odds ratio 
# 0.6526729

with(data_setA, table(GHgene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(GHgene, Ever_smoked_dummy)))
# p-value = 0.06606

with(data_setA, table(GHgene, Sex_Male))			
fisher.test(with(data_setA, table(GHgene, Sex_Male)))
# p-value = 0.5723



### data_setB
with(data_setB, table(GHgene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(GHgene, cancer_ROOC_dummy)))
# p-value = 0.1796

with(data_setB, table(GHgene, Ethnicity_White))			
fisher.test(with(data_setB, table(GHgene, Ethnicity_White)))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5983415 0.7001546
# sample estimates:
#   odds ratio 
# 0.6468427

with(data_setB, table(GHgene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(GHgene, Ever_smoked_dummy)))
# p-value = 0.1145

with(data_setB, table(GHgene, Sex_Male))			
fisher.test(with(data_setB, table(GHgene, Sex_Male)))
# p-value = 0.66

with(data_setB, table(GHgene, University_education))			
fisher.test(with(data_setB, table(GHgene, University_education)))
# p-value = 0.8947

with(data_setB, table(GHgene, Annual_income))			
fisher.test(with(data_setB, table(GHgene, Annual_income)))
# p-value = 0.5427






### can plot the continuous variables: 

# make GHgene into factor
df_to_plot <- EntireCohort_pathvar
df_to_plot$GHgene <- as.factor(df_to_plot$GHgene)

# try a violin plot

# Age_at_recruitment  
df_to_plot %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ GHgene, data = df_to_plot)
# W = 3092815355, p-value = 0.06973

# BMI  
df_to_plot %>% 
  ggplot(mapping = aes(x = GHgene, y = BMI)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = GHgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ GHgene, data = df_to_plot)
# W = 3039595188, p-value = 0.973

# Standing_height 
df_to_plot %>% 
  ggplot(mapping = aes(x = GHgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_to_plot %>% 
  ggplot(mapping = aes(x = GHgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ GHgene, data = df_to_plot)
# W = 3162672108, p-value = 3.727e-14

# make this a bar plot
# Weekly_exercise
df_to_plot %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(GHgene), scales = "free_y") +
  theme_light()
# no differences 



### data_setA

# make GHgene into factor
data_setA$GHgene <- as.factor(data_setA$GHgene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ GHgene, data = data_setA)
# W = 3.041e+09, p-value = 0.07775



### data_setB

# make GHgene into factor
data_setB$GHgene <- as.factor(data_setB$GHgene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ GHgene, data = data_setB)
# W = 2026895522, p-value = 0.008838

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = GHgene, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = GHgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ GHgene, data = data_setB)
# W = 1997114150, p-value = 0.9782

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = GHgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = GHgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ GHgene, data = data_setB)
# W = 2073459878, p-value = 1.447e-11


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(GHgene), scales = "free_y") +
  theme_light()
# no difference 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# GHgene Type ------------------------------------------------------------

# make sure the execute the code in the 'GHgenes' section

## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_GHgene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_type_1, "pathvar_GHgene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_type_1.csv --path Emily-folder/results_2/pathvar_GHgene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_GHgene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_type_2, "pathvar_GHgene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_type_2.csv --path Emily-folder/results_2/pathvar_GHgene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_GHgene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_type_3, "pathvar_GHgene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_type_3.csv --path Emily-folder/results_2/pathvar_GHgene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_GHgene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_type_4, "pathvar_GHgene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_type_4.csv --path Emily-folder/results_2/pathvar_GHgene_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_GHgene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_type_5, "pathvar_GHgene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_type_5.csv --path Emily-folder/results_2/pathvar_GHgene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_GHgene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_type_6, "pathvar_GHgene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_type_6.csv --path Emily-folder/results_2/pathvar_GHgene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_GHgene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_type_7, "pathvar_GHgene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_type_7.csv --path Emily-folder/results_2/pathvar_GHgene_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_GHgene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_type_8, "pathvar_GHgene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_type_8.csv --path Emily-folder/results_2/pathvar_GHgene_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_GHgene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_type_9, "pathvar_GHgene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_type_9.csv --path Emily-folder/results_2/pathvar_GHgene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_GHgene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_type_10, "pathvar_GHgene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_type_10.csv --path Emily-folder/results_2/pathvar_GHgene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_GHgene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_type_11, "pathvar_GHgene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_type_11.csv --path Emily-folder/results_2/pathvar_GHgene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_GHgene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_type_12, "pathvar_GHgene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_type_12.csv --path Emily-folder/results_2/pathvar_GHgene_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_GHgene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_type_13, "pathvar_GHgene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_type_13.csv --path Emily-folder/results_2/pathvar_GHgene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_GHgene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_type_14, "pathvar_GHgene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_type_14.csv --path Emily-folder/results_2/pathvar_GHgene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_GHgene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_GHgene_type_15, "pathvar_GHgene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_type_15.csv --path Emily-folder/results_2/pathvar_GHgene_type_15.csv')






# GHgene Age -------------------------------------------------------------

# make sure the execute the code in the 'GHgenes' section

### linear regression
# "GHgene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ GHgene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_GHgene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_GHgene_age_1, "pathvar_GHgene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_age_1.csv --path Emily-folder/results_2/pathvar_GHgene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ BTK+GH1+GHR+GHRH+GHRHR+GHSR+HESX1+IGF1+IGF1R+IGFALS+LHX3+LHX4+OTX2+PNPLA6+POU1F1+PROP1+STAT5B, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_GHgene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_GHgene_age_2, "pathvar_GHgene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_age_2.csv --path Emily-folder/results_2/pathvar_GHgene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "GHgene", 'BTK', 'GH1', 'GHR', 'GHRH', 'GHRHR', 'GHSR', 
                                     'HESX1', 'IGF1', 'IGF1R', 'IGFALS', 'LHX3', 'LHX4', 
                                     'OTX2', 'PNPLA6', 'POU1F1', 'PROP1', 'STAT5B')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     24




data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "GHgene", 'BTK', 'GH1', 'GHR', 'GHRH', 
                                     'GHRHR', 'GHSR', 'HESX1', 'IGF1', 'IGF1R', 'IGFALS', 'LHX3', 
                                     'LHX4', 'OTX2', 'PNPLA6', 'POU1F1', 'PROP1', 'STAT5B')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    29


### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ GHgene, data = data_setD)
# summary(LM_15)
pathvar_GHgene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_GHgene_age_3, "pathvar_GHgene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_age_3.csv --path Emily-folder/results_2/pathvar_GHgene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+GHgene, data = data_setD)
# summary(LM_16)
pathvar_GHgene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_GHgene_age_4, "pathvar_GHgene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_age_4.csv --path Emily-folder/results_2/pathvar_GHgene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ BTK+GH1+GHR+GHRH+GHRHR+GHSR+HESX1+IGF1+IGF1R+IGFALS+LHX3+LHX4+OTX2+PNPLA6+POU1F1+PROP1+STAT5B, data = data_setD)
# summary(LM_17)
pathvar_GHgene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_GHgene_age_5, "pathvar_GHgene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_age_5.csv --path Emily-folder/results_2/pathvar_GHgene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BTK+GH1+GHR+GHRH+GHRHR+GHSR+HESX1+IGF1+IGF1R+IGFALS+LHX3+LHX4+OTX2+PNPLA6+POU1F1+PROP1+STAT5B, data = data_setD)
# summary(LM_18)
pathvar_GHgene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_GHgene_age_6, "pathvar_GHgene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_age_6.csv --path Emily-folder/results_2/pathvar_GHgene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ GHgene, data = data_setE)
# summary(LM_15)
pathvar_GHgene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_GHgene_age_7, "pathvar_GHgene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_age_7.csv --path Emily-folder/results_2/pathvar_GHgene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+GHgene, data = data_setE)
# summary(LM_16)
pathvar_GHgene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_GHgene_age_8, "pathvar_GHgene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_age_8.csv --path Emily-folder/results_2/pathvar_GHgene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ BTK+GH1+GHR+GHRH+GHRHR+GHSR+HESX1+IGF1+IGF1R+IGFALS+LHX3+LHX4+OTX2+PNPLA6+POU1F1+PROP1+STAT5B, data = data_setE)
# summary(LM_17)
pathvar_GHgene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_GHgene_age_9, "pathvar_GHgene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_age_9.csv --path Emily-folder/results_2/pathvar_GHgene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+BTK+GH1+GHR+GHRH+GHRHR+GHSR+HESX1+IGF1+IGF1R+IGFALS+LHX3+LHX4+OTX2+PNPLA6+POU1F1+PROP1+STAT5B, data = data_setE)
# summary(LM_18)
pathvar_GHgene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_GHgene_age_10, "pathvar_GHgene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_GHgene_age_10.csv --path Emily-folder/results_2/pathvar_GHgene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(df_cancer_diagnosis, table(GHgene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(GHgene, Ethnicity_White)))
# p-value = 3.928e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5122663 0.7368038
# sample estimates:
#   odds ratio 
# 0.6121719

with(df_cancer_diagnosis, table(GHgene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(GHgene, Ever_smoked_dummy)))
# p-value = 0.498

with(df_cancer_diagnosis, table(GHgene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(GHgene, Sex_Male)))
# p-value = 0.2301

with(df_cancer_diagnosis, table(GHgene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(GHgene, University_education)))
# p-value = 0.3707

with(df_cancer_diagnosis, table(GHgene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(GHgene, Annual_income)))
# p-value = 0.8829


### data_setD
with(data_setD, table(GHgene, Ethnicity_White))			
fisher.test(with(data_setD, table(GHgene, Ethnicity_White)))
# p-value = 2.593e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5077212 0.7303704
# sample estimates:
#   odds ratio 
# 0.6067999 

with(data_setD, table(GHgene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(GHgene, Ever_smoked_dummy)))
# p-value = 0.4625

with(data_setD, table(GHgene, Sex_Male))			
fisher.test(with(data_setD, table(GHgene, Sex_Male)))
# p-value = 0.1715


### data_setE
with(data_setE, table(GHgene, Ethnicity_White))			
fisher.test(with(data_setE, table(GHgene, Ethnicity_White)))
# p-value = 0.00139
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5386601 0.8628165
# sample estimates:
#   odds ratio 
# 0.6775258

with(data_setE, table(GHgene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(GHgene, Ever_smoked_dummy)))
# p-value = 0.7654

with(data_setE, table(GHgene, Sex_Male))			
fisher.test(with(data_setE, table(GHgene, Sex_Male)))
# p-value = 0.08844

with(data_setE, table(GHgene, University_education))			
fisher.test(with(data_setE, table(GHgene, University_education)))
# p-value = 0.5568

with(data_setE, table(GHgene, Annual_income))			
fisher.test(with(data_setE, table(GHgene, Annual_income)))
# p-value = 0.4999



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 
# make GHgene into factor
df_cancer_diagnosis$GHgene <- as.factor(df_cancer_diagnosis$GHgene)

# try a violin plot


# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ GHgene, data = df_cancer_diagnosis)
# 



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ GHgene, data = df_cancer_diagnosis)
# W = 171908777, p-value = 0.08582



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = GHgene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = GHgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ GHgene, data = df_cancer_diagnosis)
# W = 167842106, p-value = 0.8721



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = GHgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = GHgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ GHgene, data = df_cancer_diagnosis)
# W = 171727371, p-value = 0.03019



# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(GHgene), scales = "free_y") +
  theme_light()
# no difference




### data_setD

data_setD$GHgene <- as.factor(data_setD$GHgene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ GHgene, data = data_setD)
#



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ GHgene, data = data_setD)
# W = 169539835, p-value = 0.08997




### data_setE

data_setE$GHgene <- as.factor(data_setE$GHgene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ GHgene, data = data_setE)
# 



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = GHgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ GHgene, data = data_setE)
# W = 107717452, p-value = 0.04021



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = GHgene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = GHgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ GHgene, data = data_setE)
# W = 106429212, p-value = 0.3102



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = GHgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = GHgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ GHgene, data = data_setE)
# W = 106597196, p-value = 0.2501


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(GHgene), scales = "free_y") +
  theme_light()
# no difference










# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_GHgenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_GHgenes.R')
