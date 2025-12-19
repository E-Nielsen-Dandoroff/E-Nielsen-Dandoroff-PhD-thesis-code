# 09/07/2024

# EntireCohort_pathvar_Biosynthgenes.R

# this script is for the Biosynthesis associated genes group.

# load these in each time:
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# and for saving tidy LR results_2
install.packages('broom')
library(broom)


# Biosynthgene
# INSR
# PIK3R1
# PRMT7
# PEX7
# SMPD4
# PAPSS2
# B3GLCT
# MINPP1
# BPNT2
# TMEM165
# ALDH18A1
# COASY
# GNPAT
# PTDSS1
# MSMO1
# DHCR24
# DHCR7
# PCYT1A
# MTHFS
# NANS
# EBP
# AMPD2
# LBR
# AGPS
# SLC35C1
# PEX5
# SLC39A8
# ASNS
# PSAT1
# DPH1
# EIF2AK3
# EIF2S3
# EIF5A
# RPL13
# PHGDH
# NEPRO













# Biosynthgenes ---------------------------------------------------------------

# load in the dataframe
EntireCohort_pathvar <- fread('/mnt/project/Emily-folder/pathvar/EntireCohort_pathvar.csv', data.table = FALSE)
EntireCohort_pathvar <- EntireCohort_pathvar[-1]

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_pathvar <- left_join(EntireCohort_pathvar, EntireCohort_ready, by = 'eid')

# go thr all the genes to check that they're present
grep("NEPRO", EntireCohort_pathvar$GENE)
# all seem to be present, except for:
# PTDSS1 = 0
# PEX5 = 0



# tag by GH genes
# separate the GENE values into 9 columns and then use these to make a new column
EntireCohort_pathvar <- separate(data = EntireCohort_pathvar, col = GENE, into = c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";")


EntireCohort_pathvar$Biosynthgene <- NA
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE1 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE2 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE3 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE4 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE5 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE6 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE7 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE8 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1
EntireCohort_pathvar$Biosynthgene[EntireCohort_pathvar$GENE9 %in% c('INSR', 'PIK3R1', 'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')] <- 1


EntireCohort_pathvar$Biosynthgene[is.na(EntireCohort_pathvar$Biosynthgene)] <- 0
EntireCohort_pathvar$Biosynthgene <- as.numeric(EntireCohort_pathvar$Biosynthgene)
unique(EntireCohort_pathvar$Biosynthgene)
# 0 1 

table(EntireCohort_pathvar$Biosynthgene)
# 0      1 
# 405291  64508 




# now to make a column for each of the Biosynthgenes
### 'INSR'
EntireCohort_pathvar$INSR <- NA
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE1 %in% c("INSR")] <- 1
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE2 %in% c("INSR")] <- 1
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE3 %in% c("INSR")] <- 1
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE4 %in% c("INSR")] <- 1
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE5 %in% c("INSR")] <- 1
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE6 %in% c("INSR")] <- 1
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE7 %in% c("INSR")] <- 1
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE8 %in% c("INSR")] <- 1
EntireCohort_pathvar$INSR[EntireCohort_pathvar$GENE9 %in% c("INSR")] <- 1

EntireCohort_pathvar$INSR[is.na(EntireCohort_pathvar$INSR)] <- 0
EntireCohort_pathvar$INSR <- as.numeric(EntireCohort_pathvar$INSR)
unique(EntireCohort_pathvar$INSR)
# 0 1

table(EntireCohort_pathvar$INSR)
#  0      1 
# 468326   1473





### 'PIK3R1'
EntireCohort_pathvar$PIK3R1 <- NA
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE1 %in% c("PIK3R1")] <- 1
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE2 %in% c("PIK3R1")] <- 1
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE3 %in% c("PIK3R1")] <- 1
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE4 %in% c("PIK3R1")] <- 1
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE5 %in% c("PIK3R1")] <- 1
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE6 %in% c("PIK3R1")] <- 1
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE7 %in% c("PIK3R1")] <- 1
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE8 %in% c("PIK3R1")] <- 1
EntireCohort_pathvar$PIK3R1[EntireCohort_pathvar$GENE9 %in% c("PIK3R1")] <- 1

EntireCohort_pathvar$PIK3R1[is.na(EntireCohort_pathvar$PIK3R1)] <- 0
EntireCohort_pathvar$PIK3R1 <- as.numeric(EntireCohort_pathvar$PIK3R1)
unique(EntireCohort_pathvar$PIK3R1)
# 0 1

table(EntireCohort_pathvar$PIK3R1)
#  0      1 
# 469410    389




### 'PRMT7'
EntireCohort_pathvar$PRMT7 <- NA
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE1 %in% c("PRMT7")] <- 1
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE2 %in% c("PRMT7")] <- 1
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE3 %in% c("PRMT7")] <- 1
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE4 %in% c("PRMT7")] <- 1
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE5 %in% c("PRMT7")] <- 1
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE6 %in% c("PRMT7")] <- 1
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE7 %in% c("PRMT7")] <- 1
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE8 %in% c("PRMT7")] <- 1
EntireCohort_pathvar$PRMT7[EntireCohort_pathvar$GENE9 %in% c("PRMT7")] <- 1

EntireCohort_pathvar$PRMT7[is.na(EntireCohort_pathvar$PRMT7)] <- 0
EntireCohort_pathvar$PRMT7 <- as.numeric(EntireCohort_pathvar$PRMT7)
unique(EntireCohort_pathvar$PRMT7)
# 0 1

table(EntireCohort_pathvar$PRMT7)
#  0      1 
# 468091   1708




### 'PEX7'
EntireCohort_pathvar$PEX7 <- NA
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE1 %in% c("PEX7")] <- 1
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE2 %in% c("PEX7")] <- 1
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE3 %in% c("PEX7")] <- 1
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE4 %in% c("PEX7")] <- 1
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE5 %in% c("PEX7")] <- 1
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE6 %in% c("PEX7")] <- 1
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE7 %in% c("PEX7")] <- 1
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE8 %in% c("PEX7")] <- 1
EntireCohort_pathvar$PEX7[EntireCohort_pathvar$GENE9 %in% c("PEX7")] <- 1

EntireCohort_pathvar$PEX7[is.na(EntireCohort_pathvar$PEX7)] <- 0
EntireCohort_pathvar$PEX7 <- as.numeric(EntireCohort_pathvar$PEX7)
unique(EntireCohort_pathvar$PEX7)
# 0 1

table(EntireCohort_pathvar$PEX7)
#  0      1 
# 468124   1675 




### 'SMPD4'
EntireCohort_pathvar$SMPD4 <- NA
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE1 %in% c("SMPD4")] <- 1
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE2 %in% c("SMPD4")] <- 1
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE3 %in% c("SMPD4")] <- 1
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE4 %in% c("SMPD4")] <- 1
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE5 %in% c("SMPD4")] <- 1
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE6 %in% c("SMPD4")] <- 1
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE7 %in% c("SMPD4")] <- 1
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE8 %in% c("SMPD4")] <- 1
EntireCohort_pathvar$SMPD4[EntireCohort_pathvar$GENE9 %in% c("SMPD4")] <- 1

EntireCohort_pathvar$SMPD4[is.na(EntireCohort_pathvar$SMPD4)] <- 0
EntireCohort_pathvar$SMPD4 <- as.numeric(EntireCohort_pathvar$SMPD4)
unique(EntireCohort_pathvar$SMPD4)
# 0 1

table(EntireCohort_pathvar$SMPD4)
#  0      1 
# 456726  13073 




### 'PAPSS2'
EntireCohort_pathvar$PAPSS2 <- NA
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE1 %in% c("PAPSS2")] <- 1
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE2 %in% c("PAPSS2")] <- 1
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE3 %in% c("PAPSS2")] <- 1
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE4 %in% c("PAPSS2")] <- 1
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE5 %in% c("PAPSS2")] <- 1
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE6 %in% c("PAPSS2")] <- 1
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE7 %in% c("PAPSS2")] <- 1
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE8 %in% c("PAPSS2")] <- 1
EntireCohort_pathvar$PAPSS2[EntireCohort_pathvar$GENE9 %in% c("PAPSS2")] <- 1

EntireCohort_pathvar$PAPSS2[is.na(EntireCohort_pathvar$PAPSS2)] <- 0
EntireCohort_pathvar$PAPSS2 <- as.numeric(EntireCohort_pathvar$PAPSS2)
unique(EntireCohort_pathvar$PAPSS2)
# 0 1

table(EntireCohort_pathvar$PAPSS2)
#  0      1 
# 467166   2633 




### 'B3GLCT'    
EntireCohort_pathvar$B3GLCT <- NA
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE1 %in% c("B3GLCT")] <- 1
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE2 %in% c("B3GLCT")] <- 1
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE3 %in% c("B3GLCT")] <- 1
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE4 %in% c("B3GLCT")] <- 1
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE5 %in% c("B3GLCT")] <- 1
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE6 %in% c("B3GLCT")] <- 1
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE7 %in% c("B3GLCT")] <- 1
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE8 %in% c("B3GLCT")] <- 1
EntireCohort_pathvar$B3GLCT[EntireCohort_pathvar$GENE9 %in% c("B3GLCT")] <- 1

EntireCohort_pathvar$B3GLCT[is.na(EntireCohort_pathvar$B3GLCT)] <- 0
EntireCohort_pathvar$B3GLCT <- as.numeric(EntireCohort_pathvar$B3GLCT)
unique(EntireCohort_pathvar$B3GLCT)
# 0 1

table(EntireCohort_pathvar$B3GLCT)
#  0      1 
# 467332   2467 




### 'MINPP1'
EntireCohort_pathvar$MINPP1 <- NA
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE1 %in% c("MINPP1")] <- 1
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE2 %in% c("MINPP1")] <- 1
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE3 %in% c("MINPP1")] <- 1
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE4 %in% c("MINPP1")] <- 1
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE5 %in% c("MINPP1")] <- 1
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE6 %in% c("MINPP1")] <- 1
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE7 %in% c("MINPP1")] <- 1
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE8 %in% c("MINPP1")] <- 1
EntireCohort_pathvar$MINPP1[EntireCohort_pathvar$GENE9 %in% c("MINPP1")] <- 1

EntireCohort_pathvar$MINPP1[is.na(EntireCohort_pathvar$MINPP1)] <- 0
EntireCohort_pathvar$MINPP1 <- as.numeric(EntireCohort_pathvar$MINPP1)
unique(EntireCohort_pathvar$MINPP1)
# 0 1

table(EntireCohort_pathvar$MINPP1)
#  0      1 
# 469166    633




### 'BPNT2'
EntireCohort_pathvar$BPNT2 <- NA
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE1 %in% c("BPNT2")] <- 1
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE2 %in% c("BPNT2")] <- 1
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE3 %in% c("BPNT2")] <- 1
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE4 %in% c("BPNT2")] <- 1
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE5 %in% c("BPNT2")] <- 1
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE6 %in% c("BPNT2")] <- 1
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE7 %in% c("BPNT2")] <- 1
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE8 %in% c("BPNT2")] <- 1
EntireCohort_pathvar$BPNT2[EntireCohort_pathvar$GENE9 %in% c("BPNT2")] <- 1

EntireCohort_pathvar$BPNT2[is.na(EntireCohort_pathvar$BPNT2)] <- 0
EntireCohort_pathvar$BPNT2 <- as.numeric(EntireCohort_pathvar$BPNT2)
unique(EntireCohort_pathvar$BPNT2)
# 0 1

table(EntireCohort_pathvar$BPNT2)
#  0      1 
# 469298    501 





### 'TMEM165'
EntireCohort_pathvar$TMEM165 <- NA
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE1 %in% c("TMEM165")] <- 1
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE2 %in% c("TMEM165")] <- 1
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE3 %in% c("TMEM165")] <- 1
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE4 %in% c("TMEM165")] <- 1
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE5 %in% c("TMEM165")] <- 1
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE6 %in% c("TMEM165")] <- 1
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE7 %in% c("TMEM165")] <- 1
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE8 %in% c("TMEM165")] <- 1
EntireCohort_pathvar$TMEM165[EntireCohort_pathvar$GENE9 %in% c("TMEM165")] <- 1

EntireCohort_pathvar$TMEM165[is.na(EntireCohort_pathvar$TMEM165)] <- 0
EntireCohort_pathvar$TMEM165 <- as.numeric(EntireCohort_pathvar$TMEM165)
unique(EntireCohort_pathvar$TMEM165)
# 0 1

table(EntireCohort_pathvar$TMEM165)
#  0      1 
# 469329    470




### 'ALDH18A1'
EntireCohort_pathvar$ALDH18A1 <- NA
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE1 %in% c("ALDH18A1")] <- 1
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE2 %in% c("ALDH18A1")] <- 1
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE3 %in% c("ALDH18A1")] <- 1
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE4 %in% c("ALDH18A1")] <- 1
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE5 %in% c("ALDH18A1")] <- 1
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE6 %in% c("ALDH18A1")] <- 1
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE7 %in% c("ALDH18A1")] <- 1
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE8 %in% c("ALDH18A1")] <- 1
EntireCohort_pathvar$ALDH18A1[EntireCohort_pathvar$GENE9 %in% c("ALDH18A1")] <- 1

EntireCohort_pathvar$ALDH18A1[is.na(EntireCohort_pathvar$ALDH18A1)] <- 0
EntireCohort_pathvar$ALDH18A1 <- as.numeric(EntireCohort_pathvar$ALDH18A1)
unique(EntireCohort_pathvar$ALDH18A1)
# 0 1

table(EntireCohort_pathvar$ALDH18A1)
#  0      1 
# 468481   1318 




### 'COASY'
EntireCohort_pathvar$COASY <- NA
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE1 %in% c("COASY")] <- 1
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE2 %in% c("COASY")] <- 1
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE3 %in% c("COASY")] <- 1
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE4 %in% c("COASY")] <- 1
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE5 %in% c("COASY")] <- 1
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE6 %in% c("COASY")] <- 1
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE7 %in% c("COASY")] <- 1
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE8 %in% c("COASY")] <- 1
EntireCohort_pathvar$COASY[EntireCohort_pathvar$GENE9 %in% c("COASY")] <- 1

EntireCohort_pathvar$COASY[is.na(EntireCohort_pathvar$COASY)] <- 0
EntireCohort_pathvar$COASY <- as.numeric(EntireCohort_pathvar$COASY)
unique(EntireCohort_pathvar$COASY)
# 0 1

table(EntireCohort_pathvar$COASY)
#  0      1 
# 468634   1165 




### 'GNPAT'
EntireCohort_pathvar$GNPAT <- NA
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE1 %in% c("GNPAT")] <- 1
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE2 %in% c("GNPAT")] <- 1
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE3 %in% c("GNPAT")] <- 1
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE4 %in% c("GNPAT")] <- 1
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE5 %in% c("GNPAT")] <- 1
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE6 %in% c("GNPAT")] <- 1
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE7 %in% c("GNPAT")] <- 1
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE8 %in% c("GNPAT")] <- 1
EntireCohort_pathvar$GNPAT[EntireCohort_pathvar$GENE9 %in% c("GNPAT")] <- 1

EntireCohort_pathvar$GNPAT[is.na(EntireCohort_pathvar$GNPAT)] <- 0
EntireCohort_pathvar$GNPAT <- as.numeric(EntireCohort_pathvar$GNPAT)
unique(EntireCohort_pathvar$GNPAT)
# 0 1

table(EntireCohort_pathvar$GNPAT)
#  0      1 
# 468842    957 




### 'PTDSS1'
# no variants present




### 'MSMO1'
EntireCohort_pathvar$MSMO1 <- NA
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE1 %in% c("MSMO1")] <- 1
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE2 %in% c("MSMO1")] <- 1
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE3 %in% c("MSMO1")] <- 1
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE4 %in% c("MSMO1")] <- 1
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE5 %in% c("MSMO1")] <- 1
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE6 %in% c("MSMO1")] <- 1
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE7 %in% c("MSMO1")] <- 1
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE8 %in% c("MSMO1")] <- 1
EntireCohort_pathvar$MSMO1[EntireCohort_pathvar$GENE9 %in% c("MSMO1")] <- 1

EntireCohort_pathvar$MSMO1[is.na(EntireCohort_pathvar$MSMO1)] <- 0
EntireCohort_pathvar$MSMO1 <- as.numeric(EntireCohort_pathvar$MSMO1)
unique(EntireCohort_pathvar$MSMO1)
# 0 1

table(EntireCohort_pathvar$MSMO1)
#  0      1 
# 469177    622 





### 'DHCR24'
EntireCohort_pathvar$DHCR24 <- NA
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE1 %in% c("DHCR24")] <- 1
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE2 %in% c("DHCR24")] <- 1
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE3 %in% c("DHCR24")] <- 1
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE4 %in% c("DHCR24")] <- 1
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE5 %in% c("DHCR24")] <- 1
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE6 %in% c("DHCR24")] <- 1
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE7 %in% c("DHCR24")] <- 1
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE8 %in% c("DHCR24")] <- 1
EntireCohort_pathvar$DHCR24[EntireCohort_pathvar$GENE9 %in% c("DHCR24")] <- 1

EntireCohort_pathvar$DHCR24[is.na(EntireCohort_pathvar$DHCR24)] <- 0
EntireCohort_pathvar$DHCR24 <- as.numeric(EntireCohort_pathvar$DHCR24)
unique(EntireCohort_pathvar$DHCR24)
# 0 1

table(EntireCohort_pathvar$DHCR24)
#  0      1 
# 466322   3477 





### 'DHCR7'
EntireCohort_pathvar$DHCR7 <- NA
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE1 %in% c("DHCR7")] <- 1
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE2 %in% c("DHCR7")] <- 1
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE3 %in% c("DHCR7")] <- 1
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE4 %in% c("DHCR7")] <- 1
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE5 %in% c("DHCR7")] <- 1
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE6 %in% c("DHCR7")] <- 1
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE7 %in% c("DHCR7")] <- 1
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE8 %in% c("DHCR7")] <- 1
EntireCohort_pathvar$DHCR7[EntireCohort_pathvar$GENE9 %in% c("DHCR7")] <- 1

EntireCohort_pathvar$DHCR7[is.na(EntireCohort_pathvar$DHCR7)] <- 0
EntireCohort_pathvar$DHCR7 <- as.numeric(EntireCohort_pathvar$DHCR7)
unique(EntireCohort_pathvar$DHCR7)
# 0 1

table(EntireCohort_pathvar$DHCR7)
#  0      1 
# 457683  12116 




### 'PCYT1A'
EntireCohort_pathvar$PCYT1A <- NA
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE1 %in% c("PCYT1A")] <- 1
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE2 %in% c("PCYT1A")] <- 1
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE3 %in% c("PCYT1A")] <- 1
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE4 %in% c("PCYT1A")] <- 1
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE5 %in% c("PCYT1A")] <- 1
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE6 %in% c("PCYT1A")] <- 1
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE7 %in% c("PCYT1A")] <- 1
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE8 %in% c("PCYT1A")] <- 1
EntireCohort_pathvar$PCYT1A[EntireCohort_pathvar$GENE9 %in% c("PCYT1A")] <- 1

EntireCohort_pathvar$PCYT1A[is.na(EntireCohort_pathvar$PCYT1A)] <- 0
EntireCohort_pathvar$PCYT1A <- as.numeric(EntireCohort_pathvar$PCYT1A)
unique(EntireCohort_pathvar$PCYT1A)
# 0 1

table(EntireCohort_pathvar$PCYT1A)
#  0      1 
# 468548   1251 




### 'MTHFS'
EntireCohort_pathvar$MTHFS <- NA
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE1 %in% c("MTHFS")] <- 1
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE2 %in% c("MTHFS")] <- 1
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE3 %in% c("MTHFS")] <- 1
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE4 %in% c("MTHFS")] <- 1
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE5 %in% c("MTHFS")] <- 1
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE6 %in% c("MTHFS")] <- 1
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE7 %in% c("MTHFS")] <- 1
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE8 %in% c("MTHFS")] <- 1
EntireCohort_pathvar$MTHFS[EntireCohort_pathvar$GENE9 %in% c("MTHFS")] <- 1

EntireCohort_pathvar$MTHFS[is.na(EntireCohort_pathvar$MTHFS)] <- 0
EntireCohort_pathvar$MTHFS <- as.numeric(EntireCohort_pathvar$MTHFS)
unique(EntireCohort_pathvar$MTHFS)
# 0 1

table(EntireCohort_pathvar$MTHFS)
#  0      1 
# 469150    649 




### 'NANS'
EntireCohort_pathvar$NANS <- NA
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE1 %in% c("NANS")] <- 1
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE2 %in% c("NANS")] <- 1
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE3 %in% c("NANS")] <- 1
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE4 %in% c("NANS")] <- 1
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE5 %in% c("NANS")] <- 1
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE6 %in% c("NANS")] <- 1
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE7 %in% c("NANS")] <- 1
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE8 %in% c("NANS")] <- 1
EntireCohort_pathvar$NANS[EntireCohort_pathvar$GENE9 %in% c("NANS")] <- 1

EntireCohort_pathvar$NANS[is.na(EntireCohort_pathvar$NANS)] <- 0
EntireCohort_pathvar$NANS <- as.numeric(EntireCohort_pathvar$NANS)
unique(EntireCohort_pathvar$NANS)
# 0 1

table(EntireCohort_pathvar$NANS)
#  0      1 
# 467588   2211




### 'EBP'
EntireCohort_pathvar$EBP <- NA
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE1 %in% c("EBP")] <- 1
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE2 %in% c("EBP")] <- 1
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE3 %in% c("EBP")] <- 1
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE4 %in% c("EBP")] <- 1
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE5 %in% c("EBP")] <- 1
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE6 %in% c("EBP")] <- 1
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE7 %in% c("EBP")] <- 1
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE8 %in% c("EBP")] <- 1
EntireCohort_pathvar$EBP[EntireCohort_pathvar$GENE9 %in% c("EBP")] <- 1

EntireCohort_pathvar$EBP[is.na(EntireCohort_pathvar$EBP)] <- 0
EntireCohort_pathvar$EBP <- as.numeric(EntireCohort_pathvar$EBP)
unique(EntireCohort_pathvar$EBP)
# 0 1

table(EntireCohort_pathvar$EBP)
#  0      1 
# 469556    243




### 'AMPD2' 
EntireCohort_pathvar$AMPD2 <- NA
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE1 %in% c("AMPD2")] <- 1
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE2 %in% c("AMPD2")] <- 1
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE3 %in% c("AMPD2")] <- 1
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE4 %in% c("AMPD2")] <- 1
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE5 %in% c("AMPD2")] <- 1
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE6 %in% c("AMPD2")] <- 1
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE7 %in% c("AMPD2")] <- 1
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE8 %in% c("AMPD2")] <- 1
EntireCohort_pathvar$AMPD2[EntireCohort_pathvar$GENE9 %in% c("AMPD2")] <- 1

EntireCohort_pathvar$AMPD2[is.na(EntireCohort_pathvar$AMPD2)] <- 0
EntireCohort_pathvar$AMPD2 <- as.numeric(EntireCohort_pathvar$AMPD2)
unique(EntireCohort_pathvar$AMPD2)
# 0 1

table(EntireCohort_pathvar$AMPD2)
#  0      1 
# 468772   1027 






### 'LBR' 
EntireCohort_pathvar$LBR <- NA
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE1 %in% c("LBR")] <- 1
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE2 %in% c("LBR")] <- 1
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE3 %in% c("LBR")] <- 1
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE4 %in% c("LBR")] <- 1
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE5 %in% c("LBR")] <- 1
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE6 %in% c("LBR")] <- 1
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE7 %in% c("LBR")] <- 1
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE8 %in% c("LBR")] <- 1
EntireCohort_pathvar$LBR[EntireCohort_pathvar$GENE9 %in% c("LBR")] <- 1

EntireCohort_pathvar$LBR[is.na(EntireCohort_pathvar$LBR)] <- 0
EntireCohort_pathvar$LBR <- as.numeric(EntireCohort_pathvar$LBR)
unique(EntireCohort_pathvar$LBR)
# 0 1

table(EntireCohort_pathvar$LBR)
#  0      1 
# 467126   2673





### 'AGPS' 
EntireCohort_pathvar$AGPS <- NA
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE1 %in% c("AGPS")] <- 1
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE2 %in% c("AGPS")] <- 1
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE3 %in% c("AGPS")] <- 1
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE4 %in% c("AGPS")] <- 1
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE5 %in% c("AGPS")] <- 1
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE6 %in% c("AGPS")] <- 1
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE7 %in% c("AGPS")] <- 1
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE8 %in% c("AGPS")] <- 1
EntireCohort_pathvar$AGPS[EntireCohort_pathvar$GENE9 %in% c("AGPS")] <- 1

EntireCohort_pathvar$AGPS[is.na(EntireCohort_pathvar$AGPS)] <- 0
EntireCohort_pathvar$AGPS <- as.numeric(EntireCohort_pathvar$AGPS)
unique(EntireCohort_pathvar$AGPS)
# 0 1

table(EntireCohort_pathvar$AGPS)
#  0      1 
# 469338    461 




### 'SLC35C1' 
EntireCohort_pathvar$SLC35C1 <- NA
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE1 %in% c("SLC35C1")] <- 1
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE2 %in% c("SLC35C1")] <- 1
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE3 %in% c("SLC35C1")] <- 1
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE4 %in% c("SLC35C1")] <- 1
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE5 %in% c("SLC35C1")] <- 1
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE6 %in% c("SLC35C1")] <- 1
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE7 %in% c("SLC35C1")] <- 1
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE8 %in% c("SLC35C1")] <- 1
EntireCohort_pathvar$SLC35C1[EntireCohort_pathvar$GENE9 %in% c("SLC35C1")] <- 1

EntireCohort_pathvar$SLC35C1[is.na(EntireCohort_pathvar$SLC35C1)] <- 0
EntireCohort_pathvar$SLC35C1 <- as.numeric(EntireCohort_pathvar$SLC35C1)
unique(EntireCohort_pathvar$SLC35C1)
# 0 1

table(EntireCohort_pathvar$SLC35C1)
#  0      1 
# 469230    569 




### 'PEX5'         
# no variants present




### 'SLC39A8'
EntireCohort_pathvar$SLC39A8 <- NA
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE1 %in% c("SLC39A8")] <- 1
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE2 %in% c("SLC39A8")] <- 1
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE3 %in% c("SLC39A8")] <- 1
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE4 %in% c("SLC39A8")] <- 1
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE5 %in% c("SLC39A8")] <- 1
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE6 %in% c("SLC39A8")] <- 1
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE7 %in% c("SLC39A8")] <- 1
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE8 %in% c("SLC39A8")] <- 1
EntireCohort_pathvar$SLC39A8[EntireCohort_pathvar$GENE9 %in% c("SLC39A8")] <- 1

EntireCohort_pathvar$SLC39A8[is.na(EntireCohort_pathvar$SLC39A8)] <- 0
EntireCohort_pathvar$SLC39A8 <- as.numeric(EntireCohort_pathvar$SLC39A8)
unique(EntireCohort_pathvar$SLC39A8)
# 0 1

table(EntireCohort_pathvar$SLC39A8)
#  0      1 
# 469387    412




### 'ASNS'  
EntireCohort_pathvar$ASNS <- NA
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE1 %in% c("ASNS")] <- 1
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE2 %in% c("ASNS")] <- 1
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE3 %in% c("ASNS")] <- 1
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE4 %in% c("ASNS")] <- 1
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE5 %in% c("ASNS")] <- 1
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE6 %in% c("ASNS")] <- 1
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE7 %in% c("ASNS")] <- 1
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE8 %in% c("ASNS")] <- 1
EntireCohort_pathvar$ASNS[EntireCohort_pathvar$GENE9 %in% c("ASNS")] <- 1

EntireCohort_pathvar$ASNS[is.na(EntireCohort_pathvar$ASNS)] <- 0
EntireCohort_pathvar$ASNS <- as.numeric(EntireCohort_pathvar$ASNS)
unique(EntireCohort_pathvar$ASNS)
# 0 1

table(EntireCohort_pathvar$ASNS)
#  0      1 
# 469052    747




### 'PSAT1'  
EntireCohort_pathvar$PSAT1 <- NA
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE1 %in% c("PSAT1")] <- 1
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE2 %in% c("PSAT1")] <- 1
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE3 %in% c("PSAT1")] <- 1
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE4 %in% c("PSAT1")] <- 1
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE5 %in% c("PSAT1")] <- 1
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE6 %in% c("PSAT1")] <- 1
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE7 %in% c("PSAT1")] <- 1
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE8 %in% c("PSAT1")] <- 1
EntireCohort_pathvar$PSAT1[EntireCohort_pathvar$GENE9 %in% c("PSAT1")] <- 1

EntireCohort_pathvar$PSAT1[is.na(EntireCohort_pathvar$PSAT1)] <- 0
EntireCohort_pathvar$PSAT1 <- as.numeric(EntireCohort_pathvar$PSAT1)
unique(EntireCohort_pathvar$PSAT1)
# 0 1

table(EntireCohort_pathvar$PSAT1)
#   0      1 
# 468749   1050





### 'DPH1'
EntireCohort_pathvar$DPH1 <- NA
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE1 %in% c("DPH1")] <- 1
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE2 %in% c("DPH1")] <- 1
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE3 %in% c("DPH1")] <- 1
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE4 %in% c("DPH1")] <- 1
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE5 %in% c("DPH1")] <- 1
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE6 %in% c("DPH1")] <- 1
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE7 %in% c("DPH1")] <- 1
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE8 %in% c("DPH1")] <- 1
EntireCohort_pathvar$DPH1[EntireCohort_pathvar$GENE9 %in% c("DPH1")] <- 1

EntireCohort_pathvar$DPH1[is.na(EntireCohort_pathvar$DPH1)] <- 0
EntireCohort_pathvar$DPH1 <- as.numeric(EntireCohort_pathvar$DPH1)
unique(EntireCohort_pathvar$DPH1)
# 0 1

table(EntireCohort_pathvar$DPH1)
#  0      1 
# 463052   6747





### 'EIF2AK3'
EntireCohort_pathvar$EIF2AK3 <- NA
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE1 %in% c("EIF2AK3")] <- 1
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE2 %in% c("EIF2AK3")] <- 1
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE3 %in% c("EIF2AK3")] <- 1
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE4 %in% c("EIF2AK3")] <- 1
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE5 %in% c("EIF2AK3")] <- 1
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE6 %in% c("EIF2AK3")] <- 1
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE7 %in% c("EIF2AK3")] <- 1
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE8 %in% c("EIF2AK3")] <- 1
EntireCohort_pathvar$EIF2AK3[EntireCohort_pathvar$GENE9 %in% c("EIF2AK3")] <- 1

EntireCohort_pathvar$EIF2AK3[is.na(EntireCohort_pathvar$EIF2AK3)] <- 0
EntireCohort_pathvar$EIF2AK3 <- as.numeric(EntireCohort_pathvar$EIF2AK3)
unique(EntireCohort_pathvar$EIF2AK3)
# 0 1

table(EntireCohort_pathvar$EIF2AK3)
#  0      1 
# 468939    860




### 'EIF2S3'  
EntireCohort_pathvar$EIF2S3 <- NA
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE1 %in% c("EIF2S3")] <- 1
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE2 %in% c("EIF2S3")] <- 1
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE3 %in% c("EIF2S3")] <- 1
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE4 %in% c("EIF2S3")] <- 1
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE5 %in% c("EIF2S3")] <- 1
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE6 %in% c("EIF2S3")] <- 1
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE7 %in% c("EIF2S3")] <- 1
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE8 %in% c("EIF2S3")] <- 1
EntireCohort_pathvar$EIF2S3[EntireCohort_pathvar$GENE9 %in% c("EIF2S3")] <- 1

EntireCohort_pathvar$EIF2S3[is.na(EntireCohort_pathvar$EIF2S3)] <- 0
EntireCohort_pathvar$EIF2S3 <- as.numeric(EntireCohort_pathvar$EIF2S3)
unique(EntireCohort_pathvar$EIF2S3)
# 0 1

table(EntireCohort_pathvar$EIF2S3)
#  0      1 
# 469784     15 




### 'EIF5A'         
EntireCohort_pathvar$EIF5A <- NA
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE1 %in% c("EIF5A")] <- 1
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE2 %in% c("EIF5A")] <- 1
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE3 %in% c("EIF5A")] <- 1
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE4 %in% c("EIF5A")] <- 1
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE5 %in% c("EIF5A")] <- 1
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE6 %in% c("EIF5A")] <- 1
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE7 %in% c("EIF5A")] <- 1
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE8 %in% c("EIF5A")] <- 1
EntireCohort_pathvar$EIF5A[EntireCohort_pathvar$GENE9 %in% c("EIF5A")] <- 1

EntireCohort_pathvar$EIF5A[is.na(EntireCohort_pathvar$EIF5A)] <- 0
EntireCohort_pathvar$EIF5A <- as.numeric(EntireCohort_pathvar$EIF5A)
unique(EntireCohort_pathvar$EIF5A)
# 0 1

table(EntireCohort_pathvar$EIF5A)
#  0      1 
# 469760     39 




### 'RPL13'
EntireCohort_pathvar$RPL13 <- NA
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE1 %in% c("RPL13")] <- 1
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE2 %in% c("RPL13")] <- 1
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE3 %in% c("RPL13")] <- 1
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE4 %in% c("RPL13")] <- 1
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE5 %in% c("RPL13")] <- 1
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE6 %in% c("RPL13")] <- 1
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE7 %in% c("RPL13")] <- 1
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE8 %in% c("RPL13")] <- 1
EntireCohort_pathvar$RPL13[EntireCohort_pathvar$GENE9 %in% c("RPL13")] <- 1

EntireCohort_pathvar$RPL13[is.na(EntireCohort_pathvar$RPL13)] <- 0
EntireCohort_pathvar$RPL13 <- as.numeric(EntireCohort_pathvar$RPL13)
unique(EntireCohort_pathvar$RPL13)
# 0 1

table(EntireCohort_pathvar$RPL13)
#   0      1 
# 468905    894




### 'PHGDH'
EntireCohort_pathvar$PHGDH <- NA
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE1 %in% c("PHGDH")] <- 1
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE2 %in% c("PHGDH")] <- 1
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE3 %in% c("PHGDH")] <- 1
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE4 %in% c("PHGDH")] <- 1
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE5 %in% c("PHGDH")] <- 1
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE6 %in% c("PHGDH")] <- 1
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE7 %in% c("PHGDH")] <- 1
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE8 %in% c("PHGDH")] <- 1
EntireCohort_pathvar$PHGDH[EntireCohort_pathvar$GENE9 %in% c("PHGDH")] <- 1

EntireCohort_pathvar$PHGDH[is.na(EntireCohort_pathvar$PHGDH)] <- 0
EntireCohort_pathvar$PHGDH <- as.numeric(EntireCohort_pathvar$PHGDH)
unique(EntireCohort_pathvar$PHGDH)
# 0 1

table(EntireCohort_pathvar$PHGDH)
#  0      1 
# 466675   3124




### 'NEPRO'   
EntireCohort_pathvar$NEPRO <- NA
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE1 %in% c("NEPRO")] <- 1
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE2 %in% c("NEPRO")] <- 1
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE3 %in% c("NEPRO")] <- 1
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE4 %in% c("NEPRO")] <- 1
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE5 %in% c("NEPRO")] <- 1
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE6 %in% c("NEPRO")] <- 1
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE7 %in% c("NEPRO")] <- 1
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE8 %in% c("NEPRO")] <- 1
EntireCohort_pathvar$NEPRO[EntireCohort_pathvar$GENE9 %in% c("NEPRO")] <- 1

EntireCohort_pathvar$NEPRO[is.na(EntireCohort_pathvar$NEPRO)] <- 0
EntireCohort_pathvar$NEPRO <- as.numeric(EntireCohort_pathvar$NEPRO)
unique(EntireCohort_pathvar$NEPRO)
# 0 1

table(EntireCohort_pathvar$NEPRO)
#  0      1 
# 468580   1219






# remerge the GENE names columns 
EntireCohort_pathvar <- unite(EntireCohort_pathvar, "GENE", c("GENE1","GENE2","GENE3","GENE4","GENE5","GENE6","GENE7","GENE8","GENE9"), sep = ";", na.rm = TRUE)
unique(EntireCohort_pathvar$GENE)

















# Biosynthgene Incidence ------------------------------------------------------------


# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_pathvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))
#
# # reference code for saving the data as .csv:
# pathvar_Biosynthgene_ <- tidy(LR_, conf.int = TRUE, exp = TRUE)
# write.csv(pathvar_Biosynthgene_, "pathvar_Biosynthgene_.csv", row.names=FALSE)
# system('dx upload pathvar_Biosynthgene_.csv --path Emily-folder/results_2/pathvar_Biosynthgene_.csv')


# Biosynthgene
LR_1 <- glm(formula = cancer_ROOC_dummy ~ Biosynthgene, data = EntireCohort_pathvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
pathvar_Biosynthgene_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_incidence_1, "pathvar_Biosynthgene_incidence_1.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_incidence_1.csv --path Emily-folder/results_2/pathvar_Biosynthgene_incidence_1.csv')


# GOIs
LR_2 <- glm(formula = cancer_ROOC_dummy ~ INSR+PIK3R1+PRMT7+PEX7+SMPD4+PAPSS2+B3GLCT+MINPP1+BPNT2+TMEM165+ALDH18A1+COASY+GNPAT+MSMO1+DHCR24+DHCR7+PCYT1A+MTHFS+NANS+EBP+AMPD2+LBR+AGPS+SLC35C1+SLC39A8+ASNS+PSAT1+DPH1+EIF2AK3+EIF2S3+EIF5A+RPL13+PHGDH+NEPRO, data = EntireCohort_pathvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
pathvar_Biosynthgene_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_incidence_2, "pathvar_Biosynthgene_incidence_2.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_incidence_2.csv --path Emily-folder/results_2/pathvar_Biosynthgene_incidence_2.csv')







### multivariable logistic regression

data_setA <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "Biosynthgene", 'INSR', 'PIK3R1', 'PRMT7', 'PEX7', 
                                     'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 
                                     'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 
                                     'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 
                                     'NANS', 'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 
                                     'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 
                                     'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     41

with(data_setA, table(NEPRO))
# EIF2S3 = 15
# EIF5A = 39




# Biosynthgene by itself
LR_3 <- glm(formula = cancer_ROOC_dummy ~ Biosynthgene, data = data_setA, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
pathvar_Biosynthgene_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_incidence_3, "pathvar_Biosynthgene_incidence_3.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_incidence_3.csv --path Emily-folder/results_2/pathvar_Biosynthgene_incidence_3.csv')


# age + ethnicity + smoking + sex + Biosynthgene
LR_4 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, 
            data = data_setA, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
pathvar_Biosynthgene_incidence_4 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_incidence_4, "pathvar_Biosynthgene_incidence_4.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_incidence_4.csv --path Emily-folder/results_2/pathvar_Biosynthgene_incidence_4.csv')


# GOIs by themselves
LR_5 <- glm(formula = cancer_ROOC_dummy ~ INSR+PIK3R1+PRMT7+PEX7+SMPD4+PAPSS2+B3GLCT+MINPP1+BPNT2+TMEM165+ALDH18A1+COASY+GNPAT+MSMO1+DHCR24+DHCR7+PCYT1A+MTHFS+NANS+EBP+AMPD2+LBR+AGPS+SLC35C1+SLC39A8+ASNS+PSAT1+DPH1+EIF2AK3+EIF2S3+EIF5A+RPL13+PHGDH+NEPRO, 
            data = data_setA, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
pathvar_Biosynthgene_incidence_5 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_incidence_5, "pathvar_Biosynthgene_incidence_5.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_incidence_5.csv --path Emily-folder/results_2/pathvar_Biosynthgene_incidence_5.csv')



# age + ethnicity + smoking + sex + GOIs 
LR_6 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+INSR+PIK3R1+PRMT7+PEX7+SMPD4+PAPSS2+B3GLCT+MINPP1+BPNT2+TMEM165+ALDH18A1+COASY+GNPAT+MSMO1+DHCR24+DHCR7+PCYT1A+MTHFS+NANS+EBP+AMPD2+LBR+AGPS+SLC35C1+SLC39A8+ASNS+PSAT1+DPH1+EIF2AK3+EIF2S3+EIF5A+RPL13+PHGDH+NEPRO, 
            data = data_setA, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
pathvar_Biosynthgene_incidence_6 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_incidence_6, "pathvar_Biosynthgene_incidence_6.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_incidence_6.csv --path Emily-folder/results_2/pathvar_Biosynthgene_incidence_6.csv')











### data_setB
data_setB <- EntireCohort_pathvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                     "BMI", "Standing_height", "University_education",
                                     "Annual_income", "Weekly_exercise",
                                     "Biosynthgene", 'INSR', 'PIK3R1', 'PRMT7', 'PEX7', 
                                     'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 
                                     'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 
                                     'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 
                                     'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1',
                                     'SLC39A8', 'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 
                                     'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 'NEPRO')]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     46

with(data_setB, table(NEPRO))
# EIF2S3 = 11
# EIF5A = 34



# Biosynthgene only
LR_7 <- glm(formula = cancer_ROOC_dummy ~ Biosynthgene, data = data_setB, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
pathvar_Biosynthgene_incidence_7 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_incidence_7, "pathvar_Biosynthgene_incidence_7.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_incidence_7.csv --path Emily-folder/results_2/pathvar_Biosynthgene_incidence_7.csv')


# all variables + Biosynthgene
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+Biosynthgene, 
            data = data_setB, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
pathvar_Biosynthgene_incidence_8 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_incidence_8, "pathvar_Biosynthgene_incidence_8.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_incidence_8.csv --path Emily-folder/results_2/pathvar_Biosynthgene_incidence_8.csv')


# GOIs alone
LR_9 <- glm(formula = cancer_ROOC_dummy ~ INSR+PIK3R1+PRMT7+PEX7+SMPD4+PAPSS2+B3GLCT+MINPP1+BPNT2+TMEM165+ALDH18A1+COASY+GNPAT+MSMO1+DHCR24+DHCR7+PCYT1A+MTHFS+NANS+EBP+AMPD2+LBR+AGPS+SLC35C1+SLC39A8+ASNS+PSAT1+DPH1+EIF2AK3+EIF2S3+EIF5A+RPL13+PHGDH+NEPRO, 
            data = data_setB, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
pathvar_Biosynthgene_incidence_9 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_incidence_9, "pathvar_Biosynthgene_incidence_9.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_incidence_9.csv --path Emily-folder/results_2/pathvar_Biosynthgene_incidence_9.csv')


# all variables + GOIs
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+INSR+PIK3R1+PRMT7+PEX7+SMPD4+PAPSS2+B3GLCT+MINPP1+BPNT2+TMEM165+ALDH18A1+COASY+GNPAT+MSMO1+DHCR24+DHCR7+PCYT1A+MTHFS+NANS+EBP+AMPD2+LBR+AGPS+SLC35C1+SLC39A8+ASNS+PSAT1+DPH1+EIF2AK3+EIF2S3+EIF5A+RPL13+PHGDH+NEPRO, 
             data = data_setB, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
pathvar_Biosynthgene_incidence_10 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_incidence_10, "pathvar_Biosynthgene_incidence_10.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_incidence_10.csv --path Emily-folder/results_2/pathvar_Biosynthgene_incidence_10.csv')










### Fisher test on these datasets

### starting with EntireCohort_pathvar

with(EntireCohort_pathvar, table(Biosynthgene, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_pathvar, table(Biosynthgene, cancer_ROOC_dummy)))
# p-value = 0.9761

with(EntireCohort_pathvar, table(Biosynthgene, Ethnicity_White))
fisher.test(with(EntireCohort_pathvar, table(Biosynthgene, Ethnicity_White)))
# p-value = 2.952e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.043521 1.126818
# sample estimates:
#   odds ratio 
# 1.084273

with(EntireCohort_pathvar, table(Biosynthgene, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_pathvar, table(Biosynthgene, Ever_smoked_dummy)))
# p-value = 0.1564

with(EntireCohort_pathvar, table(Biosynthgene, Sex_Male))			
fisher.test(with(EntireCohort_pathvar, table(Biosynthgene, Sex_Male)))
# p-value = 0.6613

with(EntireCohort_pathvar, table(Biosynthgene, University_education))			
fisher.test(with(EntireCohort_pathvar, table(Biosynthgene, University_education)))
# p-value = 0.01122
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9597956 0.9947736
# sample estimates:
#   odds ratio 
# 0.977145

with(EntireCohort_pathvar, table(Biosynthgene, Annual_income))			
fisher.test(with(EntireCohort_pathvar, table(Biosynthgene, Annual_income))	)
# p-value = 0.0192
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9555497 0.9959916
# sample estimates:
#   odds ratio 
# 0.9755823 

with(EntireCohort_pathvar, table(Biosynthgene, Smoker_previous))			
fisher.test(with(EntireCohort_pathvar, table(Biosynthgene, Smoker_previous)))
# p-value = 0.06268

with(EntireCohort_pathvar, table(Biosynthgene, Smoker_current))			
fisher.test(with(EntireCohort_pathvar, table(Biosynthgene, Smoker_current))	)
# p-value = 0.5708



### data_setA
with(data_setA, table(Biosynthgene, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(Biosynthgene, cancer_ROOC_dummy)))
# p-value = 0.9201

with(data_setA, table(Biosynthgene, Ethnicity_White))			
fisher.test(with(data_setA, table(Biosynthgene, Ethnicity_White)))
# p-value = 3.669e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.042734 1.126307
# sample estimates:
#   odds ratio 
# 1.083616

with(data_setA, table(Biosynthgene, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(Biosynthgene, Ever_smoked_dummy)))
# p-value = 0.147

with(data_setA, table(Biosynthgene, Sex_Male))			
fisher.test(with(data_setA, table(Biosynthgene, Sex_Male)))
# p-value = 0.6445



### data_setB
with(data_setB, table(Biosynthgene, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(Biosynthgene, cancer_ROOC_dummy)))
# p-value = 0.6881

with(data_setB, table(Biosynthgene, Ethnicity_White))			
fisher.test(with(data_setB, table(Biosynthgene, Ethnicity_White)))
# p-value = 0.002892
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.023901 1.123151
# sample estimates:
#   odds ratio 
# 1.072193

with(data_setB, table(Biosynthgene, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(Biosynthgene, Ever_smoked_dummy)))
# p-value = 0.2233

with(data_setB, table(Biosynthgene, Sex_Male))			
fisher.test(with(data_setB, table(Biosynthgene, Sex_Male)))
# p-value = 0.9474

with(data_setB, table(Biosynthgene, University_education))			
fisher.test(with(data_setB, table(Biosynthgene, University_education)))
# p-value = 0.07555

with(data_setB, table(Biosynthgene, Annual_income))			
fisher.test(with(data_setB, table(Biosynthgene, Annual_income)))
# p-value = 0.04605
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9584598 0.9996408
# sample estimates:
#   odds ratio 
# 0.9788583 






### can plot the continuous variables: 

# make Biosynthgene into factor
df_for_plots <- EntireCohort_pathvar
df_for_plots$Biosynthgene <- as.factor(df_for_plots$Biosynthgene)

# try a violin plot

# Age_at_recruitment  
df_for_plots %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ Biosynthgene, data = df_for_plots)
# W = 1.3007e+10, p-value = 0.04147

# BMI  
df_for_plots %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = BMI)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ Biosynthgene, data = df_for_plots)
# W = 1.2912e+10, p-value = 0.07475

# Standing_height 
df_for_plots %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
df_for_plots %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ Biosynthgene, data = df_for_plots)
# W = 1.3011e+10, p-value = 0.6281

# make this a bar plot
# Weekly_exercise
df_for_plots %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(Biosynthgene), scales = "free_y") +
  theme_light()
# 



### data_setA

# make Biosynthgene into factor
data_setA$Biosynthgene <- as.factor(data_setA$Biosynthgene)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ Biosynthgene, data = data_setA)
# W = 1.2789e+10, p-value = 0.03972



### data_setB

# make Biosynthgene into factor
data_setB$Biosynthgene <- as.factor(data_setB$Biosynthgene)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ Biosynthgene, data = data_setB)
# W = 8555388061, p-value = 0.05349

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ Biosynthgene, data = data_setB)
# W = 8577899132, p-value = 0.3337

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ Biosynthgene, data = data_setB)
# W = 8618700516, p-value = 0.4358


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(Biosynthgene), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes















# Biosynthgene Type ------------------------------------------------------------

# make sure the execute the code in the 'Biosynthgenes' section

## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
pathvar_Biosynthgene_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_type_1, "pathvar_Biosynthgene_type_1.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_type_1.csv --path Emily-folder/results_2/pathvar_Biosynthgene_type_1.csv')



LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
pathvar_Biosynthgene_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_type_2, "pathvar_Biosynthgene_type_2.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_type_2.csv --path Emily-folder/results_2/pathvar_Biosynthgene_type_2.csv')



LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
pathvar_Biosynthgene_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_type_3, "pathvar_Biosynthgene_type_3.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_type_3.csv --path Emily-folder/results_2/pathvar_Biosynthgene_type_3.csv')



LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
pathvar_Biosynthgene_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_type_4, "pathvar_Biosynthgene_type_4.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_type_4.csv --path Emily-folder/results_2/pathvar_Biosynthgene_type_4.csv')



LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, 
               data = EntireCohort_pathvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
pathvar_Biosynthgene_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_type_5, "pathvar_Biosynthgene_type_5.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_type_5.csv --path Emily-folder/results_2/pathvar_Biosynthgene_type_5.csv')



LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
pathvar_Biosynthgene_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_type_6, "pathvar_Biosynthgene_type_6.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_type_6.csv --path Emily-folder/results_2/pathvar_Biosynthgene_type_6.csv')



LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, 
                 data = EntireCohort_pathvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
pathvar_Biosynthgene_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_type_7, "pathvar_Biosynthgene_type_7.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_type_7.csv --path Emily-folder/results_2/pathvar_Biosynthgene_type_7.csv')



LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, 
                        data = EntireCohort_pathvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
pathvar_Biosynthgene_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_type_8, "pathvar_Biosynthgene_type_8.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_type_8.csv --path Emily-folder/results_2/pathvar_Biosynthgene_type_8.csv')



LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, 
                      data = EntireCohort_pathvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
pathvar_Biosynthgene_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_type_9, "pathvar_Biosynthgene_type_9.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_type_9.csv --path Emily-folder/results_2/pathvar_Biosynthgene_type_9.csv')



LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, 
                  data = EntireCohort_pathvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
pathvar_Biosynthgene_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_type_10, "pathvar_Biosynthgene_type_10.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_type_10.csv --path Emily-folder/results_2/pathvar_Biosynthgene_type_10.csv')



LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, 
              data = EntireCohort_pathvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
pathvar_Biosynthgene_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_type_11, "pathvar_Biosynthgene_type_11.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_type_11.csv --path Emily-folder/results_2/pathvar_Biosynthgene_type_11.csv')



LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
pathvar_Biosynthgene_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_type_12, "pathvar_Biosynthgene_type_12.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_type_12.csv --path Emily-folder/results_2/pathvar_Biosynthgene_type_12.csv')



LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
pathvar_Biosynthgene_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_type_13, "pathvar_Biosynthgene_type_13.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_type_13.csv --path Emily-folder/results_2/pathvar_Biosynthgene_type_13.csv')



LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, 
                    data = EntireCohort_pathvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
pathvar_Biosynthgene_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_type_14, "pathvar_Biosynthgene_type_14.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_type_14.csv --path Emily-folder/results_2/pathvar_Biosynthgene_type_14.csv')



LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, 
                data = EntireCohort_pathvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
pathvar_Biosynthgene_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(pathvar_Biosynthgene_type_15, "pathvar_Biosynthgene_type_15.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_type_15.csv --path Emily-folder/results_2/pathvar_Biosynthgene_type_15.csv')









# Biosynthgene Age -------------------------------------------------------------

# make sure the execute the code in the 'Biosynthgenes' section

### linear regression
# "Biosynthgene"                                      
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ Biosynthgene, data = EntireCohort_pathvar)
# summary(LM_1)
pathvar_Biosynthgene_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(pathvar_Biosynthgene_age_1, "pathvar_Biosynthgene_age_1.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_age_1.csv --path Emily-folder/results_2/pathvar_Biosynthgene_age_1.csv')


# GOIs                                      
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ INSR+PIK3R1+PRMT7+PEX7+SMPD4+PAPSS2+B3GLCT+MINPP1+BPNT2+TMEM165+ALDH18A1+COASY+GNPAT+MSMO1+DHCR24+DHCR7+PCYT1A+MTHFS+NANS+EBP+AMPD2+LBR+AGPS+SLC35C1+SLC39A8+ASNS+PSAT1+DPH1+EIF2AK3+EIF2S3+EIF5A+RPL13+PHGDH+NEPRO, data = EntireCohort_pathvar)
# summary(LM_2)
pathvar_Biosynthgene_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(pathvar_Biosynthgene_age_2, "pathvar_Biosynthgene_age_2.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_age_2.csv --path Emily-folder/results_2/pathvar_Biosynthgene_age_2.csv')


# haven't done the analysis for each individual gene






## multivariate analysis

# make dataset without NA rows

data_setD <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                     "Biosynthgene", 'INSR', 'PIK3R1', 'PRMT7', 'PEX7', 
                                     'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 'BPNT2', 'TMEM165', 
                                     'ALDH18A1', 'COASY', 'GNPAT', 'MSMO1', 'DHCR24', 
                                     'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 'EBP', 'AMPD2', 'LBR', 
                                     'AGPS', 'SLC35C1', 'SLC39A8', 'ASNS', 'PSAT1', 
                                     'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 'RPL13', 'PHGDH', 
                                     'NEPRO')]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     41

with(data_setD, table(NEPRO))
# PIK3R1 = 88
# EBP = 50
# SLC39A8 = 81
# EIF2S3 = 4
# EIF5A = 7






data_setE <- EntireCohort_pathvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                     "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                     "Standing_height", "University_education", "Annual_income",
                                     "Weekly_exercise", "Biosynthgene", 'INSR', 'PIK3R1', 
                                     'PRMT7', 'PEX7', 'SMPD4', 'PAPSS2', 'B3GLCT', 'MINPP1', 
                                     'BPNT2', 'TMEM165', 'ALDH18A1', 'COASY', 'GNPAT', 
                                     'MSMO1', 'DHCR24', 'DHCR7', 'PCYT1A', 'MTHFS', 'NANS', 
                                     'EBP', 'AMPD2', 'LBR', 'AGPS', 'SLC35C1', 'SLC39A8', 
                                     'ASNS', 'PSAT1', 'DPH1', 'EIF2AK3', 'EIF2S3', 'EIF5A', 
                                     'RPL13', 'PHGDH', 'NEPRO')]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    46

with(data_setE, table(NEPRO))
# PIK3R1 = 72
# BPNT2 = 88
# TMEM165 = 88
# EBP = 37
# AGPS = 88
# SLC39A8 = 63
# EIF2S3 = 3
# EIF5A = 6




### data_setD
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ Biosynthgene, data = data_setD)
# summary(LM_15)
pathvar_Biosynthgene_age_3 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_Biosynthgene_age_3, "pathvar_Biosynthgene_age_3.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_age_3.csv --path Emily-folder/results_2/pathvar_Biosynthgene_age_3.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+Biosynthgene, data = data_setD)
# summary(LM_16)
pathvar_Biosynthgene_age_4 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_Biosynthgene_age_4, "pathvar_Biosynthgene_age_4.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_age_4.csv --path Emily-folder/results_2/pathvar_Biosynthgene_age_4.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ INSR+PIK3R1+PRMT7+PEX7+SMPD4+PAPSS2+B3GLCT+MINPP1+BPNT2+TMEM165+ALDH18A1+COASY+GNPAT+MSMO1+DHCR24+DHCR7+PCYT1A+MTHFS+NANS+EBP+AMPD2+LBR+AGPS+SLC35C1+SLC39A8+ASNS+PSAT1+DPH1+EIF2AK3+EIF2S3+EIF5A+RPL13+PHGDH+NEPRO, data = data_setD)
# summary(LM_17)
pathvar_Biosynthgene_age_5 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_Biosynthgene_age_5, "pathvar_Biosynthgene_age_5.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_age_5.csv --path Emily-folder/results_2/pathvar_Biosynthgene_age_5.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+INSR+PIK3R1+PRMT7+PEX7+SMPD4+PAPSS2+B3GLCT+MINPP1+BPNT2+TMEM165+ALDH18A1+COASY+GNPAT+MSMO1+DHCR24+DHCR7+PCYT1A+MTHFS+NANS+EBP+AMPD2+LBR+AGPS+SLC35C1+SLC39A8+ASNS+PSAT1+DPH1+EIF2AK3+EIF2S3+EIF5A+RPL13+PHGDH+NEPRO, data = data_setD)
# summary(LM_18)
pathvar_Biosynthgene_age_6 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_Biosynthgene_age_6, "pathvar_Biosynthgene_age_6.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_age_6.csv --path Emily-folder/results_2/pathvar_Biosynthgene_age_6.csv')



### data_setE
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ Biosynthgene, data = data_setE)
# summary(LM_15)
pathvar_Biosynthgene_age_7 <- tidy(LM_15, conf.int = TRUE)
write.csv(pathvar_Biosynthgene_age_7, "pathvar_Biosynthgene_age_7.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_age_7.csv --path Emily-folder/results_2/pathvar_Biosynthgene_age_7.csv')


LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+Biosynthgene, data = data_setE)
# summary(LM_16)
pathvar_Biosynthgene_age_8 <- tidy(LM_16, conf.int = TRUE)
write.csv(pathvar_Biosynthgene_age_8, "pathvar_Biosynthgene_age_8.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_age_8.csv --path Emily-folder/results_2/pathvar_Biosynthgene_age_8.csv')


LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ INSR+PIK3R1+PRMT7+PEX7+SMPD4+PAPSS2+B3GLCT+MINPP1+BPNT2+TMEM165+ALDH18A1+COASY+GNPAT+MSMO1+DHCR24+DHCR7+PCYT1A+MTHFS+NANS+EBP+AMPD2+LBR+AGPS+SLC35C1+SLC39A8+ASNS+PSAT1+DPH1+EIF2AK3+EIF2S3+EIF5A+RPL13+PHGDH+NEPRO, data = data_setE)
# summary(LM_17)
pathvar_Biosynthgene_age_9 <- tidy(LM_17, conf.int = TRUE)
write.csv(pathvar_Biosynthgene_age_9, "pathvar_Biosynthgene_age_9.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_age_9.csv --path Emily-folder/results_2/pathvar_Biosynthgene_age_9.csv')


LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+INSR+PIK3R1+PRMT7+PEX7+SMPD4+PAPSS2+B3GLCT+MINPP1+BPNT2+TMEM165+ALDH18A1+COASY+GNPAT+MSMO1+DHCR24+DHCR7+PCYT1A+MTHFS+NANS+EBP+AMPD2+LBR+AGPS+SLC35C1+SLC39A8+ASNS+PSAT1+DPH1+EIF2AK3+EIF2S3+EIF5A+RPL13+PHGDH+NEPRO, data = data_setE)
# summary(LM_18)
pathvar_Biosynthgene_age_10 <- tidy(LM_18, conf.int = TRUE)
write.csv(pathvar_Biosynthgene_age_10, "pathvar_Biosynthgene_age_10.csv", row.names=FALSE)
system('dx upload pathvar_Biosynthgene_age_10.csv --path Emily-folder/results_2/pathvar_Biosynthgene_age_10.csv')




## fisher test for the linear regression data

# need to first subset by cancer diagnosis
df_cancer_diagnosis <- EntireCohort_pathvar
df_cancer_diagnosis <- df_cancer_diagnosis[complete.cases(df_cancer_diagnosis$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(df_cancer_diagnosis, table(Biosynthgene, Ethnicity_White))			
fisher.test(with(df_cancer_diagnosis, table(Biosynthgene, Ethnicity_White)))
# p-value = 0.01199
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.030360 1.291366
# sample estimates:
#   odds ratio 
# 1.152153

with(df_cancer_diagnosis, table(Biosynthgene, Ever_smoked_dummy))			
fisher.test(with(df_cancer_diagnosis, table(Biosynthgene, Ever_smoked_dummy)))
# p-value = 0.2073

with(df_cancer_diagnosis, table(Biosynthgene, Sex_Male))			
fisher.test(with(df_cancer_diagnosis, table(Biosynthgene, Sex_Male)))
# p-value = 0.6363

with(df_cancer_diagnosis, table(Biosynthgene, University_education))			
fisher.test(with(df_cancer_diagnosis, table(Biosynthgene, University_education)))
# p-value = 0.03362
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.9240616 0.9968907
# sample estimates:
#   odds ratio 
# 0.9598261

with(df_cancer_diagnosis, table(Biosynthgene, Annual_income))			
fisher.test(with(df_cancer_diagnosis, table(Biosynthgene, Annual_income)))
# p-value = 0.07622


### data_setD
with(data_setD, table(Biosynthgene, Ethnicity_White))			
fisher.test(with(data_setD, table(Biosynthgene, Ethnicity_White)))
# p-value = 0.009916
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.035247 1.299439
# sample estimates:
#   odds ratio 
# 1.158469 

with(data_setD, table(Biosynthgene, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(Biosynthgene, Ever_smoked_dummy)))
# p-value = 0.2269

with(data_setD, table(Biosynthgene, Sex_Male))			
fisher.test(with(data_setD, table(Biosynthgene, Sex_Male)))
# p-value = 0.7053


### data_setE
with(data_setE, table(Biosynthgene, Ethnicity_White))			
fisher.test(with(data_setE, table(Biosynthgene, Ethnicity_White)))
# p-value = 0.2753

with(data_setE, table(Biosynthgene, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(Biosynthgene, Ever_smoked_dummy)))
# p-value = 0.3092

with(data_setE, table(Biosynthgene, Sex_Male))			
fisher.test(with(data_setE, table(Biosynthgene, Sex_Male)))
# p-value = 0.6245

with(data_setE, table(Biosynthgene, University_education))			
fisher.test(with(data_setE, table(Biosynthgene, University_education)))
# p-value = 0.07934

with(data_setE, table(Biosynthgene, Annual_income))			
fisher.test(with(data_setE, table(Biosynthgene, Annual_income)))
# p-value = 0.09746



## make plots for the linear regression data

# starting with the EntireCohort_pathvar dataset subset by age at cancer diagnosis 

# make Biosynthgene into factor
df_cancer_diagnosis$Biosynthgene <- as.factor(df_cancer_diagnosis$Biosynthgene)

# try a violin plot



# Age_at_cancer_diagnosis_earliest
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ Biosynthgene, data = df_cancer_diagnosis)
# 



# Age_at_recruitment  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ Biosynthgene, data = df_cancer_diagnosis)
# W = 725444162, p-value = 0.4531



# BMI  
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = BMI)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ Biosynthgene, data = df_cancer_diagnosis)
# W = 716183726, p-value = 0.06335



# Standing_height 
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ Biosynthgene, data = df_cancer_diagnosis)
# W = 722311684, p-value = 0.5078



# make this a bar plot
# Weekly_exercise
df_cancer_diagnosis %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(Biosynthgene), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$Biosynthgene <- as.factor(data_setD$Biosynthgene)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ Biosynthgene, data = data_setD)
# 



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ Biosynthgene, data = data_setD)
# W = 713813388, p-value = 0.5146




### data_setE

data_setE$Biosynthgene <- as.factor(data_setE$Biosynthgene)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ Biosynthgene, data = data_setE)
# 



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ Biosynthgene, data = data_setE)
# W = 459215432, p-value = 0.2577



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ Biosynthgene, data = data_setE)
# W = 456424957, p-value = 0.02784



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = Biosynthgene, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ Biosynthgene, data = data_setE)
# W = 459949744, p-value = 0.3955


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(Biosynthgene), scales = "free_y") +
  theme_light()
# 













# save the script ---------------------------------------------------------


system('dx upload EntireCohort_pathvar_Biosynthgenes.R --path Emily-folder/pathvar/EntireCohort_pathvar_Biosynthgenes.R')
