## Cutpoint analysis for in silico predictive tools


# execute each time -------------------------------------------------------


library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggpubr)


# execute once ------------------------------------------------------------

# removing duplicate rows

Extra_variants_data <- read_csv("Extra_variants_2.csv")
View(Extra_variants_data)

# need to change the 'Type' data to factor - to allow grouping based on Pathogenic or Benign
Extra_variants_data$Type <- as.factor(Extra_variants_data$Type)


# have discovered a way to remove duplicate rows!
# note: have done this all once 

unique_data <- distinct(Extra_variants_data, Type, chr, pos, ref, alt, rs_dbSNP, genename,
                        M_CAP_score, M_CAP_rankscore, REVEL_score, REVEL_rankscore,
                        BayesDel_addAF_score, BayesDel_noAF_score, CADD_raw, CADD_phred,
                        DANN_score, fathmm_MKL_coding_score, Eigen_raw_coding,
                        Eigen_phred_coding, GERP_NR, GERP_RS, phyloP100way_vertebrate,
                        phyloP470way_mammalian, phyloP17way_primate, phastCons100way_vertebrate,
                        phastCons470way_mammalian, phastCons17way_primate, SiPhy_29way_logOdds,
                        bStatistic)

View(unique_data)

# save the dataframe
write.csv(unique_data, "Extra_variants_2_unique.csv")
# then edit manually to make sure there aren't still double ups with missing M_CAP/REVEL values!


# # remove GMNN! not relevant for hypomorphic/LoF
# Variant_data <- Variant_data[!(Variant_data$genename == "GMNN"),]
# write.csv(Variant_data, "Extra_variants_3_unique.csv")


# score threshold analysis -----------------------------------------------------

Variant_data <- read_csv("Extra_variants_3_unique.csv")
View(Variant_data)

# need to change the 'Type' data to factor - to allow grouping based on Pathogenic or Benign
Variant_data$Type <- as.factor(Variant_data$Type)
Variant_data$Syndrome <- as.factor(Variant_data$Syndrome)




library(OptimalCutpoints)
library(cutpointr)


#################### not sure if this is necessary 
# subsetting out MGS and Benign
Pathogenic_only <- Variant_data %>% subset(Type == "Pathogenic")
# n = 129
Benign_only <- Variant_data %>% subset(Type == "Benign")
# n = 134






# OptimalCutpoints --------------------------------------------------------
# Using OptimalCutpoints package to generate new thresholds

# Optimal Cutpoints needs data to be in data.frame format:
Variant_data_copy <- as.data.frame(Variant_data)


# the SpEqualSe cutpoint added as comment under each block of code:
OC_BayesDel_AF <- Variant_data_copy %>% 
  optimal.cutpoints(X = "BayesDel_addAF_score", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# -0.0292

OC_BayesDel_noAF <- Variant_data_copy %>% 
  optimal.cutpoints(X = "BayesDel_noAF_score", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# -0.1086

OC_bStatistic <- Variant_data_copy %>% 
  optimal.cutpoints(X = "bStatistic", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 793.0000

OC_CADD_raw <- Variant_data_copy %>% 
  optimal.cutpoints(X = "CADD_raw", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 3.4618

OC_CADD_phred <- Variant_data_copy %>% 
  optimal.cutpoints(X = "CADD_phred", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 24.6000

OC_DANN <- Variant_data_copy %>% 
  optimal.cutpoints(X = "DANN_score", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 0.9933

OC_Eigen_raw <- Variant_data_copy %>% 
  optimal.cutpoints(X = "Eigen_raw_coding", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 0.4058

OC_Eigen_phred <- Variant_data_copy %>% 
  optimal.cutpoints(X = "Eigen_phred_coding", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 4.3792

OC_fathmm_MKL <- Variant_data_copy %>% 
  optimal.cutpoints(X = "fathmm_MKL_coding_score", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 0.8794

OC_GERP_RS <- Variant_data_copy %>% 
  optimal.cutpoints(X = "GERP_RS", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 4.2800

OC_GERP_NR <- Variant_data_copy %>% 
  optimal.cutpoints(X = "GERP_NR", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 5.2500

OC_M_CAP_score <- Variant_data_copy %>% 
  optimal.cutpoints(X = "M_CAP_score", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 0.0505

OC_M_CAP_rankscore <- Variant_data_copy %>% 
  optimal.cutpoints(X = "M_CAP_rankscore", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 0.6425

OC_phastCons_17P <- Variant_data_copy %>% 
  optimal.cutpoints(X = "phastCons17way_primate", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 0.9210

OC_phastCons_470M <- Variant_data_copy %>% 
  optimal.cutpoints(X = "phastCons470way_mammalian", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 1.0000

OC_phastCons_100V <- Variant_data_copy %>% 
  optimal.cutpoints(X = "phastCons100way_vertebrate", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 1.0000

OC_phyloP_17P <- Variant_data_copy %>% 
  optimal.cutpoints(X = "phyloP17way_primate", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 0.6090

OC_phyloP_470M <- Variant_data_copy %>% 
  optimal.cutpoints(X = "phyloP470way_mammalian", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 5.8610

OC_phyloP_100V <- Variant_data_copy %>% 
  optimal.cutpoints(X = "phyloP100way_vertebrate", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 4.1470

OC_REVEL_score <- Variant_data_copy %>% 
  optimal.cutpoints(X = "REVEL_score", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 0.3140

OC_REVEL_rankscore <- Variant_data_copy %>% 
  optimal.cutpoints(X = "REVEL_rankscore", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 0.6359

OC_SiPhy <- Variant_data_copy %>% 
  optimal.cutpoints(X = "SiPhy_29way_logOdds", status = "Type", tag.healthy = "Benign", methods = "SpEqualSe")
# 11.9800








# cutpointr ---------------------------------------------------------------
# Using cutpointr to generate cutpoints

CR_BayesDel_AF <- Variant_data_copy %>% cutpointr(x = BayesDel_addAF_score, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# -0.125064

CR_BayesDel_noAF <- Variant_data_copy %>% cutpointr(x = BayesDel_noAF_score, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# -0.125788

CR_bStatistic <- Variant_data_copy %>% cutpointr(x = bStatistic, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 868

CR_CADD_raw <- Variant_data_copy %>% cutpointr(x = CADD_raw, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 3.26586

CR_CADD_phred <- Variant_data_copy %>% cutpointr(x = CADD_phred, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 24

CR_DANN <- Variant_data_copy %>% cutpointr(x = DANN_score, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 0.992171

CR_Eigen_raw <- Variant_data_copy %>% cutpointr(x = Eigen_raw_coding, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 0.435474

CR_Eigen_phred <- Variant_data_copy %>% cutpointr(x = Eigen_phred_coding, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 4.57111

CR_fathmm_MKL <- Variant_data_copy %>% cutpointr(x = fathmm_MKL_coding_score, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 0.71135

CR_GERP_RS <- Variant_data_copy %>% cutpointr(x = GERP_RS, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 2.51

CR_GERP_NR <- Variant_data_copy %>% cutpointr(x = GERP_NR, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 4.9 

CR_M_CAP_score <- Variant_data_copy %>% cutpointr(x = M_CAP_score, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 0.036617 

CR_M_CAP_rankscore <- Variant_data_copy %>% cutpointr(x = M_CAP_rankscore, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 0.57104 

CR_phastCons_17P <- Variant_data_copy %>% cutpointr(x = phastCons17way_primate, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 0.941 

CR_phastCons_470M <- Variant_data_copy %>% cutpointr(x = phastCons470way_mammalian, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 1

CR_phastCons_100V <- Variant_data_copy %>% cutpointr(x = phastCons100way_vertebrate, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 0.994 

CR_phyloP_17P <- Variant_data_copy %>% cutpointr(x = phyloP17way_primate, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 0.323 

CR_phyloP_470M <- Variant_data_copy %>% cutpointr(x = phyloP470way_mammalian, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 5.719 

CR_phyloP_100V <- Variant_data_copy %>% cutpointr(x = phyloP100way_vertebrate, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 4.147 

CR_REVEL_score <- Variant_data_copy %>% cutpointr(x = REVEL_score, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 0.308 

CR_REVEL_rankscore <- Variant_data_copy %>% cutpointr(x = REVEL_rankscore, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 0.62947 

CR_SiPhy <- Variant_data_copy %>% cutpointr(x = SiPhy_29way_logOdds, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# 8.0188 







# final cutpoints!!! ------------------------------------------------------


# previously used code for cutpointr
#CR_BayesDel_AF <- Variant_data_copy %>% cutpointr(x = BayesDel_addAF_score, class = Type, pos_class = "Pathogenic", na.rm = TRUE)
# I used default method

# the two bootstrapping methods are: 'maximize_boot_metric' and 'minimize_boot_metric'
# another argument called 'boot_cut' needs to be specified when using either method
# also set summary_func to mean
# metric = accuracy
# boot_stratify = TRUE (to take equal amount of pos and neg during bootstrapping)



Variant_data <- read_csv("Extra_variants_3_unique.csv")
View(Variant_data)

# need to change the 'Type' data to factor - to allow grouping based on Pathogenic or Benign
Variant_data$Type <- as.factor(Variant_data$Type)
Variant_data$Syndrome <- as.factor(Variant_data$Syndrome)

#requires dataframe
Variant_data_copy <- as.data.frame(Variant_data)

# load in cutpointr
library(cutpointr)


# first trying with method=maximize_boot_metric
BayesDel_AF_maxbm <- Variant_data_copy %>% 
  cutpointr(x = BayesDel_addAF_score,
            class = Type,
            pos_class = "Pathogenic",
            method = maximize_boot_metric,
            boot_cut = 200, 
            summary_func = mean,
            metric = accuracy,
            boot_stratify = TRUE,
            na.rm = TRUE)

summary(BayesDel_AF_maxbm)
# Method: maximize_boot_metric 
# Predictor: BayesDel_addAF_score 
# Outcome: Type 
# Direction: >= 
#   
#   AUC   n n_pos n_neg
# 0.9862 262   128   134
# 
# optimal_cutpoint accuracy    acc sensitivity specificity  tp fn fp  tn
# -0.1002   0.9466 0.9466      0.9688      0.9254 124  4 10 124
# 
# Predictor summary: 
#   Data      Min.          5%    1st Qu.      Median        Mean    3rd Qu.          95%     Max.        SD NAs
# Overall -0.836670 -0.73430230 -0.5273238 -0.05503715 -0.08681692  0.3180443  0.618043000 0.625005 0.4526594   0
# Benign -0.836670 -0.78338650 -0.6449435 -0.51285300 -0.46863107 -0.3132670 -0.000442265 0.173542 0.2298310   0
# Pathogenic -0.326563 -0.05708148  0.1273445  0.33510350  0.31289477  0.5258318  0.625005000 0.625005 0.2256817   0






# CADD Phred 
CADD_phred_maxbm <- Variant_data_copy %>% 
  cutpointr(x = CADD_phred,
            class = Type,
            pos_class = "Pathogenic",
            method = maximize_boot_metric,
            boot_cut = 200, 
            summary_func = mean,
            metric = accuracy,
            boot_stratify = TRUE,
            na.rm = TRUE)

summary(CADD_phred_maxbm)
# Method: maximize_boot_metric 
# Predictor: CADD_phred 
# Outcome: Type 
# Direction: >= 
#   
#   AUC   n n_pos n_neg
# 0.937 262   128   134
# 
# optimal_cutpoint accuracy    acc sensitivity specificity  tp fn fp  tn
# 24.6818   0.8626 0.8626      0.8594      0.8657 110 18 18 116
# 
# Predictor summary: 
#   Data   Min.      5% 1st Qu. Median     Mean 3rd Qu.    95% Max.        SD NAs
# Overall  0.001  0.9453  16.720  24.50 22.49947  28.975 36.950   57 10.383783   0
# Benign  0.001  0.2341   7.668  16.87 15.52516  23.375 27.255   32  8.995921   0
# Pathogenic 13.950 23.0700  26.175  28.85 29.80070  33.000 39.650   57  5.638201   0




# also tried minimum_boot_metric, which very clearly is not appropriate for what I want!
# e.g. the CADD phred cutpoint was 0.0151 :o












# Determining filtering method -----------------------------------
# 14/11/23

# model 1: only one threshold must be met
# model 2: both thresholds must be met



model_1 <- subset(Variant_data_copy, BayesDel_addAF_score >= -0.1002 | CADD_phred >= 24.6818)
# 150 obs

model_2 <- subset(Variant_data_copy, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# 112 obs



summary(Variant_data_copy$Type)
# Benign Pathogenic 
# 134        128

summary(model_1$Type)
# Benign Pathogenic 
# 22        128

summary(model_2$Type)
# Benign Pathogenic 
# 6        106

# sensitivity (true positives)
# specificity (1- false negatives)

(1-(22/134))*100 
# model 1 = 100% sensitivity, and 83.58% specificity 
183.58/2
# 91.79 averaged

(106/128)*100
(1-(6/134))*100 
# model 2 = 82.81% sensitivity, and 95.52% specificity
(82.81+95.52)/2
# 89.165 averaged

# I think I would prefer to be more conservative (model 2)
# because both scoring methods have their flaws,
# and I'm looking at recessive inheritance, 
# so would prefer a cohort of high confidence disruptive variants (and less false positives)



summary(model_1)
# BayesDel_addAF_score
# Min.   :-0.50073    
# 1st Qu.: 0.07852    
# Median : 0.25499    
# Mean   : 0.24597    
# 3rd Qu.: 0.50425    
# Max.   : 0.62501

# CADD_phred   
# Min.   :13.95  
# 1st Qu.:25.68  
# Median :28.30  
# Mean   :29.32  
# 3rd Qu.:32.00  
# Max.   :57.00


summary(model_2)
# BayesDel_addAF_score
# Min.   :-0.0758     
# 1st Qu.: 0.1742     
# Median : 0.3486     
# Mean   : 0.3420     
# 3rd Qu.: 0.5484     
# Max.   : 0.6250

# CADD_phred   
# Min.   :25.10  
# 1st Qu.:27.10  
# Median :29.25  
# Mean   :30.88  
# 3rd Qu.:33.00  
# Max.   :57.00


# could plot variant types (could need to load in a different dataset for this though)


