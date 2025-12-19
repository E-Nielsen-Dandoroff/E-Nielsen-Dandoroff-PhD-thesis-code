## EntireCohort_MGSvar_analysis

# This script intends to take the Emily-folder/EntireCohort_ready.csv 
# dataset and conduct MGSvar statistical analysis.

# The code within this script was originally within:
# 'Emily-folder/MGSvar/EntireCohort_MGS_analysis.R' 
# 'Emily-folder/MGSvar/EntireCohort_MGSvar_script2.R' 


# Load packages ----------------------------------------------------------------

# load the following packages each time:
install.packages('readr')
library(readr)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
# and for saving tidy LR results_2
install.packages('broom')
library(broom)





# Annotate EntireCohort with MGSvar data ---------------------------------------
# the following code intends to annotate the EntireCohort_ready.csv dataframe
# with MGSvar data. 


# This section of code only requires execution once!
# Go straight to the next section of code for statistical analysis

####
# last updated & executed: 12/09/2024
####


# load EntireCohort dataframe
system('dx download Emily-folder/EntireCohort_ready.csv')
EntireCohort <- read_csv("EntireCohort_ready.csv")
EntireCohort <- EntireCohort[-1]

# load the lists of patients with MGS variants
system('dx download -r Emily-folder/MGSvar/IDs_VARs')

CDC45_MGSvar_IDs_VARs <- read_csv("IDs_VARs/CDC45_MGSvar_IDs_VARs.csv")
# 1629 obs
CDT1_MGSvar_IDs_VARs <- read_csv("IDs_VARs/CDT1_MGSvar_IDs_VARs.csv")
# 653 obs
DONSON_MGSvar_IDs_VARs <- read_csv("IDs_VARs/DONSON_MGSvar_IDs_VARs.csv")
#53 obs
GINS2_MGSvar_IDs_VARs <- read_csv("IDs_VARs/GINS2_MGSvar_IDs_VARs.csv")
# 21 obs
GINS3_MGSvar_IDs_VARs <- read_csv("IDs_VARs/GINS3_MGSvar_IDs_VARs.csv")
# 21 obs
MCM7_MGSvar_IDs_VARs <- read_csv("IDs_VARs/MCM7_MGSvar_IDs_VARs.csv")
# 54 obs
ORC1_MGSvar_IDs_VARs <- read_csv("IDs_VARs/ORC1_MGSvar_IDs_VARs.csv")
# 654 obs
ORC4_MGSvar_IDs_VARs <- read_csv("IDs_VARs/ORC4_MGSvar_IDs_VARs.csv")
# 3 obs
ORC6_MGSvar_IDs_VARs <- read_csv("IDs_VARs/ORC6_MGSvar_IDs_VARs.csv")
# 870 obs


# data cleaning
# remove the first column from all the GOI dataframes
CDC45_MGSvar_IDs_VARs <- CDC45_MGSvar_IDs_VARs[-1]
CDT1_MGSvar_IDs_VARs <- CDT1_MGSvar_IDs_VARs[-1]
DONSON_MGSvar_IDs_VARs <- DONSON_MGSvar_IDs_VARs[-1]
GINS2_MGSvar_IDs_VARs <- GINS2_MGSvar_IDs_VARs[-1]
GINS3_MGSvar_IDs_VARs <- GINS3_MGSvar_IDs_VARs[-1]
MCM7_MGSvar_IDs_VARs <- MCM7_MGSvar_IDs_VARs[-1]
ORC1_MGSvar_IDs_VARs <- ORC1_MGSvar_IDs_VARs[-1]
ORC4_MGSvar_IDs_VARs <- ORC4_MGSvar_IDs_VARs[-1]
ORC6_MGSvar_IDs_VARs <- ORC6_MGSvar_IDs_VARs[-1]

# make a list of the relevent eids
eid <- c(CDC45_MGSvar_IDs_VARs$FID, 
         CDT1_MGSvar_IDs_VARs$FID,
         DONSON_MGSvar_IDs_VARs$FID,
         GINS2_MGSvar_IDs_VARs$FID,
         GINS3_MGSvar_IDs_VARs$FID,
         MCM7_MGSvar_IDs_VARs$FID,
         ORC1_MGSvar_IDs_VARs$FID,
         ORC4_MGSvar_IDs_VARs$FID,
         ORC6_MGSvar_IDs_VARs$FID)

# make this into a dataframe
eids_df <- data.frame(eid, 1:1)
names(eids_df)[names(eids_df) == 'X1.1'] <- 'MGSvar'

# remove any duplicated eids
eids_df[duplicated(eids_df$eid),]
# 1649 5241008      1
# 1862 2955680      1
# 1958 1212794      1
# 2120 3209758      1
# 2443 3985909      1
# 2447 3202892      1
# 2472 2995708      1
# 2995 5984681      1
# 3080 4652048      1
# 3180 1084117      1
# 3473 3650868      1
# 3772 3678144      1
# 3816 2743983      1
# 3895 1633515      1
# 3915 5182400      1


eids_df <- eids_df[-c(1649,1862,1958,2120,2443,2447,2472,2995,3080,3180,3473,3772,3816,3895,3915), ]
# check (should be 0)
eids_df[duplicated(eids_df$eid),]


# add MGSvar column to the EntireCohort_1 dataframe by merging eids_df
EntireCohort_MGSvar <- left_join(EntireCohort, eids_df, by="eid")
summary(EntireCohort_MGSvar)
# looks correct





# before merging var info, should figure out what to do with the eIDs in common between the GOI dataframes

# use full_join() between the GOI dataframes
# then merge VAR columns using something like this: DF$C = paste(DF$A, DF$B, sep=";")

# first need to make sure we retain GOI information
# add a column to each specifying which GOI (use '1' as 'yes', in prep for glm())
CDC45_MGSvar_IDs_VARs$CDC45 <- 1
CDT1_MGSvar_IDs_VARs$CDT1 <- 1
DONSON_MGSvar_IDs_VARs$DONSON <- 1
GINS2_MGSvar_IDs_VARs$GINS2 <- 1
GINS3_MGSvar_IDs_VARs$GINS3 <- 1
MCM7_MGSvar_IDs_VARs$MCM7 <- 1
ORC1_MGSvar_IDs_VARs$ORC1 <- 1
ORC4_MGSvar_IDs_VARs$ORC4 <- 1
ORC6_MGSvar_IDs_VARs$ORC6 <-1 

# full_join() to keep all data
GOI_1 <- full_join(CDC45_MGSvar_IDs_VARs, CDT1_MGSvar_IDs_VARs, by="FID")
GOI_2 <- full_join(GOI_1, DONSON_MGSvar_IDs_VARs, by="FID")
GOI_3 <- full_join(GOI_2, GINS2_MGSvar_IDs_VARs, by="FID")
GOI_4 <- full_join(GOI_3, GINS3_MGSvar_IDs_VARs, by="FID")
GOI_5 <- full_join(GOI_4, MCM7_MGSvar_IDs_VARs, by="FID")
GOI_6 <- full_join(GOI_5, ORC1_MGSvar_IDs_VARs, by="FID")
GOI_7 <- full_join(GOI_6, ORC4_MGSvar_IDs_VARs, by="FID")
GOI_8 <- full_join(GOI_7, ORC6_MGSvar_IDs_VARs, by="FID")

# now to merge/sum columns together
# The GOI name columns can remain as they are
# but the 'Sums' columns should be summed together to give an indication of variant load


GOI_8$variant_load <- rowSums(GOI_8[, c("Sums", "Sums.x", "Sums.y", 
                                        "Sums.x.x", "Sums.y.y", 
                                        "Sums.x.x.x", "Sums.y.y.y", 
                                        "Sums.x.x.x.x", "Sums.y.y.y.y")], 
                              na.rm = TRUE)

summary(GOI_8$variant_load)
# min = 1
# max = 2

# also want to merge the VAR data using something like this: DF$C = paste(DF$A, DF$B, sep=";")

GOI_8$variants <- paste(GOI_8$VAR.x, GOI_8$VAR.y, sep = ";")

# bar <- apply(cbind(1:4, foo), 1, 
#              function(x) paste(x[!is.na(x)], collapse = ", "))

GOI_test <- GOI_8

# use a tidyverse method
GOI_8 <- GOI_8 %>% unite(., col = "variants",  
                         VAR, VAR.x, VAR.y, VAR.x.x, VAR.y.y, VAR.x.x.x, VAR.y.y.y, VAR.x.x.x.x, VAR.y.y.y.y, 
                         na.rm=TRUE, sep = ";")
# perfect :)

# now I just need to remove all of the Sums columns and this dataframe will be ready to combine with EntireCohort
GOI_final <- GOI_8[c("FID", "variants", "variant_load", "CDC45", "CDT1", "DONSON", "GINS2", "GINS3", "MCM7", "ORC1", "ORC4", "ORC6")]

summary(GOI_final)
# great!


# now to add the GOI-final data to the EntireCohort dataframe!
GOI_final
names(GOI_final)[names(GOI_final) == 'FID'] <- 'eid'

EntireCohort_MGSvar <- left_join(EntireCohort_MGSvar, GOI_final, by="eid")
# looks good to me


## a couple more things to do before saving the data
# change GOI column NA values to 0, in prep for glm() analysis

EntireCohort_MGSvar$CDC45[is.na(EntireCohort_MGSvar$CDC45)] <- 0
EntireCohort_MGSvar$CDT1[is.na(EntireCohort_MGSvar$CDT1)] <- 0
EntireCohort_MGSvar$DONSON[is.na(EntireCohort_MGSvar$DONSON)] <- 0
EntireCohort_MGSvar$GINS2[is.na(EntireCohort_MGSvar$GINS2)] <- 0
EntireCohort_MGSvar$GINS3[is.na(EntireCohort_MGSvar$GINS3)] <- 0
EntireCohort_MGSvar$MCM7[is.na(EntireCohort_MGSvar$MCM7)] <- 0
EntireCohort_MGSvar$ORC1[is.na(EntireCohort_MGSvar$ORC1)] <- 0
EntireCohort_MGSvar$ORC4[is.na(EntireCohort_MGSvar$ORC4)] <- 0
EntireCohort_MGSvar$ORC6[is.na(EntireCohort_MGSvar$ORC6)] <- 0
# check
summary(EntireCohort_MGSvar)

# also want to do this for the variant load and MGSvar columns
EntireCohort_MGSvar$variant_load[is.na(EntireCohort_MGSvar$variant_load)] <- 0
EntireCohort_MGSvar$MGSvar[is.na(EntireCohort_MGSvar$MGSvar)] <- 0




# save this to MGSvar folder
# first save it to R
write.csv(EntireCohort_MGSvar, "EntireCohort_MGSvar.csv")
# now to save this to UKB project
system('dx upload EntireCohort_MGSvar.csv --path Emily-folder/MGSvar/EntireCohort_MGSvar.csv')










# MGSvar Incidence --------------------------------------------------------

# Using data.table to load in the EntireCohort_MGS dataset
EntireCohort_MGSvar <- fread('/mnt/project/Emily-folder/MGSvar/EntireCohort_MGSvar.csv', data.table = FALSE)

# The dataframe gets loaded in with an extra column called 'V1' which I don't want
EntireCohort_MGSvar <- EntireCohort_MGSvar[-1]
# removed :)

EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)
EntireCohort_ready <- EntireCohort_ready[,c(2, 57:72)]

EntireCohort_MGSvar <- left_join(EntireCohort_MGSvar, EntireCohort_ready, by = 'eid')







# look at all datafields
colnames(EntireCohort_MGSvar)
## fields that I won't use:
# [1] "eid"
# [3] "Ethnic_background"
# [4] "Ever_smoked"                                 
# [5] "Smoking_status"
# [6] "Genetic_relatedness_pairing"                 
# [7] "Year_of_birth"                               
# [8] "Reported_occurrencces_of_cancer"
# [12] "Sex"
# [14] "Qualifications" 
# [15] "Average_household_income"                    
# [16] "Number_days_per_week_10min_physical_activity"

## fields that I will use in an analysis but aren't relevant to compare with cancer incidence (y/n)
# [9] "Age_at_cancer_diagnosis_earliest"                     
# [10] "ICD_code" 
# [19] "variants" 

## fields that I can compare with cancer incidence (y/n)    # [34] "cancer_ROOC_dummy"
# [2] "Age_at_recruitment"
# [13] "BMI"
# [17] "Standing_height"                             
# [18] "MGSvar"                                      
# [20] "variant_load"                                
# [21] "CDC45"                                       
# [22] "CDT1"                                        
# [23] "DONSON"                                      
# [24] "GINS2"                                       
# [25] "GINS3"                                       
# [26] "MCM7"                                        
# [27] "ORC1"                                        
# [28] "ORC4"                                        
# [29] "ORC6"                                        
# [30] "Ethnicity_White"                             
# [31] "Ever_smoked_dummy"                           
# [32] "Smoker_previous"                             
# [33] "Smoker_current"                              
# [35] "Sex_Male"                                    
# [36] "University_education"                        
# [37] "Annual_income"                               
# [38] "Weekly_exercise"









# look at whether these fields significantly effect cancer incidence (y/n)
# univariate logistic regression to assess relevance of confounding variables:

# cancer_ROOC_dummy is a 0/1 field - so can do logistic regression

# # reference code for logistic regression:
# LR_ <- glm(formula = cancer_ROOC_dummy ~ p_1 + p_2 + p_3, data = EntireCohort_MGSvar, family = binomial)
# summary(LR_)
# exp(cbind(OR = coef(LR_), confint(LR_)))

## age at recruitment
LR_1 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment, data = EntireCohort_MGSvar, family = binomial)
# summary(LR_1)
# exp(cbind(OR = coef(LR_1), confint(LR_1)))
ALL_incidence_1 <- tidy(LR_1, conf.int = TRUE, exp = TRUE)
write.csv(ALL_incidence_1, "ALL_incidence_1.csv", row.names=FALSE)
system('dx upload ALL_incidence_1.csv --path Emily-folder/results_2/ALL_incidence_1.csv')


## BMI
LR_2 <- glm(formula = cancer_ROOC_dummy ~ BMI, data = EntireCohort_MGSvar, family = binomial)
# summary(LR_2)
# exp(cbind(OR = coef(LR_2), confint(LR_2)))
ALL_incidence_2 <- tidy(LR_2, conf.int = TRUE, exp = TRUE)
write.csv(ALL_incidence_2, "ALL_incidence_2.csv", row.names=FALSE)
system('dx upload ALL_incidence_2.csv --path Emily-folder/results_2/ALL_incidence_2.csv')


## Standing_height
LR_3 <- glm(formula = cancer_ROOC_dummy ~ Standing_height, data = EntireCohort_MGSvar, family = binomial)
# summary(LR_3)
# exp(cbind(OR = coef(LR_3), confint(LR_3)))
ALL_incidence_3 <- tidy(LR_3, conf.int = TRUE, exp = TRUE)
write.csv(ALL_incidence_3, "ALL_incidence_3.csv", row.names=FALSE)
system('dx upload ALL_incidence_3.csv --path Emily-folder/results_2/ALL_incidence_3.csv')


## MGSvar
LR_4 <- glm(formula = cancer_ROOC_dummy ~ MGSvar, data = EntireCohort_MGSvar, family = binomial)
# summary(LR_4)
# exp(cbind(OR = coef(LR_4), confint(LR_4)))
MGSvar_incidence_1 <- tidy(LR_4, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_incidence_1, "MGSvar_incidence_1.csv", row.names=FALSE)
system('dx upload MGSvar_incidence_1.csv --path Emily-folder/results_2/MGSvar_incidence_1.csv')


## variant_load
LR_5 <- glm(formula = cancer_ROOC_dummy ~ variant_load, data = EntireCohort_MGSvar, family = binomial)
# summary(LR_5)
# exp(cbind(OR = coef(LR_5), confint(LR_5)))
MGSvar_incidence_2 <- tidy(LR_5, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_incidence_2, "MGSvar_incidence_2.csv", row.names=FALSE)
system('dx upload MGSvar_incidence_2.csv --path Emily-folder/results_2/MGSvar_incidence_2.csv')


## the GOIs: "CDC45", "CDT1", "DONSON", "GINS2", "GINS3", "MCM7", "ORC1", "ORC4", "ORC6"
LR_6 <- glm(formula = cancer_ROOC_dummy ~ CDC45+CDT1+DONSON+GINS2+GINS3+MCM7+ORC1+ORC4+ORC6, data = EntireCohort_MGSvar, family = binomial)
# summary(LR_6)
# exp(cbind(OR = coef(LR_6), confint(LR_6)))
MGSvar_incidence_3 <- tidy(LR_6, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_incidence_3, "MGSvar_incidence_3.csv", row.names=FALSE)
system('dx upload MGSvar_incidence_3.csv --path Emily-folder/results_2/MGSvar_incidence_3.csv')


# for now I'm not going to save these as individual results_2 - as I can see from the code that 
# the significance doesn't change when together and apart!
## the GOIs: "CDC45"
LR_6_CDC45 <- glm(formula = cancer_ROOC_dummy ~ CDC45, data = EntireCohort_MGSvar, family = binomial)
summary(LR_6_CDC45)
exp(cbind(OR = coef(LR_6_CDC45), confint(LR_6_CDC45)))
## the GOIs: "CDT1"
LR_6_CDT1 <- glm(formula = cancer_ROOC_dummy ~ CDT1, data = EntireCohort_MGSvar, family = binomial)
summary(LR_6_CDT1)
exp(cbind(OR = coef(LR_6_CDT1), confint(LR_6_CDT1)))
## the GOIs: "DONSON"
LR_6_DONSON <- glm(formula = cancer_ROOC_dummy ~ DONSON, data = EntireCohort_MGSvar, family = binomial)
summary(LR_6_DONSON)
exp(cbind(OR = coef(LR_6_DONSON), confint(LR_6_DONSON)))
## the GOIs: "GINS2"
LR_6_GINS2 <- glm(formula = cancer_ROOC_dummy ~ GINS2, data = EntireCohort_MGSvar, family = binomial)
summary(LR_6_GINS2)
exp(cbind(OR = coef(LR_6_GINS2), confint(LR_6_GINS2)))
## the GOIs: "GINS3"
LR_6_GINS3 <- glm(formula = cancer_ROOC_dummy ~ GINS3, data = EntireCohort_MGSvar, family = binomial)
summary(LR_6_GINS3)
exp(cbind(OR = coef(LR_6_GINS3), confint(LR_6_GINS3)))
## the GOIs: "MCM7"
LR_6_MCM7 <- glm(formula = cancer_ROOC_dummy ~ MCM7, data = EntireCohort_MGSvar, family = binomial)
summary(LR_6_MCM7)
exp(cbind(OR = coef(LR_6_MCM7), confint(LR_6_MCM7)))
## the GOIs: "ORC1"
LR_6_ORC1 <- glm(formula = cancer_ROOC_dummy ~ ORC1, data = EntireCohort_MGSvar, family = binomial)
summary(LR_6_ORC1)
exp(cbind(OR = coef(LR_6_ORC1), confint(LR_6_ORC1)))
## the GOIs: "ORC4"
LR_6_ORC4 <- glm(formula = cancer_ROOC_dummy ~ ORC4, data = EntireCohort_MGSvar, family = binomial)
summary(LR_6_ORC4)
exp(cbind(OR = coef(LR_6_ORC4), confint(LR_6_ORC4)))
## the GOIs: "ORC6"
LR_6_ORC6 <- glm(formula = cancer_ROOC_dummy ~ ORC6, data = EntireCohort_MGSvar, family = binomial)
summary(LR_6_ORC6)
exp(cbind(OR = coef(LR_6_ORC6), confint(LR_6_ORC6)))


# [30] "Ethnicity_White"                             
LR_7 <- glm(formula = cancer_ROOC_dummy ~ Ethnicity_White, data = EntireCohort_MGSvar, family = binomial)
# summary(LR_7)
# exp(cbind(OR = coef(LR_7), confint(LR_7)))
ALL_incidence_4 <- tidy(LR_7, conf.int = TRUE, exp = TRUE)
write.csv(ALL_incidence_4, "ALL_incidence_4.csv", row.names=FALSE)
system('dx upload ALL_incidence_4.csv --path Emily-folder/results_2/ALL_incidence_4.csv')


# Ever_smoked_dummy
LR_8 <- glm(formula = cancer_ROOC_dummy ~ Ever_smoked_dummy, data = EntireCohort_MGSvar, family = binomial)
# summary(LR_8)
# exp(cbind(OR = coef(LR_8), confint(LR_8)))
ALL_incidence_5 <- tidy(LR_8, conf.int = TRUE, exp = TRUE)
write.csv(ALL_incidence_5, "ALL_incidence_5.csv", row.names=FALSE)
system('dx upload ALL_incidence_5.csv --path Emily-folder/results_2/ALL_incidence_5.csv')


# smoking status: "Smoker_previous", "Smoker_current"
LR_9 <- glm(formula = cancer_ROOC_dummy ~ Smoker_previous+Smoker_current, data = EntireCohort_MGSvar, family = binomial)
# summary(LR_9)
# exp(cbind(OR = coef(LR_9), confint(LR_9)))
ALL_incidence_6 <- tidy(LR_9, conf.int = TRUE, exp = TRUE)
write.csv(ALL_incidence_6, "ALL_incidence_6.csv", row.names=FALSE)
system('dx upload ALL_incidence_6.csv --path Emily-folder/results_2/ALL_incidence_6.csv')


## "Sex_Male"
LR_10 <- glm(formula = cancer_ROOC_dummy ~ Sex_Male, data = EntireCohort_MGSvar, family = binomial)
# summary(LR_10)
# exp(cbind(OR = coef(LR_10), confint(LR_10)))
ALL_incidence_7 <- tidy(LR_10, conf.int = TRUE, exp = TRUE)
write.csv(ALL_incidence_7, "ALL_incidence_7.csv", row.names=FALSE)
system('dx upload ALL_incidence_7.csv --path Emily-folder/results_2/ALL_incidence_7.csv')


## "University_education"
LR_11 <- glm(formula = cancer_ROOC_dummy ~ University_education, data = EntireCohort_MGSvar, family = binomial)
# summary(LR_11)
# exp(cbind(OR = coef(LR_11), confint(LR_11)))
ALL_incidence_8 <- tidy(LR_11, conf.int = TRUE, exp = TRUE)
write.csv(ALL_incidence_8, "ALL_incidence_8.csv", row.names=FALSE)
system('dx upload ALL_incidence_8.csv --path Emily-folder/results_2/ALL_incidence_8.csv')


## Annual_income
LR_12 <- glm(formula = cancer_ROOC_dummy ~ Annual_income, data = EntireCohort_MGSvar, family = binomial)
# summary(LR_12)
# exp(cbind(OR = coef(LR_12), confint(LR_12)))
ALL_incidence_9 <- tidy(LR_12, conf.int = TRUE, exp = TRUE)
write.csv(ALL_incidence_9, "ALL_incidence_9.csv", row.names=FALSE)
system('dx upload ALL_incidence_9.csv --path Emily-folder/results_2/ALL_incidence_9.csv')


## Weekly_exercise
LR_13 <- glm(formula = cancer_ROOC_dummy ~ Weekly_exercise, data = EntireCohort_MGSvar, family = binomial)
# summary(LR_13)
# exp(cbind(OR = coef(LR_13), confint(LR_13)))
ALL_incidence_10 <- tidy(LR_13, conf.int = TRUE, exp = TRUE)
write.csv(ALL_incidence_10, "ALL_incidence_10.csv", row.names=FALSE)
system('dx upload ALL_incidence_10.csv --path Emily-folder/results_2/ALL_incidence_10.csv')















### Multivariable logistic regression

# should only use one smoking variable - deciding to go with 'Ever_smoked_dummy' 
# because, while previous smokers have a slightly increased association with cancer,
# both show a significant assocaiton. Im happy to group current and previous.
# Ever_smoked_dummy also has less missing data!

# Response variable = "cancer_ROOC_dummy"

# stable predictor variables:
# "Age_at_recruitment"
# "BMI"
# "Standing_height"
# "Ethnicity_White"                             
# "Ever_smoked_dummy"                                
# "Sex_Male"                                    
# "University_education"                        
# "Annual_income"                               
# "Weekly_exercise"

# predictor variables to test:
# "MGSvar"                                      
# "variant_load"                                
# "CDC45"                                       
# "CDT1"                                        
# "DONSON"                                      
# "GINS2"                                       
# "GINS3"                                       
# "MCM7"                                        
# "ORC1"                                        
# "ORC4"                                        
# "ORC6"                                     

## make 2 sets of covariates
# set A (the stronger cancer-associated factors - that do not have as much missing data)
# = Age_at_recruitment + Ethnicity_White + Ever_smoked_dummy + Sex_Male
# set B (all factors)



# want to compare each of the MGSvar fields when univariate vs multivariate
# so need to adjust the df for set A and set B - to remove missing data 

data_setA <- EntireCohort_MGSvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                    "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                    "MGSvar", "variant_load", "CDC45", "CDT1", "DONSON",
                                    "GINS2", "GINS3", "MCM7", "ORC1", "ORC4", "ORC6")]
data_setA <- data_setA[complete.cases(data_setA), ]
dim(data_setA)
# 465881     17






# predictor variables without any MGSvar variable (set A)
LR_14 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male, 
             data = data_setA, family = binomial)
# summary(LR_14)
# exp(cbind(OR = coef(LR_14), confint(LR_14)))
ALL_incidence_11 <- tidy(LR_14, conf.int = TRUE, exp = TRUE)
write.csv(ALL_incidence_11, "ALL_incidence_11.csv", row.names=FALSE)
system('dx upload ALL_incidence_11.csv --path Emily-folder/results_2/ALL_incidence_11.csv')


# MGSvar without confounding variables (set A)
LR_15 <- glm(formula = cancer_ROOC_dummy ~ MGSvar, 
             data = data_setA, family = binomial)
# summary(LR_15)
# exp(cbind(OR = coef(LR_15), confint(LR_15)))
MGSvar_incidence_4 <- tidy(LR_15, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_incidence_4, "MGSvar_incidence_4.csv", row.names=FALSE)
system('dx upload MGSvar_incidence_4.csv --path Emily-folder/results_2/MGSvar_incidence_4.csv')


# "MGSvar" (set A)
LR_16 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSvar, 
             data = data_setA, family = binomial)
# summary(LR_16)
# exp(cbind(OR = coef(LR_16), confint(LR_16)))
MGSvar_incidence_5 <- tidy(LR_16, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_incidence_5, "MGSvar_incidence_5.csv", row.names=FALSE)
system('dx upload MGSvar_incidence_5.csv --path Emily-folder/results_2/MGSvar_incidence_5.csv')



# variant_load without confounding variables (set A)
LR_17 <- glm(formula = cancer_ROOC_dummy ~ variant_load, 
             data = data_setA, family = binomial)
# summary(LR_17)
# exp(cbind(OR = coef(LR_17), confint(LR_17)))
MGSvar_incidence_6 <- tidy(LR_17, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_incidence_6, "MGSvar_incidence_6.csv", row.names=FALSE)
system('dx upload MGSvar_incidence_6.csv --path Emily-folder/results_2/MGSvar_incidence_6.csv')


# "variant_load"  (set A)
LR_18 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+variant_load, 
             data = data_setA, family = binomial)
# summary(LR_18)
# exp(cbind(OR = coef(LR_18), confint(LR_18)))
MGSvar_incidence_7 <- tidy(LR_18, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_incidence_7, "MGSvar_incidence_7.csv", row.names=FALSE)
system('dx upload MGSvar_incidence_7.csv --path Emily-folder/results_2/MGSvar_incidence_7.csv')



# GOIs without confounding variables (set A)
LR_19 <- glm(formula = cancer_ROOC_dummy ~ CDC45 + CDT1 + DONSON + GINS2 + GINS3 + MCM7 + ORC1 + ORC4 + ORC6, 
             data = data_setA, family = binomial)
# summary(LR_19)
# exp(cbind(OR = coef(LR_19), confint(LR_19)))
MGSvar_incidence_8 <- tidy(LR_19, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_incidence_8, "MGSvar_incidence_8.csv", row.names=FALSE)
system('dx upload MGSvar_incidence_8.csv --path Emily-folder/results_2/MGSvar_incidence_8.csv')



# "CDC45"                                       
# "CDT1"                                        
# "DONSON"                                      
# "GINS2"                                       
# "GINS3"                                       
# "MCM7"                                        
# "ORC1"                                        
# "ORC4"                                        
# "ORC6" (set A)
LR_20 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CDC45+CDT1+DONSON+GINS2+GINS3+MCM7+ORC1+ORC4+ORC6, 
             data = data_setA, family = binomial)
# summary(LR_20)
# exp(cbind(OR = coef(LR_20), confint(LR_20)))
MGSvar_incidence_9 <- tidy(LR_20, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_incidence_9, "MGSvar_incidence_9.csv", row.names=FALSE)
system('dx upload MGSvar_incidence_9.csv --path Emily-folder/results_2/MGSvar_incidence_9.csv')


# won't save this one - just here to look at in R
# just CDT1
LR_20_CDT1 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment + Ethnicity_White + 
                    Ever_smoked_dummy + Sex_Male + CDT1, 
                  data = data_setA, family = binomial)
# summary(LR_20_CDT1)
# exp(cbind(OR = coef(LR_20_CDT1), confint(LR_20_CDT1)))













# want to compare each of the MGSvar fields when univariate vs multivariate
# so need to adjust the df for set A and set B - to remove missing data 

data_setB <- EntireCohort_MGSvar[,c("eid", "cancer_ROOC_dummy", "Age_at_recruitment",
                                    "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male",
                                    "BMI", "Standing_height", "University_education",
                                    "Annual_income", "Weekly_exercise",
                                    "MGSvar", "variant_load", "CDC45", "CDT1", "DONSON",
                                    "GINS2", "GINS3", "MCM7", "ORC1", "ORC6")]
data_setB <- data_setB[complete.cases(data_setB), ]
dim(data_setB)
# 381159     22

# check n for each GOI
with(data_setB, table(ORC6))	
# DONSON = 34
# GINS2 = 20
# GINS3 = 18
# MCM7 = 46
# ORC4 = 2 (removed - but still in MGSvar field)


# predictor variables without any MGSvar variable (set B)
LR_21 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise, 
             data = data_setB, family = binomial)
# summary(LR_21)
# exp(cbind(OR = coef(LR_21), confint(LR_21)))
ALL_incidence_12 <- tidy(LR_21, conf.int = TRUE, exp = TRUE)
write.csv(ALL_incidence_12, "ALL_incidence_12.csv", row.names=FALSE)
system('dx upload ALL_incidence_12.csv --path Emily-folder/results_2/ALL_incidence_12.csv')



# "MGSvar" without any confounding variables (set B)
LR_22 <- glm(formula = cancer_ROOC_dummy ~ MGSvar, 
             data = data_setB, family = binomial)
# summary(LR_22)
# exp(cbind(OR = coef(LR_22), confint(LR_22)))
MGSvar_incidence_10 <- tidy(LR_22, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_incidence_10, "MGSvar_incidence_10.csv", row.names=FALSE)
system('dx upload MGSvar_incidence_10.csv --path Emily-folder/results_2/MGSvar_incidence_10.csv')



# "MGSvar" (set B)
LR_23 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+MGSvar, 
             data = data_setB, family = binomial)
# summary(LR_23)
# exp(cbind(OR = coef(LR_23), confint(LR_23)))
MGSvar_incidence_11 <- tidy(LR_23, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_incidence_11, "MGSvar_incidence_11.csv", row.names=FALSE)
system('dx upload MGSvar_incidence_11.csv --path Emily-folder/results_2/MGSvar_incidence_11.csv')




# "variant_load" without any confounding variables (set B)
LR_24 <- glm(formula = cancer_ROOC_dummy ~ variant_load, 
             data = data_setB, family = binomial)
# summary(LR_24)
# exp(cbind(OR = coef(LR_24), confint(LR_24)))
MGSvar_incidence_12 <- tidy(LR_24, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_incidence_12, "MGSvar_incidence_12.csv", row.names=FALSE)
system('dx upload MGSvar_incidence_12.csv --path Emily-folder/results_2/MGSvar_incidence_12.csv')


# "variant_load"  (set B)
LR_25 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+variant_load, 
             data = data_setB, family = binomial)
# summary(LR_25)
# exp(cbind(OR = coef(LR_25), confint(LR_25)))
MGSvar_incidence_13 <- tidy(LR_25, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_incidence_13, "MGSvar_incidence_13.csv", row.names=FALSE)
system('dx upload MGSvar_incidence_13.csv --path Emily-folder/results_2/MGSvar_incidence_13.csv')



# GOIs without any confounding variables (set B)
LR_26 <- glm(formula = cancer_ROOC_dummy ~ CDC45 + CDT1 + DONSON + GINS2 + GINS3 + MCM7 + ORC1 + ORC6, 
             data = data_setB, family = binomial)
# summary(LR_26)
# exp(cbind(OR = coef(LR_26), confint(LR_26)))
MGSvar_incidence_14 <- tidy(LR_26, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_incidence_14, "MGSvar_incidence_14.csv", row.names=FALSE)
system('dx upload MGSvar_incidence_14.csv --path Emily-folder/results_2/MGSvar_incidence_14.csv')


# "CDC45"                                       
# "CDT1"                                        
# "DONSON"                                      
# "GINS2"                                       
# "GINS3"                                       
# "MCM7"                                        
# "ORC1"                                        
# "ORC6" (set B)
LR_27 <- glm(formula = cancer_ROOC_dummy ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+CDC45+CDT1+DONSON+GINS2+GINS3+MCM7+ORC1+ORC6, 
             data = data_setB, family = binomial)
# summary(LR_27)
# exp(cbind(OR = coef(LR_27), confint(LR_27)))
MGSvar_incidence_15 <- tidy(LR_27, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_incidence_15, "MGSvar_incidence_15.csv", row.names=FALSE)
system('dx upload MGSvar_incidence_15.csv --path Emily-folder/results_2/MGSvar_incidence_15.csv')




























### Fisher test on these datasets

### starting with EntireCohort_MGSvar

with(EntireCohort_MGSvar, table(MGSvar, cancer_ROOC_dummy))
fisher.test(with(EntireCohort_MGSvar, table(MGSvar, cancer_ROOC_dummy)))
# p-value = 0.1471

with(EntireCohort_MGSvar, table(MGSvar, Ethnicity_White))
fisher.test(with(EntireCohort_MGSvar, table(MGSvar, Ethnicity_White)))
# p-value = 0.008923
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.7356580 0.9581976
# sample estimates:
#   odds ratio 
# 0.8379687

with(EntireCohort_MGSvar, table(MGSvar, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_MGSvar, table(MGSvar, Ever_smoked_dummy)))
# p-value = 0.1501

with(EntireCohort_MGSvar, table(MGSvar, Sex_Male))			
fisher.test(with(EntireCohort_MGSvar, table(MGSvar, Sex_Male)))
# p-value = 0.03841
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.8774304 0.9965520
# sample estimates:
#   odds ratio 
# 0.9351268

with(EntireCohort_MGSvar, table(MGSvar, University_education))			
fisher.test(with(EntireCohort_MGSvar, table(MGSvar, University_education)))
# p-value = 0.07436

with(EntireCohort_MGSvar, table(MGSvar, Annual_income))			
fisher.test(with(EntireCohort_MGSvar, table(MGSvar, Annual_income))	)
# p-value = 0.4618

with(EntireCohort_MGSvar, table(MGSvar, Smoker_previous))			
fisher.test(with(EntireCohort_MGSvar, table(MGSvar, Smoker_previous)))
# p-value = 0.6252

with(EntireCohort_MGSvar, table(MGSvar, Smoker_current))			
fisher.test(with(EntireCohort_MGSvar, table(MGSvar, Smoker_current))	)
# p-value = 0.6574



### data_setA
with(data_setA, table(MGSvar, cancer_ROOC_dummy))			
fisher.test(with(data_setA, table(MGSvar, cancer_ROOC_dummy)))
# p-value = 0.1402

with(data_setA, table(MGSvar, Ethnicity_White))			
fisher.test(with(data_setA, table(MGSvar, Ethnicity_White)))
# p-value = 0.007748
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.7317101 0.9536117
# sample estimates:
#   odds ratio 
# 0.8337248

with(data_setA, table(MGSvar, Ever_smoked_dummy))			
fisher.test(with(data_setA, table(MGSvar, Ever_smoked_dummy)))
# p-value = 0.1639

with(data_setA, table(MGSvar, Sex_Male))			
fisher.test(with(data_setA, table(MGSvar, Sex_Male)))
# p-value = 0.04746
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.8798296 0.9997692
# sample estimates:
#   odds ratio 
# 0.9379211



### data_setB
with(data_setB, table(MGSvar, cancer_ROOC_dummy))			
fisher.test(with(data_setB, table(MGSvar, cancer_ROOC_dummy)))
# p-value = 0.2284

with(data_setB, table(MGSvar, Ethnicity_White))			
fisher.test(with(data_setB, table(MGSvar, Ethnicity_White)))
# p-value = 0.3392

with(data_setB, table(MGSvar, Ever_smoked_dummy))			
fisher.test(with(data_setB, table(MGSvar, Ever_smoked_dummy)))
# p-value = 0.6353

with(data_setB, table(MGSvar, Sex_Male))			
fisher.test(with(data_setB, table(MGSvar, Sex_Male)))
# p-value = 0.1986

with(data_setB, table(MGSvar, University_education))			
fisher.test(with(data_setB, table(MGSvar, University_education)))
# p-value = 0.09777

with(data_setB, table(MGSvar, Annual_income))			
fisher.test(with(data_setB, table(MGSvar, Annual_income)))
# p-value = 0.5309






### can plot the continuous variables: 

# make MGSvar into factor
EntireCohort_MGSvar$MGSvar <- as.factor(EntireCohort_MGSvar$MGSvar)

# try a violin plot

# Age_at_recruitment  
EntireCohort_MGSvar %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
EntireCohort_MGSvar %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ MGSvar, data = EntireCohort_MGSvar)
# W = 909424509, p-value = 0.2877

# BMI  
EntireCohort_MGSvar %>% 
  ggplot(mapping = aes(x = MGSvar, y = BMI)) +
  geom_violin() +
  theme_light()
EntireCohort_MGSvar %>% 
  ggplot(mapping = aes(x = MGSvar, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ MGSvar, data = EntireCohort_MGSvar)
# W = 904798292, p-value = 0.4939

# Standing_height 
EntireCohort_MGSvar %>% 
  ggplot(mapping = aes(x = MGSvar, y = Standing_height)) +
  geom_violin() +
  theme_light()
EntireCohort_MGSvar %>% 
  ggplot(mapping = aes(x = MGSvar, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ MGSvar, data = EntireCohort_MGSvar)
# W = 941500846, p-value = 0.000476

# make this a bar plot
# Weekly_exercise
EntireCohort_MGSvar %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(MGSvar), scales = "free_y") +
  theme_light()
# 



### data_setA

# make MGSvar into factor
data_setA$MGSvar <- as.factor(data_setA$MGSvar)

# Age_at_recruitment  
data_setA %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setA %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ MGSvar, data = data_setA)
# W = 895083818, p-value = 0.3214



### data_setB

# make MGSvar into factor
data_setB$MGSvar <- as.factor(data_setB$MGSvar)

# Age_at_recruitment  
data_setB %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Age_at_recruitment ~ MGSvar, data = data_setB)
# W = 595401695, p-value = 0.5883

# BMI  
data_setB %>% 
  ggplot(mapping = aes(x = MGSvar, y = BMI)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = MGSvar, y = BMI)) +
  geom_boxplot() +
  theme_light()
wilcox.test(BMI ~ MGSvar, data = data_setB)
# W = 600130190, p-value = 0.8214

# Standing_height 
data_setB %>% 
  ggplot(mapping = aes(x = MGSvar, y = Standing_height)) +
  geom_violin() +
  theme_light()
data_setB %>% 
  ggplot(mapping = aes(x = MGSvar, y = Standing_height)) +
  geom_boxplot() +
  theme_light()
wilcox.test(Standing_height ~ MGSvar, data = data_setB)
# W = 614421662, p-value = 0.01095


# make this a bar plot
# Weekly_exercise
data_setB %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(MGSvar), scales = "free_y") +
  theme_light()
# 




### fisher test for GOIs
# only do this for the GOIs that look interesting!



### could also do the plots + wilcox tests for these genes

















# MGSvar Type -------------------------------------------------------------





## won't do nnet package multinomial LR because it seemed a bit inconsistant in previous analyses

# LR for each cancer category

# Neoplasm_oral
LR_oral <- glm(formula = Neoplasm_oral ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSvar, 
               data = EntireCohort_MGSvar, family = binomial)
# summary(LR_oral)
# exp(cbind(OR = coef(LR_oral), confint(LR_oral)))
MGSvar_type_1 <- tidy(LR_oral, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_type_1, "MGSvar_type_1.csv", row.names=FALSE)
system('dx upload MGSvar_type_1.csv --path Emily-folder/results_2/MGSvar_type_1.csv')


# Neoplasm_digestive
LR_digestive <- glm(formula = Neoplasm_digestive ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSvar, 
                    data = EntireCohort_MGSvar, family = binomial)
# summary(LR_digestive)
# exp(cbind(OR = coef(LR_digestive), confint(LR_digestive)))
MGSvar_type_2 <- tidy(LR_digestive, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_type_2, "MGSvar_type_2.csv", row.names=FALSE)
system('dx upload MGSvar_type_2.csv --path Emily-folder/results_2/MGSvar_type_2.csv')


# Neoplasm_respiratory
LR_respiratory <- glm(formula = Neoplasm_respiratory ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSvar, 
                      data = EntireCohort_MGSvar, family = binomial)
# summary(LR_respiratory)
# exp(cbind(OR = coef(LR_respiratory), confint(LR_respiratory)))
MGSvar_type_3 <- tidy(LR_respiratory, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_type_3, "MGSvar_type_3.csv", row.names=FALSE)
system('dx upload MGSvar_type_3.csv --path Emily-folder/results_2/MGSvar_type_3.csv')


# Neoplasm_bone
LR_bone <- glm(formula = Neoplasm_bone ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSvar, 
               data = EntireCohort_MGSvar, family = binomial)
# summary(LR_bone)
# exp(cbind(OR = coef(LR_bone), confint(LR_bone)))
MGSvar_type_4 <- tidy(LR_bone, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_type_4, "MGSvar_type_4.csv", row.names=FALSE)
system('dx upload MGSvar_type_4.csv --path Emily-folder/results_2/MGSvar_type_4.csv')


# Neoplasm_skin
LR_skin <- glm(formula = Neoplasm_skin ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSvar, 
               data = EntireCohort_MGSvar, family = binomial)
# summary(LR_skin)
# exp(cbind(OR = coef(LR_skin), confint(LR_skin)))
MGSvar_type_5 <- tidy(LR_skin, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_type_5, "MGSvar_type_5.csv", row.names=FALSE)
system('dx upload MGSvar_type_5.csv --path Emily-folder/results_2/MGSvar_type_5.csv')


# Neoplasm_mesothelium
LR_mesothelium <- glm(formula = Neoplasm_mesothelium ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSvar, 
                      data = EntireCohort_MGSvar, family = binomial)
# summary(LR_mesothelium)
# exp(cbind(OR = coef(LR_mesothelium), confint(LR_mesothelium)))
MGSvar_type_6 <- tidy(LR_mesothelium, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_type_6, "MGSvar_type_6.csv", row.names=FALSE)
system('dx upload MGSvar_type_6.csv --path Emily-folder/results_2/MGSvar_type_6.csv')


# Neoplasm_breast
LR_breast <- glm(formula = Neoplasm_breast ~ Age_at_recruitment + Ethnicity_White + 
                   Ever_smoked_dummy + Sex_Male + MGSvar, 
                 data = EntireCohort_MGSvar, family = binomial)
# summary(LR_breast)
# exp(cbind(OR = coef(LR_breast), confint(LR_breast)))
MGSvar_type_7 <- tidy(LR_breast, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_type_7, "MGSvar_type_7.csv", row.names=FALSE)
system('dx upload MGSvar_type_7.csv --path Emily-folder/results_2/MGSvar_type_7.csv')


# Neoplasm_femalegenital
LR_femalegenital <- glm(formula = Neoplasm_femalegenital ~ Age_at_recruitment + Ethnicity_White + 
                          Ever_smoked_dummy + Sex_Male + MGSvar, 
                        data = EntireCohort_MGSvar, family = binomial)
# summary(LR_femalegenital)
# exp(cbind(OR = coef(LR_femalegenital), confint(LR_femalegenital)))
MGSvar_type_8 <- tidy(LR_femalegenital, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_type_8, "MGSvar_type_8.csv", row.names=FALSE)
system('dx upload MGSvar_type_8.csv --path Emily-folder/results_2/MGSvar_type_8.csv')


# Neoplasm_malegenital
LR_malegenital <- glm(formula = Neoplasm_malegenital ~ Age_at_recruitment + Ethnicity_White + 
                        Ever_smoked_dummy + Sex_Male + MGSvar, 
                      data = EntireCohort_MGSvar, family = binomial)
# summary(LR_malegenital)
# exp(cbind(OR = coef(LR_malegenital), confint(LR_malegenital)))
MGSvar_type_9 <- tidy(LR_malegenital, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_type_9, "MGSvar_type_9.csv", row.names=FALSE)
system('dx upload MGSvar_type_9.csv --path Emily-folder/results_2/MGSvar_type_9.csv')


# Neoplasm_urinary
LR_urinary <- glm(formula = Neoplasm_urinary ~ Age_at_recruitment + Ethnicity_White + 
                    Ever_smoked_dummy + Sex_Male + MGSvar, 
                  data = EntireCohort_MGSvar, family = binomial)
# summary(LR_urinary)
# exp(cbind(OR = coef(LR_urinary), confint(LR_urinary)))
MGSvar_type_10 <- tidy(LR_urinary, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_type_10, "MGSvar_type_10.csv", row.names=FALSE)
system('dx upload MGSvar_type_10.csv --path Emily-folder/results_2/MGSvar_type_10.csv')


# Neoplasm_cns
LR_cns <- glm(formula = Neoplasm_cns ~ Age_at_recruitment + Ethnicity_White + 
                Ever_smoked_dummy + Sex_Male + MGSvar, 
              data = EntireCohort_MGSvar, family = binomial)
# summary(LR_cns)
# exp(cbind(OR = coef(LR_cns), confint(LR_cns)))
MGSvar_type_11 <- tidy(LR_cns, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_type_11, "MGSvar_type_11.csv", row.names=FALSE)
system('dx upload MGSvar_type_11.csv --path Emily-folder/results_2/MGSvar_type_11.csv')


# Neoplasm_endocrine
LR_endocrine <- glm(formula = Neoplasm_endocrine ~ Age_at_recruitment + Ethnicity_White + 
                      Ever_smoked_dummy + Sex_Male + MGSvar, 
                    data = EntireCohort_MGSvar, family = binomial)
# summary(LR_endocrine)
# exp(cbind(OR = coef(LR_endocrine), confint(LR_endocrine)))
MGSvar_type_12 <- tidy(LR_endocrine, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_type_12, "MGSvar_type_12.csv", row.names=FALSE)
system('dx upload MGSvar_type_12.csv --path Emily-folder/results_2/MGSvar_type_12.csv')


# Neoplasm_lymphatic
LR_lymphatic <- glm(formula = Neoplasm_lymphatic ~ Age_at_recruitment + Ethnicity_White + 
                      Ever_smoked_dummy + Sex_Male + MGSvar, 
                    data = EntireCohort_MGSvar, family = binomial)
# summary(LR_lymphatic)
# exp(cbind(OR = coef(LR_lymphatic), confint(LR_lymphatic)))
MGSvar_type_13 <- tidy(LR_lymphatic, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_type_13, "MGSvar_type_13.csv", row.names=FALSE)
system('dx upload MGSvar_type_13.csv --path Emily-folder/results_2/MGSvar_type_13.csv')


# Neoplasm_secondary
LR_secondary <- glm(formula = Neoplasm_secondary ~ Age_at_recruitment + Ethnicity_White + 
                      Ever_smoked_dummy + Sex_Male + MGSvar, 
                    data = EntireCohort_MGSvar, family = binomial)
# summary(LR_secondary)
# exp(cbind(OR = coef(LR_secondary), confint(LR_secondary)))
MGSvar_type_14 <- tidy(LR_secondary, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_type_14, "MGSvar_type_14.csv", row.names=FALSE)
system('dx upload MGSvar_type_14.csv --path Emily-folder/results_2/MGSvar_type_14.csv')


# Neoplasm_other
LR_other <- glm(formula = Neoplasm_other ~ Age_at_recruitment + Ethnicity_White + 
                  Ever_smoked_dummy + Sex_Male + MGSvar, 
                data = EntireCohort_MGSvar, family = binomial)
# summary(LR_other)
# exp(cbind(OR = coef(LR_other), confint(LR_other)))
MGSvar_type_15 <- tidy(LR_other, conf.int = TRUE, exp = TRUE)
write.csv(MGSvar_type_15, "MGSvar_type_15.csv", row.names=FALSE)
system('dx upload MGSvar_type_15.csv --path Emily-folder/results_2/MGSvar_type_15.csv')




































# MGSvar Age --------------------------------------------------------------




# a few checks before I start

summary(EntireCohort_MGSvar$Age_at_cancer_diagnosis_earliest)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 4.7    53.9    62.8    60.7    69.5    83.9  354948
469799 - 354948
# 114851 = size of dataset when NA observations are removed


# To check whether the response variable follows a normal distribution:
hist(EntireCohort_MGSvar$Age_at_cancer_diagnosis_earliest)
# nice and bell shaped - though skewed to the right (as expected)








## linear regression reference code:
# LM_ <- lm(Age_at_cancer_diagnosis_earliest ~ , data = EntireCohort_MGSvar)
# summary(LM_)

# "Age_at_recruitment"
LM_1 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment, data = EntireCohort_MGSvar)
# summary(LM_1)
ALL_age_1 <- tidy(LM_1, conf.int = TRUE)
write.csv(ALL_age_1, "ALL_age_1.csv", row.names=FALSE)
system('dx upload ALL_age_1.csv --path Emily-folder/results_2/ALL_age_1.csv')


# "BMI"
LM_2 <- lm(Age_at_cancer_diagnosis_earliest ~ BMI, data = EntireCohort_MGSvar)
# summary(LM_2)
ALL_age_2 <- tidy(LM_2, conf.int = TRUE)
write.csv(ALL_age_2, "ALL_age_2.csv", row.names=FALSE)
system('dx upload ALL_age_2.csv --path Emily-folder/results_2/ALL_age_2.csv')


# "Standing_height"                             
LM_3 <- lm(Age_at_cancer_diagnosis_earliest ~ Standing_height, data = EntireCohort_MGSvar)
# summary(LM_3)
ALL_age_3 <- tidy(LM_3, conf.int = TRUE)
write.csv(ALL_age_3, "ALL_age_3.csv", row.names=FALSE)
system('dx upload ALL_age_3.csv --path Emily-folder/results_2/ALL_age_3.csv')


# "MGSvar"                                      
LM_4 <- lm(Age_at_cancer_diagnosis_earliest ~ MGSvar, data = EntireCohort_MGSvar)
# summary(LM_4)
MGSvar_age_1 <- tidy(LM_4, conf.int = TRUE)
write.csv(MGSvar_age_1, "MGSvar_age_1.csv", row.names=FALSE)
system('dx upload MGSvar_age_1.csv --path Emily-folder/results_2/MGSvar_age_1.csv')


# "variant_load"                                
LM_5 <- lm(Age_at_cancer_diagnosis_earliest ~ variant_load, data = EntireCohort_MGSvar)
# summary(LM_5)
MGSvar_age_2 <- tidy(LM_5, conf.int = TRUE)
write.csv(MGSvar_age_2, "MGSvar_age_2.csv", row.names=FALSE)
system('dx upload MGSvar_age_2.csv --path Emily-folder/results_2/MGSvar_age_2.csv')


# "CDC45"                                       
# "CDT1"                                        
# "DONSON"                                      
# "GINS2"                                       
# "GINS3"                                       
# "MCM7"                                        
# "ORC1"                                        
# "ORC4"                                        
# "ORC6"                                        
LM_6 <- lm(Age_at_cancer_diagnosis_earliest ~ CDC45 + CDT1 + DONSON + GINS2 + GINS3 + MCM7 + ORC1 + ORC4 + ORC6, data = EntireCohort_MGSvar)
# summary(LM_6)
MGSvar_age_3 <- tidy(LM_6, conf.int = TRUE)
write.csv(MGSvar_age_3, "MGSvar_age_3.csv", row.names=FALSE)
system('dx upload MGSvar_age_3.csv --path Emily-folder/results_2/MGSvar_age_3.csv')

# won't save the individual gene ones - doesn't change compared to when they're added together
# "CDC45"                                           
LM_6_CDC45 <- lm(Age_at_cancer_diagnosis_earliest ~ CDC45, data = EntireCohort_MGSvar)
summary(LM_6_CDC45)
# "CDT1"                                           
LM_6_CDT1 <- lm(Age_at_cancer_diagnosis_earliest ~ CDT1, data = EntireCohort_MGSvar)
summary(LM_6_CDT1)
# "DONSON"                                        
LM_6_DONSON <- lm(Age_at_cancer_diagnosis_earliest ~ DONSON, data = EntireCohort_MGSvar)
summary(LM_6_DONSON)
# "GINS2"                                          
LM_6_GINS2 <- lm(Age_at_cancer_diagnosis_earliest ~ GINS2, data = EntireCohort_MGSvar)
summary(LM_6_GINS2)
# "GINS3"                                          
LM_6_GINS3 <- lm(Age_at_cancer_diagnosis_earliest ~ GINS3, data = EntireCohort_MGSvar)
summary(LM_6_GINS3)
# "MCM7"                                         
LM_6_MCM7 <- lm(Age_at_cancer_diagnosis_earliest ~ MCM7, data = EntireCohort_MGSvar)
summary(LM_6_MCM7)
# "ORC1"                                        
LM_6_ORC1 <- lm(Age_at_cancer_diagnosis_earliest ~ ORC1, data = EntireCohort_MGSvar)
summary(LM_6_ORC1)
# "ORC4"                                          
LM_6_ORC4 <- lm(Age_at_cancer_diagnosis_earliest ~ ORC4, data = EntireCohort_MGSvar)
summary(LM_6_ORC4)
# "ORC6"                                        
LM_6_ORC6 <- lm(Age_at_cancer_diagnosis_earliest ~ ORC6, data = EntireCohort_MGSvar)
summary(LM_6_ORC6)


# "Ethnicity_White"                             
LM_7 <- lm(Age_at_cancer_diagnosis_earliest ~ Ethnicity_White, data = EntireCohort_MGSvar)
# summary(LM_7)
ALL_age_4 <- tidy(LM_7, conf.int = TRUE)
write.csv(ALL_age_4, "ALL_age_4.csv", row.names=FALSE)
system('dx upload ALL_age_4.csv --path Emily-folder/results_2/ALL_age_4.csv')


# "Ever_smoked_dummy"                           
LM_8 <- lm(Age_at_cancer_diagnosis_earliest ~ Ever_smoked_dummy, data = EntireCohort_MGSvar)
# summary(LM_8)
ALL_age_5 <- tidy(LM_8, conf.int = TRUE)
write.csv(ALL_age_5, "ALL_age_5.csv", row.names=FALSE)
system('dx upload ALL_age_5.csv --path Emily-folder/results_2/ALL_age_5.csv')


# "Sex_Male"                                    
LM_9 <- lm(Age_at_cancer_diagnosis_earliest ~ Sex_Male, data = EntireCohort_MGSvar)
# summary(LM_9)
ALL_age_6 <- tidy(LM_9, conf.int = TRUE)
write.csv(ALL_age_6, "ALL_age_6.csv", row.names=FALSE)
system('dx upload ALL_age_6.csv --path Emily-folder/results_2/ALL_age_6.csv')


# "University_education"                        
LM_10 <- lm(Age_at_cancer_diagnosis_earliest ~ University_education, data = EntireCohort_MGSvar)
# summary(LM_10)
ALL_age_7 <- tidy(LM_10, conf.int = TRUE)
write.csv(ALL_age_7, "ALL_age_7.csv", row.names=FALSE)
system('dx upload ALL_age_7.csv --path Emily-folder/results_2/ALL_age_7.csv')


# "Annual_income"                               
LM_11 <- lm(Age_at_cancer_diagnosis_earliest ~ Annual_income, data = EntireCohort_MGSvar)
# summary(LM_11)
ALL_age_8 <- tidy(LM_11, conf.int = TRUE)
write.csv(ALL_age_8, "ALL_age_8.csv", row.names=FALSE)
system('dx upload ALL_age_8.csv --path Emily-folder/results_2/ALL_age_8.csv')


# "Weekly_exercise"
LM_12 <- lm(Age_at_cancer_diagnosis_earliest ~ Weekly_exercise, data = EntireCohort_MGSvar)
# summary(LM_12)
ALL_age_9 <- tidy(LM_12, conf.int = TRUE)
write.csv(ALL_age_9, "ALL_age_9.csv", row.names=FALSE)
system('dx upload ALL_age_9.csv --path Emily-folder/results_2/ALL_age_9.csv')






## multivariable analysis



# make dataset without NA rows

data_setD <- EntireCohort_MGSvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                    "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", 
                                    "MGSvar", "variant_load", "CDC45", "CDT1", "DONSON", 
                                    "GINS2", "GINS3", "MCM7", "ORC1", "ORC4", "ORC6")]
data_setD <- data_setD[complete.cases(data_setD), ]
dim(data_setD)
# 109936     17


with(data_setD, table(ORC6))
# ORC4 = 1 (removed)


data_setE <- EntireCohort_MGSvar[,c("eid", "Age_at_cancer_diagnosis_earliest", "Age_at_recruitment", 
                                    "Ethnicity_White", "Ever_smoked_dummy", "Sex_Male", "BMI",
                                    "Standing_height", "University_education", "Annual_income",
                                    "Weekly_exercise", "MGSvar", "variant_load", "CDC45", "CDT1", 
                                    "DONSON", "GINS2", "GINS3", "MCM7", "ORC1", "ORC4", "ORC6")]
data_setE <- data_setE[complete.cases(data_setE), ]
dim(data_setE)
# 88470    22


with(data_setE, table(ORC6))
# ORC4 = 0 (removed)




#data_setD

# (age, ethnicity, sex, smoking)
LM_13 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male, data = data_setD)
# summary(LM_13)
ALL_age_10 <- tidy(LM_13, conf.int = TRUE)
write.csv(ALL_age_10, "ALL_age_10.csv", row.names=FALSE)
system('dx upload ALL_age_10.csv --path Emily-folder/results_2/ALL_age_10.csv')



# MGSvar
LM_14 <- lm(Age_at_cancer_diagnosis_earliest ~ MGSvar, data = data_setD)
# summary(LM_14)
MGSvar_age_4 <- tidy(LM_14, conf.int = TRUE)
write.csv(MGSvar_age_4, "MGSvar_age_4.csv", row.names=FALSE)
system('dx upload MGSvar_age_4.csv --path Emily-folder/results_2/MGSvar_age_4.csv')


# (age, ethnicity, sex, smoking) + MGSvar
LM_15 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+MGSvar, data = data_setD)
# summary(LM_15)
# confint(LM_15)
MGSvar_age_5 <- tidy(LM_15, conf.int = TRUE)
write.csv(MGSvar_age_5, "MGSvar_age_5.csv", row.names=FALSE)
system('dx upload MGSvar_age_5.csv --path Emily-folder/results_2/MGSvar_age_5.csv')



#  variant_load
LM_16 <- lm(Age_at_cancer_diagnosis_earliest ~ variant_load, data = data_setD)
# summary(LM_16)
MGSvar_age_6 <- tidy(LM_16, conf.int = TRUE)
write.csv(MGSvar_age_6, "MGSvar_age_6.csv", row.names=FALSE)
system('dx upload MGSvar_age_6.csv --path Emily-folder/results_2/MGSvar_age_6.csv')


# (age, ethnicity, sex, smoking) + variant_load
LM_17 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+variant_load, data = data_setD)
# summary(LM_17)
MGSvar_age_7 <- tidy(LM_17, conf.int = TRUE)
write.csv(MGSvar_age_7, "MGSvar_age_7.csv", row.names=FALSE)
system('dx upload MGSvar_age_7.csv --path Emily-folder/results_2/MGSvar_age_7.csv')



#  GOIs
LM_18 <- lm(Age_at_cancer_diagnosis_earliest ~ CDC45 + CDT1 + DONSON + GINS2 + GINS3 + MCM7 + ORC1 + ORC6, data = data_setD)
# summary(LM_18)
MGSvar_age_8 <- tidy(LM_18, conf.int = TRUE)
write.csv(MGSvar_age_8, "MGSvar_age_8.csv", row.names=FALSE)
system('dx upload MGSvar_age_8.csv --path Emily-folder/results_2/MGSvar_age_8.csv')


# (age, ethnicity, sex, smoking) + GOIs
LM_19 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+CDC45+CDT1+DONSON+GINS2+GINS3+MCM7+ORC1+ORC6, data = data_setD)
# summary(LM_19)
# confint(LM_19)
MGSvar_age_9 <- tidy(LM_19, conf.int = TRUE)
write.csv(MGSvar_age_9, "MGSvar_age_9.csv", row.names=FALSE)
system('dx upload MGSvar_age_9.csv --path Emily-folder/results_2/MGSvar_age_9.csv')







# data_setE

# (age, ethnicity, sex, smoking, BMI, height, education, income, exercise)
LM_20 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise, data = data_setE)
# summary(LM_20)
ALL_age_11 <- tidy(LM_20, conf.int = TRUE)
write.csv(ALL_age_11, "ALL_age_11.csv", row.names=FALSE)
system('dx upload ALL_age_11.csv --path Emily-folder/results_2/ALL_age_11.csv')



# MGSvar
LM_21 <- lm(Age_at_cancer_diagnosis_earliest ~ MGSvar, data = data_setE)
# summary(LM_21)
MGSvar_age_10 <- tidy(LM_21, conf.int = TRUE)
write.csv(MGSvar_age_10, "MGSvar_age_10.csv", row.names=FALSE)
system('dx upload MGSvar_age_10.csv --path Emily-folder/results_2/MGSvar_age_10.csv')


# (age, ethnicity, sex, smoking, BMI, height, education, income, exercise) + MGSvar
LM_22 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+MGSvar, data = data_setE)
# summary(LM_22)
MGSvar_age_11 <- tidy(LM_22, conf.int = TRUE)
write.csv(MGSvar_age_11, "MGSvar_age_11.csv", row.names=FALSE)
system('dx upload MGSvar_age_11.csv --path Emily-folder/results_2/MGSvar_age_11.csv')



# variant_load
LM_23 <- lm(Age_at_cancer_diagnosis_earliest ~ variant_load, data = data_setE)
# summary(LM_23)
MGSvar_age_12 <- tidy(LM_23, conf.int = TRUE)
write.csv(MGSvar_age_12, "MGSvar_age_12.csv", row.names=FALSE)
system('dx upload MGSvar_age_12.csv --path Emily-folder/results_2/MGSvar_age_12.csv')


# (age, ethnicity, sex, smoking, BMI, height, education, income, exercise) + variant_load
LM_24 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+variant_load, data = data_setE)
# summary(LM_24)
MGSvar_age_13 <- tidy(LM_24, conf.int = TRUE)
write.csv(MGSvar_age_13, "MGSvar_age_13.csv", row.names=FALSE)
system('dx upload MGSvar_age_13.csv --path Emily-folder/results_2/MGSvar_age_13.csv')



#  GOIs
LM_25 <- lm(Age_at_cancer_diagnosis_earliest ~ CDC45 + CDT1 + DONSON + GINS2 + GINS3 + MCM7 + ORC1 + ORC6, data = data_setE)
# summary(LM_25)
MGSvar_age_14 <- tidy(LM_25, conf.int = TRUE)
write.csv(MGSvar_age_14, "MGSvar_age_14.csv", row.names=FALSE)
system('dx upload MGSvar_age_14.csv --path Emily-folder/results_2/MGSvar_age_14.csv')


# (age, ethnicity, sex, smoking, BMI, height, education, income, exercise) + GOIs
LM_26 <- lm(Age_at_cancer_diagnosis_earliest ~ Age_at_recruitment+Ethnicity_White+Ever_smoked_dummy+Sex_Male+BMI+Standing_height+University_education+Annual_income+Weekly_exercise+CDC45+CDT1+DONSON+GINS2+GINS3+MCM7+ORC1+ORC6, data = data_setE)
# summary(LM_26)
MGSvar_age_15 <- tidy(LM_26, conf.int = TRUE)
write.csv(MGSvar_age_15, "MGSvar_age_15.csv", row.names=FALSE)
system('dx upload MGSvar_age_15.csv --path Emily-folder/results_2/MGSvar_age_15.csv')













## fisher test for the linear regression data

# need to first subset by cancer diagnosis
EntireCohort_MGSvar <- EntireCohort_MGSvar[complete.cases(EntireCohort_MGSvar$Age_at_cancer_diagnosis_earliest), ]
# 110872


with(EntireCohort_MGSvar, table(MGSvar, Ethnicity_White))			
fisher.test(with(EntireCohort_MGSvar, table(MGSvar, Ethnicity_White)))
# p-value = 0.1614

with(EntireCohort_MGSvar, table(MGSvar, Ever_smoked_dummy))			
fisher.test(with(EntireCohort_MGSvar, table(MGSvar, Ever_smoked_dummy)))
# p-value = 0.03464
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.008920 1.325919
# sample estimates:
#   odds ratio 
# 1.155767

with(EntireCohort_MGSvar, table(MGSvar, Sex_Male))			
fisher.test(with(EntireCohort_MGSvar, table(MGSvar, Sex_Male)))
# p-value = 0.4762

with(EntireCohort_MGSvar, table(MGSvar, University_education))			
fisher.test(with(EntireCohort_MGSvar, table(MGSvar, University_education)))
# p-value = 0.9154

with(EntireCohort_MGSvar, table(MGSvar, Annual_income))			
fisher.test(with(EntireCohort_MGSvar, table(MGSvar, Annual_income)))
# p-value = 0.4582


### data_setD
with(data_setD, table(MGSvar, Ethnicity_White))			
fisher.test(with(data_setD, table(MGSvar, Ethnicity_White)))
# p-value = 0.1587

with(data_setD, table(MGSvar, Ever_smoked_dummy))			
fisher.test(with(data_setD, table(MGSvar, Ever_smoked_dummy)))
# p-value = 0.02903
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.013914 1.333363
# sample estimates:
#   odds ratio 
# 1.161834

with(data_setD, table(MGSvar, Sex_Male))			
fisher.test(with(data_setD, table(MGSvar, Sex_Male)))
# p-value = 0.537


### data_setE
with(data_setE, table(MGSvar, Ethnicity_White))			
fisher.test(with(data_setE, table(MGSvar, Ethnicity_White)))
# p-value = 1

with(data_setE, table(MGSvar, Ever_smoked_dummy))			
fisher.test(with(data_setE, table(MGSvar, Ever_smoked_dummy)))
# p-value = 0.009113
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.051171 1.435849
# sample estimates:
#   odds ratio 
# 1.22725

with(data_setE, table(MGSvar, Sex_Male))			
fisher.test(with(data_setE, table(MGSvar, Sex_Male)))
# p-value = 0.5856

with(data_setE, table(MGSvar, University_education))			
fisher.test(with(data_setE, table(MGSvar, University_education)))
# p-value = 0.9386

with(data_setE, table(MGSvar, Annual_income))			
fisher.test(with(data_setE, table(MGSvar, Annual_income)))
# p-value = 0.3757



## make plots for the linear regression data

# starting with the EntireCohort_MGSvar dataset subset by age at cancer diagnosis 
EntireCohort_MGSvar <- EntireCohort_MGSvar[complete.cases(EntireCohort_MGSvar$Age_at_cancer_diagnosis_earliest), ]
# 110872

# make MGSvar into factor
EntireCohort_MGSvar$MGSvar <- as.factor(EntireCohort_MGSvar$MGSvar)

# try a violin plot
# Age_at_cancer_diagnosis_earliest
EntireCohort_MGSvar %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

EntireCohort_MGSvar %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ MGSvar, data = EntireCohort_MGSvar)
# 



# Age_at_recruitment  
EntireCohort_MGSvar %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

EntireCohort_MGSvar %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ MGSvar, data = EntireCohort_MGSvar)
# W = 53139946, p-value = 0.9132



# BMI  
EntireCohort_MGSvar %>% 
  ggplot(mapping = aes(x = MGSvar, y = BMI)) +
  geom_violin() +
  theme_light()

EntireCohort_MGSvar %>% 
  ggplot(mapping = aes(x = MGSvar, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ MGSvar, data = EntireCohort_MGSvar)
# W = 53542684, p-value = 0.5424



# Standing_height 
EntireCohort_MGSvar %>% 
  ggplot(mapping = aes(x = MGSvar, y = Standing_height)) +
  geom_violin() +
  theme_light()

EntireCohort_MGSvar %>% 
  ggplot(mapping = aes(x = MGSvar, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ MGSvar, data = EntireCohort_MGSvar)
# W = 53297362, p-value = 0.7659



# make this a bar plot
# Weekly_exercise
EntireCohort_MGSvar %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(MGSvar), scales = "free_y") +
  theme_light()
# 




### data_setD

data_setD$MGSvar <- as.factor(data_setD$MGSvar)

# Age_at_cancer_diagnosis_earliest
data_setD %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ MGSvar, data = data_setD)
# 



# Age_at_recruitment
data_setD %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setD %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ MGSvar, data = data_setD)
# W = 52261828, p-value = 0.8744




### data_setE

data_setE$MGSvar <- as.factor(data_setE$MGSvar)

# Age_at_cancer_diagnosis_earliest
data_setE %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_cancer_diagnosis_earliest)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_cancer_diagnosis_earliest)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_cancer_diagnosis_earliest ~ MGSvar, data = data_setE)
# 



# Age_at_recruitment
data_setE %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_recruitment)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = MGSvar, y = Age_at_recruitment)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Age_at_recruitment ~ MGSvar, data = data_setE)
# W = 33176761, p-value = 0.6415



# BMI
data_setE %>% 
  ggplot(mapping = aes(x = MGSvar, y = BMI)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = MGSvar, y = BMI)) +
  geom_boxplot() +
  theme_light()

wilcox.test(BMI ~ MGSvar, data = data_setE)
# W = 34243668, p-value = 0.2924



# Standing_height
data_setE %>% 
  ggplot(mapping = aes(x = MGSvar, y = Standing_height)) +
  geom_violin() +
  theme_light()

data_setE %>% 
  ggplot(mapping = aes(x = MGSvar, y = Standing_height)) +
  geom_boxplot() +
  theme_light()

wilcox.test(Standing_height ~ MGSvar, data = data_setE)
# W = 33764530, p-value = 0.7104


# make this a bar plot
# Weekly_exercise
data_setE %>% 
  ggplot(mapping = aes(x = Weekly_exercise)) +
  geom_bar() +
  facet_grid(rows = vars(MGSvar), scales = "free_y") +
  theme_light()
# 









# Save the script!!! -----------------------------------------------------------

system('dx upload EntireCohort_MGSvar_analysis.R --path Emily-folder/MGSvar/EntireCohort_MGSvar_analysis.R')

