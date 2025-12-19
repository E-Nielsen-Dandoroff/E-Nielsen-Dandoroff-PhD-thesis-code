## EntireCohort_script.R

# This script intends to take the Emily-folder/EntireCohort_Cancer_participant_2.csv 
# dataset and prepare it for statistical analysis.

# The code within this script was originally within:
# 'Emily-folder/MGSvar/EntireCohort_MGSvar_script2.R' (EntireCohort_Cancer_participant_2.csv -> EntireCohort_MGS.csv)
# 'Emily-folder/MGSvar/EntireCohort_MGS_analysis.R' (EntireCohort_MGS.csv -> EntireCohort_MGS_expanded_ICD.csv)

# new code added for cancer-related data (age at diagnosis and ICD categories)
# from script saved to my personal computer ('Backup copies of scripts/instance_data_investigation.R')




####
# This script only requires execution once
# last updated & executed: 12/09/2024
####





# Load packages and dataframe --------------------------------------------------

install.packages('readr')
library(readr)
library(data.table)
library(dplyr)
library(ggplot2)
install.packages('reshape2')
library(reshape2)
install.packages('ggmosaic')
library(ggmosaic)

# library(tidyr)
## load these in if required - if not remove

system('dx download Emily-folder/EntireCohort_Cancer_participant_2.csv')
EntireCohort <- read_csv("EntireCohort_Cancer_participant_2.csv")



# Adjust and clean EntireCohort dataframe --------------------------------------

# rename datafields
EntireCohort <- rename(EntireCohort, 
                       Age_at_recruitment = p21022,
                       Ethnic_background = p21000_i0,
                       Ever_smoked = p20160_i0,
                       Smoking_status = p20116_i0,
                       Genetic_relatedness_pairing = p22011_a0,
                       Year_of_birth = p34,
                       Reported_occurrencces_of_cancer = p40009_i0,
                       Age_at_cancer_diagnosis = p40008_i0,
                       Type_of_cancer_ICD10 = p40006_i0,
                       Type_of_cancer_ICD9 = p40013_i0,
                       Sex = p31,
                       BMI = p21001_i0,
                       Qualifications = p6138_i0,
                       Average_household_income = p738_i0,
                       Number_days_per_week_10min_physical_activity = p884_i0,
                       Standing_height = p50_i0)

# go through each data field and check/adjust

## eid
# leave!


## Age at recruitment  
summary(EntireCohort$Age_at_recruitment)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   37.00   50.00   58.00   56.53   63.00   73.00       1 
# leave!


## Ethnic Background
# this one is a little tricky, there are so many categories!
# lets have a look at the categories UKB has:
unique(EntireCohort$Ethnic_background)
# [1] "British"                    "Irish"                      "Any other white background"
# [4] "Caribbean"                  "Any other Asian background" "Other ethnic group"        
# [7] "Indian"                     "White"                      "African"                   
# [10] "Chinese"                    NA                           "Prefer not to answer"      
# [13] "Pakistani"                  "Any other mixed background" "White and Asian"           
# [16] "Any other Black background" "White and Black Caribbean"  "White and Black African"   
# [19] "Mixed"                      "Bangladeshi"                "Asian or Asian British"    
# [22] "Do not know"                "Black or Black British"

# other researchers have categorised as 'white' and 'non-white' (using dummy variables!)
# Ethnicity_White = "British", "Irish", "Any other white background", "White"
# 0 values in this dummy variable will be all others (other than NA and "Prefer not to answer")
# I will make "Prefer not to answer" into NA (NA will remain as missing data)
EntireCohort$Ethnicity_White <- EntireCohort$Ethnic_background
EntireCohort$Ethnicity_White[EntireCohort$Ethnicity_White == "Prefer not to answer"] <- NA
# check
unique(EntireCohort$Ethnicity_White)
# I will make "Do not know" into NA (NA will remain as missing data)
EntireCohort$Ethnicity_White[EntireCohort$Ethnicity_White == "Do not know"] <- NA

# now to make "British", "Irish", "Any other white background", "White" into '1'
EntireCohort$Ethnicity_White[EntireCohort$Ethnicity_White == "British"] <- 1
EntireCohort$Ethnicity_White[EntireCohort$Ethnicity_White == "Irish"] <- 1
EntireCohort$Ethnicity_White[EntireCohort$Ethnicity_White == "Any other white background"] <- 1
EntireCohort$Ethnicity_White[EntireCohort$Ethnicity_White == "White"] <- 1
# now to make all others (except for NA) into '0'
EntireCohort$Ethnicity_White[EntireCohort$Ethnicity_White != "1"] <- 0
# now to change to int
EntireCohort$Ethnicity_White <- as.integer(EntireCohort$Ethnicity_White)
summary(EntireCohort$Ethnicity_White)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.0000  1.0000  1.0000  0.9459  1.0000  1.0000    2777
# there is a huge proportion of white people in this dataset!

# Ethnicity data is now ready to use in statistical analysis
# Make sure not to use the 'Ethnic_background' data field!!!


## Ever_smoked
# make into dummy variables
unique(EntireCohort$Ever_smoked)
# "Yes" "No"  NA 
# assign 1=Yes, 0=No
EntireCohort <- EntireCohort %>% 
  mutate(Ever_smoked_dummy = if_else(Ever_smoked == "Yes", 1, 0))
# Note: use this new column for statistical analysis and not 'Ever_smoked'!!!
unique(EntireCohort$Ever_smoked_dummy)
# 1  0 NA


## Smoking_status
# make into dummy variables
unique(EntireCohort$Smoking_status)
# "Previous"             "Current"              "Never"                NA                     "Prefer not to answer"
# missing data = NA and "Prefer not to answer" (make into NA)
# Null = "Never"
# Smoker_previous = "Previous"
# Smoker_current = "Current"

# Smoker_previous:
EntireCohort$Smoker_previous <- EntireCohort$Smoking_status
EntireCohort$Smoker_previous[EntireCohort$Smoker_previous=="Previous"] <- 1
EntireCohort$Smoker_previous[EntireCohort$Smoker_previous=="Current"] <- 0
EntireCohort$Smoker_previous[EntireCohort$Smoker_previous=="Never"] <- 0
# change "Prefer not to answer" into NA
EntireCohort$Smoker_previous[EntireCohort$Smoker_previous=="Prefer not to answer"] <- NA
# make into int
EntireCohort$Smoker_previous <- as.integer(EntireCohort$Smoker_previous)
# check
unique(EntireCohort$Smoker_previous)
# 1  0 NA

# Smoker_current:
EntireCohort$Smoker_current <- EntireCohort$Smoking_status
EntireCohort$Smoker_current[EntireCohort$Smoker_current=="Previous"] <- 0
EntireCohort$Smoker_current[EntireCohort$Smoker_current=="Current"] <- 1
EntireCohort$Smoker_current[EntireCohort$Smoker_current=="Never"] <- 0
# change "Prefer not to answer" into NA
EntireCohort$Smoker_current[EntireCohort$Smoker_current=="Prefer not to answer"] <- NA
# make into int
EntireCohort$Smoker_current <- as.integer(EntireCohort$Smoker_current)
# check
unique(EntireCohort$Smoker_current)
# 0  1 NA
# Note: use these new columns for statistical analysis and not 'Smoking_status'!!!



## Reported occurances of cancer
# want to make dummy variables out of this for cancer y/n
# this field doesn't have multiple instances - which is very helpfull
EntireCohort <- EntireCohort %>% 
  mutate(cancer_ROOC_dummy = if_else(Reported_occurrencces_of_cancer>0, 1, 0))
# check
summary(EntireCohort$cancer_ROOC_dummy)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#       1       1       1       1       1       1  383415 

# # information about this datafield from UKB:
# This is the number of reports of cancer received from UK National registries. Participants who have never been reported as having cancer do not contribute to this field unless a report was received then subsequently withdrawn.
# Please note that this field represents simply the number of cancer records per participant; due to probable (pseudo)duplicates between different sources these counts likely do not accurately reflect the number of unique cancer diagnoses. Please see Resource 115558 for further details.
# The counts indicate records where there is a non-missing date of diagnosis; a record will be counted even where other diagnosis information such as site and/or morphology are missing or invalid.
# # this was the most relevant cancer registry data field for my purposes

# even though there are no true NA results - we must assume that all NA = 0 reported cancer incidences
EntireCohort$cancer_ROOC_dummy[is.na(EntireCohort$cancer_ROOC_dummy)] <- 0
summary(EntireCohort$cancer_ROOC_dummy)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.0000  0.2368  0.0000  1.0000
# Can use this new field as an indication of cancer y/n
# only use this new datafield in the statistical analyses! don't use 'reported occurances of cancer'!!!



## sex
# also want to make this into a dummy variable
unique(EntireCohort$Sex)
# "Female" "Male" 
# make a dummy field for this with 0=female and 1=male
EntireCohort <- EntireCohort %>% 
  mutate(Sex_Male = if_else(Sex == "Male", 1, 0))
# check
summary(EntireCohort$Sex_Male)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   0.000   0.456   1.000   1.000
# so use this new field instead of 'Sex' when it comes to analysis!!!


## BMI
summary(EntireCohort$BMI)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   12.12   24.14   26.74   27.43   29.91   74.68    3101
# there are quite a few NA values, but not so many that the field isn't useful.
# a quick google let me know that a BMI of 12.12 isn't healthy but is possible,
# and 74.68 seems like a legitimate value as well.
# I will leave these as they are.


## Qualifications
# this is a categorical field
unique(EntireCohort$Qualifications)
# there are 66 different categories, mainly because it seems to be a survey question
# that can have multiple responses selected.
# one is NA

# The UKB's official information for this field:
# 1	College or University degree
# 2	A levels/AS levels or equivalent
# 3	O levels/GCSEs or equivalent
# 4	CSEs or equivalent
# 5	NVQ or HND or HNC or equivalent
# 6	Other professional qualifications eg: nursing, teaching
# -7	None of the above
# -3	Prefer not to answer

# previous publications using UKB data have grouped by 'university' and 'non-university'
# want NA and "Prefer not to answer" to be missing data (NA)
# University_education = "College or University degree" = 1, and all other options = 0
EntireCohort$University_education <- EntireCohort$Qualifications
# make "Prefer not to answer" into NA
EntireCohort$University_education[EntireCohort$University_education=="Prefer not to answer"] <- NA
# make "College or University degree" = 1
EntireCohort$University_education[grep("College or University degree", EntireCohort$University_education)] <- 1
EntireCohort$University_education[EntireCohort$University_education != "1"] <- 0
# change to int
EntireCohort$University_education <- as.integer(EntireCohort$University_education)
# check
unique(EntireCohort$University_education)
# 1  0 NA
summary(EntireCohort$University_education)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   0.000   0.000   0.327   1.000   1.000   10130
## make sure to use 'University_education' and not 'Qualifications' in analysis !!!


## average household income
# this will also be categorical
unique(EntireCohort$Average_household_income)
# "52,000 to 100,000"    "Prefer not to answer" "18,000 to 30,999"     "Do not know"         
# "31,000 to 51,999"     "Less than 18,000"     "Greater than 100,000" NA

# published papers have further categorised, one set a cutoff at 52,000, which seems reasonable
# missing data = NA, "Prefer not to answer", and "Do not know" (make into NA)
# Annual_income = "Greater than 100,000" and "52,000 to 100,000" = 1
# incomes lower than 52,000 will be the reference = 0
EntireCohort$Annual_income <- EntireCohort$Average_household_income
# change "Do not know" and "Prefer not to answer" into NA
EntireCohort$Annual_income[EntireCohort$Annual_income=="Do not know"] <- NA
EntireCohort$Annual_income[EntireCohort$Annual_income=="Prefer not to answer"] <- NA
# Annual_income = "Greater than 100,000" and "52,000 to 100,000" = 1
EntireCohort$Annual_income[EntireCohort$Annual_income=="Greater than 100,000"] <- 1
EntireCohort$Annual_income[EntireCohort$Annual_income=="52,000 to 100,000"] <- 1
# others = 0
EntireCohort$Annual_income[EntireCohort$Annual_income != "1"] <- 0
# change to int
EntireCohort$Annual_income <- as.integer(EntireCohort$Annual_income)
#check
unique(EntireCohort$Annual_income)
# 1 NA  0
summary(EntireCohort$Annual_income)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    0.00    0.00    0.00    0.26    1.00    1.00   77140
# make sure to use the new field and not 'Average_household_income' !!!


## number of days per week 10min physical activity
unique(EntireCohort$Number_days_per_week_10min_physical_activity)
# [1] "6"                    "4"                    "5"                    "2"                   
# [5] "0"                    "Do not know"          "7"                    "3"                   
# [9] "1"                    "Prefer not to answer" NA
# slightly alter this to remove the strings and convert to numerical
# convert "Do not know" and "Prefer not to answer" to NA values
EntireCohort$Weekly_exercise <- EntireCohort$Number_days_per_week_10min_physical_activity
EntireCohort$Weekly_exercise[EntireCohort$Weekly_exercise == "Do not know"] <- NA
EntireCohort$Weekly_exercise[EntireCohort$Weekly_exercise == "Prefer not to answer"] <- NA
# change to numerical
unique(EntireCohort$Weekly_exercise)
EntireCohort$Weekly_exercise <- as.numeric(EntireCohort$Weekly_exercise)
# make sure to use this field instead of 'Number_days_per_week_10min_physical_activity'!!!


## standing height
summary(EntireCohort$Standing_height)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    75.0   161.8   168.0   168.4   175.0   209.0    2535
sort(EntireCohort$Standing_height)
# 75.0 100.0 100.0 118.0 121.0 122.0 122.0 122.0 122.0 125.0 125.0 126.0 126.0 127.0 127.0 128.0 128.0 128.0
# Since this is a measurement that UKB takes themselves I'm just going to assume that this is correct.
# There are people that small, although a bit suspicious.
# going to leave this data field as it is.











# Remove participants without genotype data ------------------------------------


# also need to remove the participants that do not have genotype information.
# within the 'Emily-folder/pathvar/pathvar_script_END.R' script I created a list of 
# IDs of people that have been genotyped, now I just need to remove the rows with IDS that do not match.

# load in the list of IDs to keep
participant_IDs <- fread('/mnt/project/Emily-folder/genotyped_participants.csv', data.table = FALSE)
participant_IDs <- participant_IDs[-1]

# use inner_join() to keep only the IDs that match between the two dfs
participant_IDs <- rename(participant_IDs, eid = x)

EntireCohort_ready <- inner_join(participant_IDs, EntireCohort, by = 'eid')
# EntireCohort had 502394 observations
# participant_IDs had 469835 observations
# EntireCohort_ready has 469799 observations

# I want to check out the observations that have been excluded in this new df
dif1 <- setdiff(EntireCohort$eid, participant_IDs$eid)
502394 - 32595
# 469799

dif1 <- setdiff(participant_IDs$eid, EntireCohort$eid)
# negative numbers 1:36
# these may be people who have retracted data ???
# either way, they do not match with the EntireCohort data, so shouldn't be included.



# all done with the data adjusting and cleaning :)
# I will save this dataframe as .csv so that its ready to go for the subsequent analyses.

#save it to R
write.csv(EntireCohort_ready, "EntireCohort_ready.csv")
# now to save this to UKB project
system('dx upload EntireCohort_ready.csv --path Emily-folder/EntireCohort_ready.csv')

## load this EntireCohort_ready.csv dataset in when starting a new analysis











# age at cancer diagnosis ----------------------------------------

## adding age at first cancer diagnosis information.
# script initially included code with just the instance 0 data, but have since realised that 
# instances are not ordered chronologically.
# the following code determines the youngest reported age at cancer diagnosis for each UKB individual.

system('dx download Emily-folder/EntireCohort_ready.csv')
EntireCohort_ready <- read_csv("EntireCohort_ready.csv")
EntireCohort_ready <- EntireCohort_ready[,-1]


# I have made a CSV file with all of the age at cancer diagnosis instances:
system('dx download Emily-folder/Age_at_cancer_diagnosis_instances.csv')

# load in the data - and make sure that all fields are loaded as numbers and not logical!
Age_at_cancer_diagnosis_instances <- read_csv("Age_at_cancer_diagnosis_instances.csv", 
                                              col_types = cols(p40008_i4 = col_double(), 
                                                               p40008_i5 = col_double(), p40008_i6 = col_double(), 
                                                               p40008_i7 = col_double(), p40008_i8 = col_double(), 
                                                               p40008_i9 = col_double(), p40008_i10 = col_double(), 
                                                               p40008_i11 = col_double(), p40008_i12 = col_double(), 
                                                               p40008_i13 = col_double(), p40008_i14 = col_double(), 
                                                               p40008_i15 = col_double(), p40008_i16 = col_double(), 
                                                               p40008_i17 = col_double(), p40008_i18 = col_double(), 
                                                               p40008_i19 = col_double(), p40008_i20 = col_double(), 
                                                               p40008_i21 = col_double()))

# first: clean up the table a little by removing any row without age at cancer diagnosis info
# there's probably a cleaner way to do this - but the easiest way I can think of would be to sum the columns
# and remove rows that have sum NA
Age_at_cancer_diagnosis_instances$Sums <- rowSums(Age_at_cancer_diagnosis_instances[ , 2:23], na.rm = TRUE)
# remove rows that have 0 score for Sums
Age_at_cancer_diagnosis_instances <- Age_at_cancer_diagnosis_instances[Age_at_cancer_diagnosis_instances$Sums !=0, ]
# reduced cohort size from 502364 to 123776

# shouldn't remove any of the instances - but can get rid of sum column now
Age_at_cancer_diagnosis_instances <- Age_at_cancer_diagnosis_instances[,1:23]

# change NA to an obscenely high number
Age_at_cancer_diagnosis_instances[is.na(Age_at_cancer_diagnosis_instances)] <- 999



Age_at_cancer_diagnosis_instances$Age_at_cancer_diagnosis_earliest <- pmin(Age_at_cancer_diagnosis_instances$p40008_i0, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i1, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i2, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i3, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i4, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i5, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i6, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i7, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i8, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i9, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i10, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i11, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i12, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i13, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i14, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i15, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i16, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i17, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i18, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i19, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i20, 
                                                                           Age_at_cancer_diagnosis_instances$p40008_i21)

summary(Age_at_cancer_diagnosis_instances$Age_at_cancer_diagnosis_earliest)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.20   53.90   62.80   60.67   69.50   84.50 

# keep just the youngest cancer diagnosis ages and eid
Age_at_cancer_diagnosis_earliest <- Age_at_cancer_diagnosis_instances[,c(1,24)]

summary(Age_at_cancer_diagnosis_earliest$Age_at_cancer_diagnosis_earliest)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.20   53.90   62.80   60.67   69.50   84.50


total <- left_join(EntireCohort_ready, Age_at_cancer_diagnosis_earliest, by="eid")


#save it to R
write.csv(total, "EntireCohort_ready.csv")
# now to save this to UKB project
system('dx upload EntireCohort_ready.csv --path Emily-folder/EntireCohort_ready.csv')














# cancer type, including all ----------------------------------------------

## adding neoplasm type information.
# had been previously using only one cancer diagnsis per individual, but I have discovered that these 
# instances are not in chronological order, so it is not possible to determine first neoplasm diagnosis.
# the following code adds neoplasm type fields for each of the neoplasm types used in this research.


# I have made a CSV file with all of the type of cancer ICD codes instances:
system('dx download Emily-folder/Type_of_cancer_instances.csv')

# load data (need to make sure that the ICD codes get loaded in as characters (data goes missing when a column defaults to logical))
Type_of_cancer_instances <- read_csv("Type_of_cancer_instances.csv", 
                                     col_types = cols(p40006_i2 = col_character(), 
                                                      p40006_i4 = col_character(), p40006_i5 = col_character(), 
                                                      p40006_i6 = col_character(), p40006_i7 = col_character(), 
                                                      p40006_i8 = col_character(), p40006_i9 = col_character(), 
                                                      p40006_i10 = col_character(), p40006_i11 = col_character(), 
                                                      p40006_i12 = col_character(), p40006_i13 = col_character(), 
                                                      p40006_i14 = col_character(), p40006_i15 = col_character(), 
                                                      p40006_i16 = col_character(), p40006_i17 = col_character(), 
                                                      p40006_i18 = col_character(), p40006_i19 = col_character(), 
                                                      p40006_i20 = col_character(), p40006_i21 = col_character(), 
                                                      p40013_i3 = col_character(), p40013_i4 = col_character(), 
                                                      p40013_i5 = col_character(), p40013_i6 = col_character(), 
                                                      p40013_i7 = col_character(), p40013_i8 = col_character(), 
                                                      p40013_i9 = col_character(), p40013_i10 = col_character(), 
                                                      p40013_i11 = col_character(), p40013_i12 = col_character(), 
                                                      p40013_i13 = col_character(), p40013_i14 = col_character()))

# check whether any instances lack any data
unique(Type_of_cancer_instances$p40006_i0)
unique(Type_of_cancer_instances$p40006_i1)
unique(Type_of_cancer_instances$p40006_i2)
unique(Type_of_cancer_instances$p40006_i3)
unique(Type_of_cancer_instances$p40006_i4)
unique(Type_of_cancer_instances$p40006_i5)
unique(Type_of_cancer_instances$p40006_i6)
unique(Type_of_cancer_instances$p40006_i7)
unique(Type_of_cancer_instances$p40006_i8)
unique(Type_of_cancer_instances$p40006_i9)
unique(Type_of_cancer_instances$p40006_i10)
unique(Type_of_cancer_instances$p40006_i11)
unique(Type_of_cancer_instances$p40006_i12)
unique(Type_of_cancer_instances$p40006_i13)
unique(Type_of_cancer_instances$p40006_i14)
unique(Type_of_cancer_instances$p40006_i15)
unique(Type_of_cancer_instances$p40006_i16)
unique(Type_of_cancer_instances$p40006_i17)
unique(Type_of_cancer_instances$p40006_i18)
unique(Type_of_cancer_instances$p40006_i19)
unique(Type_of_cancer_instances$p40006_i20)
unique(Type_of_cancer_instances$p40006_i21)
unique(Type_of_cancer_instances$p40013_i0)
unique(Type_of_cancer_instances$p40013_i1)
unique(Type_of_cancer_instances$p40013_i2)
unique(Type_of_cancer_instances$p40013_i3)
unique(Type_of_cancer_instances$p40013_i4)
unique(Type_of_cancer_instances$p40013_i5)
unique(Type_of_cancer_instances$p40013_i6)
unique(Type_of_cancer_instances$p40013_i7)
unique(Type_of_cancer_instances$p40013_i8)
unique(Type_of_cancer_instances$p40013_i10)
unique(Type_of_cancer_instances$p40013_i11)
unique(Type_of_cancer_instances$p40013_i12)
unique(Type_of_cancer_instances$p40013_i14)
# all contain data except for the following two:
unique(Type_of_cancer_instances$p40013_i9)
unique(Type_of_cancer_instances$p40013_i13)

# remove the columns that lack any data
Type_of_cancer_instances <- Type_of_cancer_instances[-c(33,37)]

# first: clean up the table a little by removing any row without cancer code info
# can't do what I did for age (because you can't sum strings)
# can do sum of NA values - and remove any row that has NA_sums = 35
Type_of_cancer_instances$NA_sums <- rowSums(is.na(Type_of_cancer_instances))
# have a look at summary stats
summary(Type_of_cancer_instances$NA_sums)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 13.00   35.00   35.00   34.67   35.00   35.00

# remove rows that only have NA values for ICD code instances
Type_of_cancer_instances <- subset(Type_of_cancer_instances, NA_sums < 35)
# from 502364 to 123772 obs





# trying to figure out a way to replace the ICD codes with the category names without taking too much code repetition!
# have figured out a way:
# make vectors for ICD10 and ICD9 codes to find and replace
# then make vectors for the column names to search within
# use these vectors to search through the df and replace with whichever string category is required!

Type_of_cancer_instances_categories <- Type_of_cancer_instances[1:36]

# make vectors for ICD10 and ICD9 codes to find and replace
### oral
ICD10_oral <- c("C00", "C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12", "C13",
                "C14", "D00.0", "D10", "D11", "D37.0")
ICD9_oral <- c("1401", "1409", "1410", "1412", "1413", "1419", "1420", "1421", "1428", "1429", "1431", "1449",
               "1450", "1451", "1452", "1453", "1455", "1459", "1460", "1461", "1470", "1479", "2300", "2350", "2351")
### digestive
ICD10_digestive <- c("C15", "C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "D00.1",
                     "D00.2", "D01", "D12", "D13", "D37.1", "D37.2", "D37.3", "D37.4", "D37.5", "D37.6", "D37.7", "D37.9" )
ICD9_digestive <- c("1505", "1509", "1510", "1512", "1514", "1515", "1519", "1520", "1521", "1522", "1529", "1530",
                    "1531", "1532", "1533", "1534", "1535", "1536", "1537", "1538", "1539", "1540", "1541", "1542",
                    "1543", "1548", "1550", "1561", "1562", "1570", "1572", "1574", "1579", "2303", "2304", "2352",
                    "2353", "2355", "2390", "150 ")
### respiratory
ICD10_respiratory <- c("C30", "C31", "C32", "C33", "C34", "C37", "C38", "C39", "D02", "D15", "D38")
ICD9_respiratory <- c("1600", "1601", "1602", "1610", "1611", "1613", "1619", "1620", "1622", "1623", "1624", "1625",
                      "1628", "1629", "1639", "1640", "1642", "1649", "2310", "2357", "162 ")
### bone
ICD10_bone <- c("C40", "C41", "D16", "D48.0")
ICD9_bone <- c("1700", "1702", "1703", "1704", "1705", "1706", "1707", "1708", "1709", "2380")
### skin
ICD10_skin <- c("C43", "C44", "D03", "D04", "D48.5")
ICD9_skin <- c("1720", "1721", "1722", "1723", "1724", "1725", "1726", "1727", "1728", "1729", "1730", "1731", "1732",
               "1733", "1734", "1735", "1736", "1737", "1738", "1739", "2321", "2322", "2323", "2324", "2325", "2326",
               "2327", "2328", "2329", "2382")
### mesothelium
ICD10_mesothelium <- c("C45", "C46", "C47", "C48", "C49", "D21", "D48.1", "D48.2", "D48.3", "D48.4")
ICD9_mesothelium <- c("1580", "1599", "1710", "1712", "1713", "1714", "1715", "1716", "1717", "1718", "1719", "2354",
                      "2381", "171 ")
### breast
ICD10_breast <- c("C50", "D05", "D48.6")
ICD9_breast <- c("1740", "1741", "1742", "1743", "1744", "1745", "1746", "1748", "1749", "2330", "2383", "175 ")
### femalegenital
ICD10_femalegenital <- c("C51", "C52", "C53", "C54", "C55", "C56", "C57", "C58", "D06", "D07.0", "D07.1", "D07.2",
                         "D07.3", "D27", "D39")
ICD9_femalegenital <- c("1800", "1801", "1808", "1809", "1820", "1828", "1830", "1832", "1840", "1841", "1844", "1848",
                        "1849", "2331", "2332", "2333", "2360", "2361", "2362", "2363", "181 ", "182 ", "179 ")
### malegenital
ICD10_malegenital <- c("C60", "C61", "C62", "C63", "D07.4", "D07.5", "D07.6", "D29", "D40")
ICD9_malegenital <- c("1869", "1871", "1874", "1875", "1877", "2335", "2336", "2364", "2365", "185 ", "187 ")
### urinary
ICD10_urinary <- c("C64", "C65", "C66", "C67", "C68", "D09.0", "D09.1", "D30.", "D41")
ICD9_urinary <- c("1880", "1882", "1884", "1886", "1888", "1889", "1890", "1891", "1892", "1893", "1899", "2233", "2337",
                  "2339", "2367", "2369")
### cns
ICD10_cns <- c("C69", "C70", "C71", "C72", "D09.2", "D32", "D33", "D42", "D43")
ICD9_cns <- c("1900", "1901", "1903", "1906", "1908", "1909", "1910", "1911", "1912", "1914", "1916", "1917", "1918",
              "1919", "1920", "1921", "1922", "2250", "2251", "2252", "2253", "2254", "2258", "2259", "2340", "2375",
              "2376", "2377", "2379", "2396", "225 ", "192 ")
### endocrine
ICD10_endocrine <- c("C73", "C74", "C75", "D09.3", "D34", "D35", "D44")
ICD9_endocrine <- c("1940", "1941", "1943", "1944", "1949", "2273", "2279", "2370", "2373", "193 ", "226 ")
### lymphatic
ICD10_lymphatic <- c("C42", "C81", "C82", "C83", "C84", "C85", "C86", "C88", "C90", "C91", "C92", "C93", "C94", "C95",
                     "C96", "D45", "D47")
ICD9_lymphatic <- c("1691", "2000", "2001", "2008", "2014", "2015", "2016", "2017", "2019", "2020", "2021", "2022",
                    "2024", "2028", "2030", "2040", "2041", "2049", "2050", "2051", "2059", "2061", "2081", "2089",
                    "2384", "2386", "2387", "207 ")
### secondary
ICD10_secondary <- c("C77", "C78", "C79")
ICD9_secondary <- c("1960", "1961", "1962", "1963", "1965", "1968", "1969", "1975", "1976", "1983", "1984", "1985",
                    "1986", "1988")
### other
ICD10_other <- c("C76", "C80", "C97", "D09.7", "D09.9", "D18", "D36", "D46", "D48.7", "D48.9")
ICD9_other <- c("1950", "1952", "1954", "1955", "1958", "1990", "1991", "2349", "2388", "2395", "2397", "195 ")




# make vectors for the column names to search within
ICD10_repl <- c("p40006_i0", "p40006_i1", "p40006_i2", "p40006_i3", "p40006_i4", "p40006_i5", 
                "p40006_i6", "p40006_i7", "p40006_i8", "p40006_i9", "p40006_i10", "p40006_i11",
                "p40006_i12", "p40006_i13", "p40006_i14", "p40006_i15", "p40006_i16", "p40006_i17",
                "p40006_i18", "p40006_i19", "p40006_i20", "p40006_i21")
ICD9_repl <- c("p40013_i0", "p40013_i1", "p40013_i2", "p40013_i3", "p40013_i4", "p40013_i5", 
               "p40013_i6", "p40013_i7", "p40013_i8", "p40013_i10", "p40013_i11", "p40013_i12", 
               "p40013_i14")


# use these vectors to search through the df and replace with whichever string category is required!
### oral
Type_of_cancer_instances_categories[, ICD10_repl] <- lapply(Type_of_cancer_instances_categories[, ICD10_repl],
                                                            function(x) replace(x, grepl(paste(ICD10_oral, collapse = '|'), x), "oral"))
Type_of_cancer_instances_categories[, ICD9_repl] <- lapply(Type_of_cancer_instances_categories[, ICD9_repl],
                                                           function(x) replace(x, grepl(paste(ICD9_oral, collapse = '|'), x), "oral"))
### digestive
Type_of_cancer_instances_categories[, ICD10_repl] <- lapply(Type_of_cancer_instances_categories[, ICD10_repl],
                                                            function(x) replace(x, grepl(paste(ICD10_digestive, collapse = '|'), x), "digestive"))
Type_of_cancer_instances_categories[, ICD9_repl] <- lapply(Type_of_cancer_instances_categories[, ICD9_repl],
                                                           function(x) replace(x, grepl(paste(ICD9_digestive, collapse = '|'), x), "digestive"))
### respiratory
Type_of_cancer_instances_categories[, ICD10_repl] <- lapply(Type_of_cancer_instances_categories[, ICD10_repl],
                                                            function(x) replace(x, grepl(paste(ICD10_respiratory, collapse = '|'), x), "respiratory"))
Type_of_cancer_instances_categories[, ICD9_repl] <- lapply(Type_of_cancer_instances_categories[, ICD9_repl],
                                                           function(x) replace(x, grepl(paste(ICD9_respiratory, collapse = '|'), x), "respiratory"))
### bone
Type_of_cancer_instances_categories[, ICD10_repl] <- lapply(Type_of_cancer_instances_categories[, ICD10_repl],
                                                            function(x) replace(x, grepl(paste(ICD10_bone, collapse = '|'), x), "bone"))
Type_of_cancer_instances_categories[, ICD9_repl] <- lapply(Type_of_cancer_instances_categories[, ICD9_repl],
                                                           function(x) replace(x, grepl(paste(ICD9_bone, collapse = '|'), x), "bone"))
### skin
Type_of_cancer_instances_categories[, ICD10_repl] <- lapply(Type_of_cancer_instances_categories[, ICD10_repl],
                                                            function(x) replace(x, grepl(paste(ICD10_skin, collapse = '|'), x), "skin"))
Type_of_cancer_instances_categories[, ICD9_repl] <- lapply(Type_of_cancer_instances_categories[, ICD9_repl],
                                                           function(x) replace(x, grepl(paste(ICD9_skin, collapse = '|'), x), "skin"))
### mesothelium
Type_of_cancer_instances_categories[, ICD10_repl] <- lapply(Type_of_cancer_instances_categories[, ICD10_repl],
                                                            function(x) replace(x, grepl(paste(ICD10_mesothelium, collapse = '|'), x), "mesothelium"))
Type_of_cancer_instances_categories[, ICD9_repl] <- lapply(Type_of_cancer_instances_categories[, ICD9_repl],
                                                           function(x) replace(x, grepl(paste(ICD9_mesothelium, collapse = '|'), x), "mesothelium"))
### breast
Type_of_cancer_instances_categories[, ICD10_repl] <- lapply(Type_of_cancer_instances_categories[, ICD10_repl],
                                                            function(x) replace(x, grepl(paste(ICD10_breast, collapse = '|'), x), "breast"))
Type_of_cancer_instances_categories[, ICD9_repl] <- lapply(Type_of_cancer_instances_categories[, ICD9_repl],
                                                           function(x) replace(x, grepl(paste(ICD9_breast, collapse = '|'), x), "breast"))
### femalegenital
Type_of_cancer_instances_categories[, ICD10_repl] <- lapply(Type_of_cancer_instances_categories[, ICD10_repl],
                                                            function(x) replace(x, grepl(paste(ICD10_femalegenital, collapse = '|'), x), "femalegenital"))
Type_of_cancer_instances_categories[, ICD9_repl] <- lapply(Type_of_cancer_instances_categories[, ICD9_repl],
                                                           function(x) replace(x, grepl(paste(ICD9_femalegenital, collapse = '|'), x), "femalegenital"))
### malegenital
Type_of_cancer_instances_categories[, ICD10_repl] <- lapply(Type_of_cancer_instances_categories[, ICD10_repl],
                                                            function(x) replace(x, grepl(paste(ICD10_malegenital, collapse = '|'), x), "malegenital"))
Type_of_cancer_instances_categories[, ICD9_repl] <- lapply(Type_of_cancer_instances_categories[, ICD9_repl],
                                                           function(x) replace(x, grepl(paste(ICD9_malegenital, collapse = '|'), x), "malegenital"))
### urinary
Type_of_cancer_instances_categories[, ICD10_repl] <- lapply(Type_of_cancer_instances_categories[, ICD10_repl],
                                                            function(x) replace(x, grepl(paste(ICD10_urinary, collapse = '|'), x), "urinary"))
Type_of_cancer_instances_categories[, ICD9_repl] <- lapply(Type_of_cancer_instances_categories[, ICD9_repl],
                                                           function(x) replace(x, grepl(paste(ICD9_urinary, collapse = '|'), x), "urinary"))
### cns
Type_of_cancer_instances_categories[, ICD10_repl] <- lapply(Type_of_cancer_instances_categories[, ICD10_repl],
                                                            function(x) replace(x, grepl(paste(ICD10_cns, collapse = '|'), x), "cns"))
Type_of_cancer_instances_categories[, ICD9_repl] <- lapply(Type_of_cancer_instances_categories[, ICD9_repl],
                                                           function(x) replace(x, grepl(paste(ICD9_cns, collapse = '|'), x), "cns"))
### endocrine
Type_of_cancer_instances_categories[, ICD10_repl] <- lapply(Type_of_cancer_instances_categories[, ICD10_repl],
                                                            function(x) replace(x, grepl(paste(ICD10_endocrine, collapse = '|'), x), "endocrine"))
Type_of_cancer_instances_categories[, ICD9_repl] <- lapply(Type_of_cancer_instances_categories[, ICD9_repl],
                                                           function(x) replace(x, grepl(paste(ICD9_endocrine, collapse = '|'), x), "endocrine"))
### lymphatic
Type_of_cancer_instances_categories[, ICD10_repl] <- lapply(Type_of_cancer_instances_categories[, ICD10_repl],
                                                            function(x) replace(x, grepl(paste(ICD10_lymphatic, collapse = '|'), x), "lymphatic"))
Type_of_cancer_instances_categories[, ICD9_repl] <- lapply(Type_of_cancer_instances_categories[, ICD9_repl],
                                                           function(x) replace(x, grepl(paste(ICD9_lymphatic, collapse = '|'), x), "lymphatic"))
### secondary
Type_of_cancer_instances_categories[, ICD10_repl] <- lapply(Type_of_cancer_instances_categories[, ICD10_repl],
                                                            function(x) replace(x, grepl(paste(ICD10_secondary, collapse = '|'), x), "secondary"))
Type_of_cancer_instances_categories[, ICD9_repl] <- lapply(Type_of_cancer_instances_categories[, ICD9_repl],
                                                           function(x) replace(x, grepl(paste(ICD9_secondary, collapse = '|'), x), "secondary"))
### other
Type_of_cancer_instances_categories[, ICD10_repl] <- lapply(Type_of_cancer_instances_categories[, ICD10_repl],
                                                            function(x) replace(x, grepl(paste(ICD10_other, collapse = '|'), x), "other"))
Type_of_cancer_instances_categories[, ICD9_repl] <- lapply(Type_of_cancer_instances_categories[, ICD9_repl],
                                                           function(x) replace(x, grepl(paste(ICD9_other, collapse = '|'), x), "other"))


# remove molar pregnancy and other unrelated related terms

ICD10_molar <- c("O01.9", "O01.1", "O01.0")
# use codes to make these variables NA
Type_of_cancer_instances_categories[, ICD10_repl] <- lapply(Type_of_cancer_instances_categories[, ICD10_repl],
                                                            function(x) replace(x, grepl(paste(ICD10_molar, collapse = '|'), x), NA))

ICD9_molar <- c("6342", "630 ", "2532")
# use codes to make these variables NA
Type_of_cancer_instances_categories[, ICD9_repl] <- lapply(Type_of_cancer_instances_categories[, ICD9_repl],
                                                           function(x) replace(x, grepl(paste(ICD9_molar, collapse = '|'), x), NA))



# check that there are no guys I missed out!
unique(Type_of_cancer_instances_categories$p40006_i0)
unique(Type_of_cancer_instances_categories$p40006_i1)
unique(Type_of_cancer_instances_categories$p40006_i2)
unique(Type_of_cancer_instances_categories$p40006_i3)
unique(Type_of_cancer_instances_categories$p40006_i4)
unique(Type_of_cancer_instances_categories$p40006_i5)
unique(Type_of_cancer_instances_categories$p40006_i6)
unique(Type_of_cancer_instances_categories$p40006_i7)
unique(Type_of_cancer_instances_categories$p40006_i8)
unique(Type_of_cancer_instances_categories$p40006_i9)
unique(Type_of_cancer_instances_categories$p40006_i10)
unique(Type_of_cancer_instances_categories$p40006_i11)
unique(Type_of_cancer_instances_categories$p40006_i12)
unique(Type_of_cancer_instances_categories$p40006_i13)
unique(Type_of_cancer_instances_categories$p40006_i14)
unique(Type_of_cancer_instances_categories$p40006_i15)
unique(Type_of_cancer_instances_categories$p40006_i16)
unique(Type_of_cancer_instances_categories$p40006_i17)
unique(Type_of_cancer_instances_categories$p40006_i18)
unique(Type_of_cancer_instances_categories$p40006_i19)
unique(Type_of_cancer_instances_categories$p40006_i20)
unique(Type_of_cancer_instances_categories$p40006_i21)

unique(Type_of_cancer_instances_categories$p40013_i0)
unique(Type_of_cancer_instances_categories$p40013_i1)
unique(Type_of_cancer_instances_categories$p40013_i2)
unique(Type_of_cancer_instances_categories$p40013_i3)
unique(Type_of_cancer_instances_categories$p40013_i4)
unique(Type_of_cancer_instances_categories$p40013_i5)
unique(Type_of_cancer_instances_categories$p40013_i6)
unique(Type_of_cancer_instances_categories$p40013_i7)
unique(Type_of_cancer_instances_categories$p40013_i8)
unique(Type_of_cancer_instances_categories$p40013_i10)
unique(Type_of_cancer_instances_categories$p40013_i11)
unique(Type_of_cancer_instances_categories$p40013_i12)
unique(Type_of_cancer_instances_categories$p40013_i14)


# now check sum of NA again

Type_of_cancer_instances_categories$NA_sums <- rowSums(is.na(Type_of_cancer_instances_categories))
# have a look at summary stats
summary(Type_of_cancer_instances_categories$NA_sums)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 13.00   33.00   34.00   33.64   34.00   35.00 

# remove the rows with 35 NA values
Type_of_cancer_instances_categories <- subset(Type_of_cancer_instances_categories, NA_sums < 35)
# from 123772 to 123703 obs




# # now check that the codes match the categories
# bone_eid_i1 <- Type_of_cancer_instances_categories$eid[Type_of_cancer_instances_categories$p40006_i1 == "bone"]
# bone_eid_i1 <- bone_eid_i1[!is.na(bone_eid_i1)]
# bone_eid_i1
# 
# Type_of_cancer_instances$p40006_i1[Type_of_cancer_instances$eid == 4639044]
# # [1] 1567168 4286754 1895315 2085757 5685203 5761416 5091956 4159459 2645067 1892553 
# # 4129901 2385373 5890858 3985560 1616388 2675602 3711643 5414158 2147611
# # [20] 1805627 5500020 3895663 1208167 1286691 4365907 5429698 4639044
# # these all match with what I expect - won't do a more thorough check as everything seems fine



## now to figure out how to keep all cancer diagnosis data

# oral
Type_of_cancer_instances_categories$Neoplasm_oral <- 
  apply(Type_of_cancer_instances_categories[, c('p40006_i0', 'p40006_i1', 'p40006_i2', 'p40006_i3', 'p40006_i4', 'p40006_i5',
                                                'p40006_i6', 'p40006_i7', 'p40006_i8', 'p40006_i9', 'p40006_i10', 'p40006_i11',
                                                'p40006_i12', 'p40006_i13', 'p40006_i14', 'p40006_i15', 'p40006_i16', 'p40006_i17',
                                                'p40006_i18', 'p40006_i19', 'p40006_i20', 'p40006_i21', 'p40013_i0', 'p40013_i1', 'p40013_i2', 
                                                'p40013_i3', 'p40013_i4', 'p40013_i5', 'p40013_i6', 'p40013_i7', 'p40013_i8', 
                                                'p40013_i10', 'p40013_i11', 'p40013_i12', 'p40013_i14')], 1, function(row) any(row == "oral"))
# bone 
Type_of_cancer_instances_categories$Neoplasm_bone <- 
  apply(Type_of_cancer_instances_categories[, c('p40006_i0', 'p40006_i1', 'p40006_i2', 'p40006_i3', 'p40006_i4', 'p40006_i5',
                                                'p40006_i6', 'p40006_i7', 'p40006_i8', 'p40006_i9', 'p40006_i10', 'p40006_i11',
                                                'p40006_i12', 'p40006_i13', 'p40006_i14', 'p40006_i15', 'p40006_i16', 'p40006_i17',
                                                'p40006_i18', 'p40006_i19', 'p40006_i20', 'p40006_i21', 'p40013_i0', 'p40013_i1', 'p40013_i2', 
                                                'p40013_i3', 'p40013_i4', 'p40013_i5', 'p40013_i6', 'p40013_i7', 'p40013_i8', 
                                                'p40013_i10', 'p40013_i11', 'p40013_i12', 'p40013_i14')], 1, function(row) any(row == "bone"))
# breast   
Type_of_cancer_instances_categories$Neoplasm_breast <- 
  apply(Type_of_cancer_instances_categories[, c('p40006_i0', 'p40006_i1', 'p40006_i2', 'p40006_i3', 'p40006_i4', 'p40006_i5',
                                                'p40006_i6', 'p40006_i7', 'p40006_i8', 'p40006_i9', 'p40006_i10', 'p40006_i11',
                                                'p40006_i12', 'p40006_i13', 'p40006_i14', 'p40006_i15', 'p40006_i16', 'p40006_i17',
                                                'p40006_i18', 'p40006_i19', 'p40006_i20', 'p40006_i21', 'p40013_i0', 'p40013_i1', 'p40013_i2', 
                                                'p40013_i3', 'p40013_i4', 'p40013_i5', 'p40013_i6', 'p40013_i7', 'p40013_i8', 
                                                'p40013_i10', 'p40013_i11', 'p40013_i12', 'p40013_i14')], 1, function(row) any(row == "breast"))
# cns     
Type_of_cancer_instances_categories$Neoplasm_cns <- 
  apply(Type_of_cancer_instances_categories[, c('p40006_i0', 'p40006_i1', 'p40006_i2', 'p40006_i3', 'p40006_i4', 'p40006_i5',
                                                'p40006_i6', 'p40006_i7', 'p40006_i8', 'p40006_i9', 'p40006_i10', 'p40006_i11',
                                                'p40006_i12', 'p40006_i13', 'p40006_i14', 'p40006_i15', 'p40006_i16', 'p40006_i17',
                                                'p40006_i18', 'p40006_i19', 'p40006_i20', 'p40006_i21', 'p40013_i0', 'p40013_i1', 'p40013_i2', 
                                                'p40013_i3', 'p40013_i4', 'p40013_i5', 'p40013_i6', 'p40013_i7', 'p40013_i8', 
                                                'p40013_i10', 'p40013_i11', 'p40013_i12', 'p40013_i14')], 1, function(row) any(row == "cns"))
# digestive 
Type_of_cancer_instances_categories$Neoplasm_digestive <- 
  apply(Type_of_cancer_instances_categories[, c('p40006_i0', 'p40006_i1', 'p40006_i2', 'p40006_i3', 'p40006_i4', 'p40006_i5',
                                                'p40006_i6', 'p40006_i7', 'p40006_i8', 'p40006_i9', 'p40006_i10', 'p40006_i11',
                                                'p40006_i12', 'p40006_i13', 'p40006_i14', 'p40006_i15', 'p40006_i16', 'p40006_i17',
                                                'p40006_i18', 'p40006_i19', 'p40006_i20', 'p40006_i21', 'p40013_i0', 'p40013_i1', 'p40013_i2', 
                                                'p40013_i3', 'p40013_i4', 'p40013_i5', 'p40013_i6', 'p40013_i7', 'p40013_i8', 
                                                'p40013_i10', 'p40013_i11', 'p40013_i12', 'p40013_i14')], 1, function(row) any(row == "digestive"))
# endocrine   
Type_of_cancer_instances_categories$Neoplasm_endocrine <- 
  apply(Type_of_cancer_instances_categories[, c('p40006_i0', 'p40006_i1', 'p40006_i2', 'p40006_i3', 'p40006_i4', 'p40006_i5',
                                                'p40006_i6', 'p40006_i7', 'p40006_i8', 'p40006_i9', 'p40006_i10', 'p40006_i11',
                                                'p40006_i12', 'p40006_i13', 'p40006_i14', 'p40006_i15', 'p40006_i16', 'p40006_i17',
                                                'p40006_i18', 'p40006_i19', 'p40006_i20', 'p40006_i21', 'p40013_i0', 'p40013_i1', 'p40013_i2', 
                                                'p40013_i3', 'p40013_i4', 'p40013_i5', 'p40013_i6', 'p40013_i7', 'p40013_i8', 
                                                'p40013_i10', 'p40013_i11', 'p40013_i12', 'p40013_i14')], 1, function(row) any(row == "endocrine"))
# femalegenital 
Type_of_cancer_instances_categories$Neoplasm_femalegenital <- 
  apply(Type_of_cancer_instances_categories[, c('p40006_i0', 'p40006_i1', 'p40006_i2', 'p40006_i3', 'p40006_i4', 'p40006_i5',
                                                'p40006_i6', 'p40006_i7', 'p40006_i8', 'p40006_i9', 'p40006_i10', 'p40006_i11',
                                                'p40006_i12', 'p40006_i13', 'p40006_i14', 'p40006_i15', 'p40006_i16', 'p40006_i17',
                                                'p40006_i18', 'p40006_i19', 'p40006_i20', 'p40006_i21', 'p40013_i0', 'p40013_i1', 'p40013_i2', 
                                                'p40013_i3', 'p40013_i4', 'p40013_i5', 'p40013_i6', 'p40013_i7', 'p40013_i8', 
                                                'p40013_i10', 'p40013_i11', 'p40013_i12', 'p40013_i14')], 1, function(row) any(row == "femalegenital"))
# lymphatic    
Type_of_cancer_instances_categories$Neoplasm_lymphatic <- 
  apply(Type_of_cancer_instances_categories[, c('p40006_i0', 'p40006_i1', 'p40006_i2', 'p40006_i3', 'p40006_i4', 'p40006_i5',
                                                'p40006_i6', 'p40006_i7', 'p40006_i8', 'p40006_i9', 'p40006_i10', 'p40006_i11',
                                                'p40006_i12', 'p40006_i13', 'p40006_i14', 'p40006_i15', 'p40006_i16', 'p40006_i17',
                                                'p40006_i18', 'p40006_i19', 'p40006_i20', 'p40006_i21', 'p40013_i0', 'p40013_i1', 'p40013_i2', 
                                                'p40013_i3', 'p40013_i4', 'p40013_i5', 'p40013_i6', 'p40013_i7', 'p40013_i8', 
                                                'p40013_i10', 'p40013_i11', 'p40013_i12', 'p40013_i14')], 1, function(row) any(row == "lymphatic"))
# malegenital  
Type_of_cancer_instances_categories$Neoplasm_malegenital <- 
  apply(Type_of_cancer_instances_categories[, c('p40006_i0', 'p40006_i1', 'p40006_i2', 'p40006_i3', 'p40006_i4', 'p40006_i5',
                                                'p40006_i6', 'p40006_i7', 'p40006_i8', 'p40006_i9', 'p40006_i10', 'p40006_i11',
                                                'p40006_i12', 'p40006_i13', 'p40006_i14', 'p40006_i15', 'p40006_i16', 'p40006_i17',
                                                'p40006_i18', 'p40006_i19', 'p40006_i20', 'p40006_i21', 'p40013_i0', 'p40013_i1', 'p40013_i2', 
                                                'p40013_i3', 'p40013_i4', 'p40013_i5', 'p40013_i6', 'p40013_i7', 'p40013_i8', 
                                                'p40013_i10', 'p40013_i11', 'p40013_i12', 'p40013_i14')], 1, function(row) any(row == "malegenital"))
# mesothelium  
Type_of_cancer_instances_categories$Neoplasm_mesothelium <- 
  apply(Type_of_cancer_instances_categories[, c('p40006_i0', 'p40006_i1', 'p40006_i2', 'p40006_i3', 'p40006_i4', 'p40006_i5',
                                                'p40006_i6', 'p40006_i7', 'p40006_i8', 'p40006_i9', 'p40006_i10', 'p40006_i11',
                                                'p40006_i12', 'p40006_i13', 'p40006_i14', 'p40006_i15', 'p40006_i16', 'p40006_i17',
                                                'p40006_i18', 'p40006_i19', 'p40006_i20', 'p40006_i21', 'p40013_i0', 'p40013_i1', 'p40013_i2', 
                                                'p40013_i3', 'p40013_i4', 'p40013_i5', 'p40013_i6', 'p40013_i7', 'p40013_i8', 
                                                'p40013_i10', 'p40013_i11', 'p40013_i12', 'p40013_i14')], 1, function(row) any(row == "mesothelium"))
# other  
Type_of_cancer_instances_categories$Neoplasm_other <- 
  apply(Type_of_cancer_instances_categories[, c('p40006_i0', 'p40006_i1', 'p40006_i2', 'p40006_i3', 'p40006_i4', 'p40006_i5',
                                                'p40006_i6', 'p40006_i7', 'p40006_i8', 'p40006_i9', 'p40006_i10', 'p40006_i11',
                                                'p40006_i12', 'p40006_i13', 'p40006_i14', 'p40006_i15', 'p40006_i16', 'p40006_i17',
                                                'p40006_i18', 'p40006_i19', 'p40006_i20', 'p40006_i21', 'p40013_i0', 'p40013_i1', 'p40013_i2', 
                                                'p40013_i3', 'p40013_i4', 'p40013_i5', 'p40013_i6', 'p40013_i7', 'p40013_i8', 
                                                'p40013_i10', 'p40013_i11', 'p40013_i12', 'p40013_i14')], 1, function(row) any(row == "other"))
# respiratory 
Type_of_cancer_instances_categories$Neoplasm_respiratory <- 
  apply(Type_of_cancer_instances_categories[, c('p40006_i0', 'p40006_i1', 'p40006_i2', 'p40006_i3', 'p40006_i4', 'p40006_i5',
                                                'p40006_i6', 'p40006_i7', 'p40006_i8', 'p40006_i9', 'p40006_i10', 'p40006_i11',
                                                'p40006_i12', 'p40006_i13', 'p40006_i14', 'p40006_i15', 'p40006_i16', 'p40006_i17',
                                                'p40006_i18', 'p40006_i19', 'p40006_i20', 'p40006_i21', 'p40013_i0', 'p40013_i1', 'p40013_i2', 
                                                'p40013_i3', 'p40013_i4', 'p40013_i5', 'p40013_i6', 'p40013_i7', 'p40013_i8', 
                                                'p40013_i10', 'p40013_i11', 'p40013_i12', 'p40013_i14')], 1, function(row) any(row == "respiratory"))

# secondary   
Type_of_cancer_instances_categories$Neoplasm_secondary <- 
  apply(Type_of_cancer_instances_categories[, c('p40006_i0', 'p40006_i1', 'p40006_i2', 'p40006_i3', 'p40006_i4', 'p40006_i5',
                                                'p40006_i6', 'p40006_i7', 'p40006_i8', 'p40006_i9', 'p40006_i10', 'p40006_i11',
                                                'p40006_i12', 'p40006_i13', 'p40006_i14', 'p40006_i15', 'p40006_i16', 'p40006_i17',
                                                'p40006_i18', 'p40006_i19', 'p40006_i20', 'p40006_i21', 'p40013_i0', 'p40013_i1', 'p40013_i2', 
                                                'p40013_i3', 'p40013_i4', 'p40013_i5', 'p40013_i6', 'p40013_i7', 'p40013_i8', 
                                                'p40013_i10', 'p40013_i11', 'p40013_i12', 'p40013_i14')], 1, function(row) any(row == "secondary"))

# skin    
Type_of_cancer_instances_categories$Neoplasm_skin <- 
  apply(Type_of_cancer_instances_categories[, c('p40006_i0', 'p40006_i1', 'p40006_i2', 'p40006_i3', 'p40006_i4', 'p40006_i5',
                                                'p40006_i6', 'p40006_i7', 'p40006_i8', 'p40006_i9', 'p40006_i10', 'p40006_i11',
                                                'p40006_i12', 'p40006_i13', 'p40006_i14', 'p40006_i15', 'p40006_i16', 'p40006_i17',
                                                'p40006_i18', 'p40006_i19', 'p40006_i20', 'p40006_i21', 'p40013_i0', 'p40013_i1', 'p40013_i2', 
                                                'p40013_i3', 'p40013_i4', 'p40013_i5', 'p40013_i6', 'p40013_i7', 'p40013_i8', 
                                                'p40013_i10', 'p40013_i11', 'p40013_i12', 'p40013_i14')], 1, function(row) any(row == "skin"))

# urinary        
Type_of_cancer_instances_categories$Neoplasm_urinary <- 
  apply(Type_of_cancer_instances_categories[, c('p40006_i0', 'p40006_i1', 'p40006_i2', 'p40006_i3', 'p40006_i4', 'p40006_i5',
                                                'p40006_i6', 'p40006_i7', 'p40006_i8', 'p40006_i9', 'p40006_i10', 'p40006_i11',
                                                'p40006_i12', 'p40006_i13', 'p40006_i14', 'p40006_i15', 'p40006_i16', 'p40006_i17',
                                                'p40006_i18', 'p40006_i19', 'p40006_i20', 'p40006_i21', 'p40013_i0', 'p40013_i1', 'p40013_i2', 
                                                'p40013_i3', 'p40013_i4', 'p40013_i5', 'p40013_i6', 'p40013_i7', 'p40013_i8', 
                                                'p40013_i10', 'p40013_i11', 'p40013_i12', 'p40013_i14')], 1, function(row) any(row == "urinary"))



Type_of_cancer_updated <- Type_of_cancer_instances_categories[,c(1, 38:52)]

summary(Type_of_cancer_updated)
# eid          Neoplasm_urinary Neoplasm_oral  Neoplasm_bone  Neoplasm_breast Neoplasm_cns   Neoplasm_digestive Neoplasm_endocrine
# Min.   :1000020   Mode:logical     Mode:logical   Mode:logical   Mode:logical    Mode:logical   Mode:logical       Mode:logical      
# 1st Qu.:2259784   TRUE:6889        TRUE:1922      TRUE:169       TRUE:22456      TRUE:3142      TRUE:16048         TRUE:1551         
# Median :3509763   NA's:116814      NA's:121781    NA's:123534    NA's:101247     NA's:120561    NA's:107655        NA's:122152       
#  Mean   :3512347                                                                                                                      
#  3rd Qu.:4769455                                                                                                                      
#  Max.   :6024097                                                                                                                      
#  Neoplasm_femalegenital Neoplasm_lymphatic Neoplasm_malegenital Neoplasm_mesothelium Neoplasm_other Neoplasm_respiratory Neoplasm_secondary
#  Mode:logical           Mode:logical       Mode:logical         Mode:logical         Mode:logical   Mode:logical         Mode:logical      
#  TRUE:12372             TRUE:8594          TRUE:17055           TRUE:1426            TRUE:939       TRUE:6003            TRUE:589          
#  NA's:111331            NA's:115109        NA's:106648          NA's:122277          NA's:122764    NA's:117700          NA's:123114       
# 
# 
# 
# Neoplasm_skin 
# Mode:logical  
# TRUE:45067    
# NA's:78636

Type_of_cancer_updated$Neoplasm_urinary <- ifelse(Type_of_cancer_updated$Neoplasm_urinary =="TRUE",1,0)
Type_of_cancer_updated$Neoplasm_oral <- ifelse(Type_of_cancer_updated$Neoplasm_oral =="TRUE",1,0)
Type_of_cancer_updated$Neoplasm_bone <- ifelse(Type_of_cancer_updated$Neoplasm_bone =="TRUE",1,0)
Type_of_cancer_updated$Neoplasm_breast <- ifelse(Type_of_cancer_updated$Neoplasm_breast =="TRUE",1,0)
Type_of_cancer_updated$Neoplasm_cns <- ifelse(Type_of_cancer_updated$Neoplasm_cns =="TRUE",1,0)
Type_of_cancer_updated$Neoplasm_digestive <- ifelse(Type_of_cancer_updated$Neoplasm_digestive =="TRUE",1,0)
Type_of_cancer_updated$Neoplasm_endocrine <- ifelse(Type_of_cancer_updated$Neoplasm_endocrine =="TRUE",1,0)
Type_of_cancer_updated$Neoplasm_femalegenital <- ifelse(Type_of_cancer_updated$Neoplasm_femalegenital =="TRUE",1,0)
Type_of_cancer_updated$Neoplasm_lymphatic <- ifelse(Type_of_cancer_updated$Neoplasm_lymphatic =="TRUE",1,0)
Type_of_cancer_updated$Neoplasm_malegenital <- ifelse(Type_of_cancer_updated$Neoplasm_malegenital =="TRUE",1,0)
Type_of_cancer_updated$Neoplasm_mesothelium <- ifelse(Type_of_cancer_updated$Neoplasm_mesothelium =="TRUE",1,0)
Type_of_cancer_updated$Neoplasm_other <- ifelse(Type_of_cancer_updated$Neoplasm_other =="TRUE",1,0)
Type_of_cancer_updated$Neoplasm_respiratory <- ifelse(Type_of_cancer_updated$Neoplasm_respiratory =="TRUE",1,0)
Type_of_cancer_updated$Neoplasm_secondary <- ifelse(Type_of_cancer_updated$Neoplasm_secondary =="TRUE",1,0)
Type_of_cancer_updated$Neoplasm_skin <- ifelse(Type_of_cancer_updated$Neoplasm_skin =="TRUE",1,0)

summary(Type_of_cancer_updated)
# eid          Neoplasm_urinary Neoplasm_oral    Neoplasm_bone    Neoplasm_breast   Neoplasm_cns    Neoplasm_digestive
# Min.   :1000020   Min.   :1        Min.   :1        Min.   :1        Min.   :1        Min.   :1        Min.   :1         
# 1st Qu.:2259784   1st Qu.:1        1st Qu.:1        1st Qu.:1        1st Qu.:1        1st Qu.:1        1st Qu.:1         
# Median :3509763   Median :1        Median :1        Median :1        Median :1        Median :1        Median :1         
# Mean   :3512347   Mean   :1        Mean   :1        Mean   :1        Mean   :1        Mean   :1        Mean   :1         
# 3rd Qu.:4769455   3rd Qu.:1        3rd Qu.:1        3rd Qu.:1        3rd Qu.:1        3rd Qu.:1        3rd Qu.:1         
# Max.   :6024097   Max.   :1        Max.   :1        Max.   :1        Max.   :1        Max.   :1        Max.   :1         
# NA's   :116814   NA's   :121781   NA's   :123534   NA's   :101247   NA's   :120561   NA's   :107655    
# Neoplasm_endocrine Neoplasm_femalegenital Neoplasm_lymphatic Neoplasm_malegenital Neoplasm_mesothelium Neoplasm_other  
# Min.   :1          Min.   :1              Min.   :1          Min.   :1            Min.   :1            Min.   :1       
# 1st Qu.:1          1st Qu.:1              1st Qu.:1          1st Qu.:1            1st Qu.:1            1st Qu.:1       
# Median :1          Median :1              Median :1          Median :1            Median :1            Median :1       
# Mean   :1          Mean   :1              Mean   :1          Mean   :1            Mean   :1            Mean   :1       
# 3rd Qu.:1          3rd Qu.:1              3rd Qu.:1          3rd Qu.:1            3rd Qu.:1            3rd Qu.:1       
# Max.   :1          Max.   :1              Max.   :1          Max.   :1            Max.   :1            Max.   :1       
# NA's   :122152     NA's   :111331         NA's   :115109     NA's   :106648       NA's   :122277       NA's   :122764  
# Neoplasm_respiratory Neoplasm_secondary Neoplasm_skin  
# Min.   :1            Min.   :1          Min.   :1      
# 1st Qu.:1            1st Qu.:1          1st Qu.:1      
# Median :1            Median :1          Median :1      
# Mean   :1            Mean   :1          Mean   :1      
# 3rd Qu.:1            3rd Qu.:1          3rd Qu.:1      
# Max.   :1            Max.   :1          Max.   :1      
# NA's   :117700       NA's   :123114     NA's   :78636

Type_of_cancer_updated[is.na(Type_of_cancer_updated)] <- 0
# eid          Neoplasm_urinary  Neoplasm_oral     Neoplasm_bone      Neoplasm_breast   Neoplasm_cns    Neoplasm_digestive
# Min.   :1000020   Min.   :0.00000   Min.   :0.00000   Min.   :0.000000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000    
# 1st Qu.:2259784   1st Qu.:0.00000   1st Qu.:0.00000   1st Qu.:0.000000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000    
# Median :3509763   Median :0.00000   Median :0.00000   Median :0.000000   Median :0.0000   Median :0.0000   Median :0.0000    
# Mean   :3512347   Mean   :0.05569   Mean   :0.01554   Mean   :0.001366   Mean   :0.1815   Mean   :0.0254   Mean   :0.1297    
# 3rd Qu.:4769455   3rd Qu.:0.00000   3rd Qu.:0.00000   3rd Qu.:0.000000   3rd Qu.:0.0000   3rd Qu.:0.0000   3rd Qu.:0.0000    
# Max.   :6024097   Max.   :1.00000   Max.   :1.00000   Max.   :1.000000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000    
# Neoplasm_endocrine Neoplasm_femalegenital Neoplasm_lymphatic Neoplasm_malegenital Neoplasm_mesothelium Neoplasm_other    
# Min.   :0.00000    Min.   :0.0            Min.   :0.00000    Min.   :0.0000       Min.   :0.00000      Min.   :0.000000  
# 1st Qu.:0.00000    1st Qu.:0.0            1st Qu.:0.00000    1st Qu.:0.0000       1st Qu.:0.00000      1st Qu.:0.000000  
# Median :0.00000    Median :0.0            Median :0.00000    Median :0.0000       Median :0.00000      Median :0.000000  
# Mean   :0.01254    Mean   :0.1            Mean   :0.06947    Mean   :0.1379       Mean   :0.01153      Mean   :0.007591  
# 3rd Qu.:0.00000    3rd Qu.:0.0            3rd Qu.:0.00000    3rd Qu.:0.0000       3rd Qu.:0.00000      3rd Qu.:0.000000  
# Max.   :1.00000    Max.   :1.0            Max.   :1.00000    Max.   :1.0000       Max.   :1.00000      Max.   :1.000000  
# Neoplasm_respiratory Neoplasm_secondary Neoplasm_skin   
# Min.   :0.00000      Min.   :0.000000   Min.   :0.0000  
# 1st Qu.:0.00000      1st Qu.:0.000000   1st Qu.:0.0000  
# Median :0.00000      Median :0.000000   Median :0.0000  
# Mean   :0.04853      Mean   :0.004761   Mean   :0.3643  
# 3rd Qu.:0.00000      3rd Qu.:0.000000   3rd Qu.:1.0000  
# Max.   :1.00000      Max.   :1.000000   Max.   :1.0000



# add these columns to the EntireCohort dataset:
system('dx download Emily-folder/EntireCohort_ready.csv')
EntireCohort_ready <- read_csv("EntireCohort_ready.csv")
EntireCohort_ready <- EntireCohort_ready[,-1]

EntireCohort_ready_2 <- left_join(EntireCohort_ready, Type_of_cancer_updated, by = 'eid')

#save it to R
write.csv(EntireCohort_ready_2, "EntireCohort_ready.csv")
# now to save this to UKB project
system('dx upload EntireCohort_ready.csv --path Emily-folder/EntireCohort_ready.csv')



# save script! -----------------------------------------------------------------


system('dx upload EntireCohort_script.R --path Emily-folder/EntireCohort_script.R')

