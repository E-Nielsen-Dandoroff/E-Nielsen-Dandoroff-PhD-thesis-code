# 12 01 2024
# updated ID extraction to retain variant information



# using Emily-folder/MGSvar/MGSvar_script2_END.R as a reference


# Load packages
library(data.table)
library(dplyr)

# DONSON IDs --------------------------------------------------------------


# Load PLINK format SNP names for MGSvar SNPs
donsonMGSsnps = readLines('/mnt/project/plink-snp-data-20230405/donson-MGSvar-snp-ids.txt')
head(donsonMGSsnps)
# as expected: "chr:pos:ref:alt" "21:33584744:G:A" "21:33586090:A:G" "21:33584804:C:T" "21:33583643:T:C" "21:33578374:G:A"

# Get rid of first entry ("chr:pos:ref:alt")
donsonMGSsnps = donsonMGSsnps[-1]
head(donsonMGSsnps)

# Use data.table (faster) to load the raw data from PLINK
donsonMGSraw = fread('/mnt/project/plink-snp-data-20230405/donson-MGSvar-snp-data/donson-MGSvar-snps.raw', data.table=FALSE)

#checking
class(donsonMGSraw)
# "data.frame"
dim(donsonMGSraw)
# 469835     13

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
donsonMGSraw[1:5,1:8]

# We only want the MGS SNPs
# First need to get the SNP id formats matching
head(donsonMGSsnps)
# "21:33584744:G:A" "21:33586090:A:G" "21:33584804:C:T" "21:33583643:T:C" "21:33578374:G:A" "21:33584705:G:A"
colnames(donsonMGSraw)[7:10]
# "21:33578374:G:A_A" "21:33578374:G:T_T" "21:33581355:G:A_A" "21:33581355:G:T_T"

# Create new variable with extra "_X" added:
donsonMGSsnp_ids = strsplit(donsonMGSsnps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(donsonMGSsnp_ids)
#"21:33584744:G:A_A" "21:33586090:A:G_G" "21:33584804:C:T_T" "21:33583643:T:C_C" "21:33578374:G:A_A" "21:33584705:G:A_A"



## checking frequency data 


# removing the first 6 columns for the next part of the analysis
donsongenotypes <- donsonMGSraw[, -c(1:6)]

# How many 0's, 1's and 2's are there across alls SNPs and samples?
# minor allele homoygotes
sum(donsongenotypes == 2, na.rm=TRUE)
# 0 (perfect)

# heterozygotes
sum(donsongenotypes == 1, na.rm=TRUE)
# 56 (looks good!)

# major allele homoygotes
sum(donsongenotypes == 0, na.rm=TRUE)
# 3288780

# Get counts of 1's and 2's (i.e., any minor allele) per SNP
donson_allele_counts = colSums(donsongenotypes!=0, na.rm=TRUE)
donson_allele_counts
# 21:33578374:G:A_A 21:33578374:G:T_T 21:33581355:G:A_A 21:33581355:G:T_T 21:33584705:G:A_A 21:33584744:G:A_A 21:33586090:A:G_G 
# 8                 1                 2                 2                20                19                 4 
# looks perfect! perfectly matches the preliminary look!

# i don't think I need to check against plink - set of variants is small enough that i can see that everything is good :)









## removing variants that I don't care about
# because plink extracts all recorded variants at each location i'm interested in!

# want to duplicate this dataframe to mess around with it
donsonMGSvar <- donsonMGSraw

# Which columns of the raw data match SNPs of interest?
match_ids_donson = match(donsonMGSsnp_ids, colnames(donsonMGSraw))
match_ids_donson
# 12 13 NA NA  7 11  9

#positions of extra variants 
all_names_donson <- colnames(donsonMGSraw)
all_names_donson

setdiff(all_names_donson, donsonMGSsnp_ids)
# "FID"               "IID"               "PAT"               "MAT"               "SEX"               "PHENOTYPE"        
# "21:33578374:G:T_T" "21:33581355:G:T_T"

# find location of the extra variant columns to delete only these
extra_variants_donson <- setdiff(all_names_donson, donsonMGSsnp_ids)
match(extra_variants_donson, colnames(donsonMGSvar))
# 1  2  3  4  5  6  8 10
# 1:6 are the other data columns
donsonMGSvar <- donsonMGSvar[, -c(8, 10)]

match(extra_variants_donson, colnames(donsonMGSvar))
# 1  2  3  4  5  6 NA NA

# have removed the unwanted variants from donsonMGSraw_2


###
# extracting list of patient IDs
# add a new column that contains the sum of SNP columns for each row(participant)
donsonMGSvar$Sums <- rowSums(donsonMGSvar[ , 7:11], na.rm = TRUE)
# remove the rows that have $Sum = 0
donsonMGSvar <- donsonMGSvar[donsonMGSvar$Sums !=0, ]
# 654 obs. (perfect!)
unique(donsonMGSvar[,12])
# 1 (each individual only has the 1 variant)

# make a list of participant IDs w relevant SNPs
DONSON_MGSvar_IDs <- donsonMGSvar[, 1]

### figured out that its essential to change to character - or will save as chunks of entries
DONSON_MGSvar_IDs <- as.character(DONSON_MGSvar_IDs)
###

#save 'ORC1_MGSvar_IDs'
write(DONSON_MGSvar_IDs, "DONSON_MGSvar_IDs.txt")

# now to save this to UKB project
system('dx upload DONSON_MGSvar_IDs.txt --path Emily-folder/MGSvar/DONSON_MGSvar_IDs.txt')




### making MGSvar_IDs_VARs.csv
# working with donsonMGSvar
# no individuals have more than one of the varaiants
# but this may not be the case for other times I want to use this code

# make a copy to mess around with
donsonMGSvar_test1 <- donsonMGSvar

donsonMGSvar_test1$VAR <- names(donsonMGSvar_test1[,c(7:11)][max.col(donsonMGSvar_test1[,c(7:11)])])
# it works!!!
# but it adds a number to the end of some of the variants, not ideal
# this code ammends this:
donsonMGSvar_test1$VAR <- sub('\\..*', '', donsonMGSvar_test1$VAR)

# condense it down to just the ID number and the variant information
donsonMGSvar_test2 <- donsonMGSvar_test1[,c(1, 12, 13)]
# can save this now to annotate the entire cohort dataframe

#save it to R
write.csv(donsonMGSvar_test2, "DONSON_MGSvar_IDs_VARs.csv")
# now to save this to UKB project
system('dx upload DONSON_MGSvar_IDs_VARs.csv --path Emily-folder/MGSvar/DONSON_MGSvar_IDs_VARs.csv')















# GINS2 and GINS3 IDs -----------------------------------------------------


# Load PLINK format SNP names for MGSvar SNPs
ginsMGSsnps = readLines('/mnt/project/plink-snp-data-20230405/gins-MGSvar-snp-ids.txt')
# 16 total
head(ginsMGSsnps)
# need to remove "chr:pos:ref:alt"
ginsMGSsnps <- ginsMGSsnps[-1]

# Use data.table (faster) to load the raw data from PLINK
ginsMGSraw = fread('/mnt/project/plink-snp-data-20230405/gins-MGSvar-snp-data/gins-MGSvar-snps.raw', data.table=FALSE)

class(ginsMGSraw)
# "data.frame"
dim(ginsMGSraw)
# 469835      9


# We only want the MGS SNPs
# First need to get the SNP id formats matching
# Create new variable with extra "_X" added:
ginsMGSsnp_ids = strsplit(ginsMGSsnps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(ginsMGSsnp_ids)

# Which columns of the raw data match SNPs of interest?
match_ids_gins = match(ginsMGSsnp_ids, colnames(ginsMGSraw))
match_ids_gins
# 8 NA  7 NA


#positions of extra variants 
all_names_gins <- colnames(ginsMGSraw)
all_names_gins
# "FID"               "IID"               "PAT"               "MAT"               "SEX"               "PHENOTYPE"        
# "16:58403156:G:A_A" "16:85678631:C:A_A" "16:85678631:C:T_T"

setdiff(all_names_gins, ginsMGSsnp_ids)
# "FID"               "IID"               "PAT"               "MAT"               "SEX"               "PHENOTYPE"        
# "16:85678631:C:T_T"

extra_variants_gins <- setdiff(all_names_gins, ginsMGSsnp_ids)
match(extra_variants_gins, colnames(ginsMGSraw))
# 1 2 3 4 5 6 9
# 1:6 are meant to be there

ginsMGSvar <- ginsMGSraw[, -c(9)]
match(extra_variants_gins, colnames(ginsMGSvar))
#great!

# removing the first 6 columns for the next part of the analysis
ginsgenotypes <- ginsMGSvar[, -c(1:6)]

# How many 0's, 1's and 2's are there across alls SNPs and samples?
# minor allele homoygotes
sum(ginsgenotypes == 2, na.rm=TRUE)
# 0 (perfect)

# heterozygotes
sum(ginsgenotypes == 1, na.rm=TRUE)
# 42 (looks good!)

# major allele homoygotes
sum(ginsgenotypes == 0, na.rm=TRUE)
# 39606

# Get counts of 1's and 2's (i.e., any minor allele) per SNP
gins_allele_counts = colSums(ginsgenotypes!=0, na.rm=TRUE)
gins_allele_counts
# 16:58403156:G:A_A 16:85678631:C:A_A 
# 21                21
# great


#
# separating by gene
install.packages("readr")
library(readr)
system('dx download Emily-folder/MGSvar/gene_MGSvar_csv/GINS2_GINS3_MGSvar.csv')
GINS2_GINS3_MGSvar <- read.csv("GINS2_GINS3_MGSvar.csv")
View(GINS2_GINS3_MGSvar)

# first making the variant IDs match
GINS2_GINS3_MGSvar$variants <- paste(GINS2_GINS3_MGSvar$chr, 
                                     GINS2_GINS3_MGSvar$pos,
                                     GINS2_GINS3_MGSvar$ref,
                                     GINS2_GINS3_MGSvar$alt,
                                     sep = ":")

GINS2_GINS3_MGSvar$variants <- paste(GINS2_GINS3_MGSvar$variants,
                                     GINS2_GINS3_MGSvar$alt,
                                     sep = "_")
GINS2_GINS3_MGSvar$variants
# "16:85678631:C:A_A" "16:58392672:A:G_G" "16:58403156:G:A_A" "16:58392671:G:A_A"
# looks good!

GINS2_variants <- GINS2_GINS3_MGSvar$variants[GINS2_GINS3_MGSvar$gene == "GINS2"]
GINS2_variants
# "16:85678631:C:A_A"

GINS3_variants <- GINS2_GINS3_MGSvar$variants[GINS2_GINS3_MGSvar$gene == "GINS3"]
GINS3_variants
# "16:58392672:A:G_G" "16:58403156:G:A_A" "16:58392671:G:A_A"

#df <- df[!(names(df) %in% rem)]
GINS2_MGSvar <- ginsMGSvar[!(names(ginsMGSvar) %in% GINS3_variants)]

GINS3_MGSvar <- ginsMGSvar[!(names(ginsMGSvar) %in% GINS2_variants)]

# to check:
colnames(GINS2_MGSvar)
# "FID"               "IID"               "PAT"               "MAT"               "SEX"               "PHENOTYPE"        
# "16:85678631:C:A_A"
GINS2_variants
# "16:85678631:C:A_A"

colnames(GINS3_MGSvar)
# "FID"               "IID"               "PAT"               "MAT"               "SEX"               "PHENOTYPE"        
# "16:58403156:G:A_A"
GINS3_variants
# "16:58392672:A:G_G" "16:58403156:G:A_A" "16:58392671:G:A_A"





###
# extracting list of patient IDs
# add a new column that contains the sum of SNP columns for each row(participant)

# the method I had previously used relied on summing multiple columns in a row
# but there is only one variant in each of these dfs
# so can just duplicate the column and rename
GINS2_MGSvar$Sums <- GINS2_MGSvar$`16:85678631:C:A_A`
GINS3_MGSvar$Sums <- GINS3_MGSvar$`16:58403156:G:A_A`


# remove the rows that have $Sum = 0
GINS2_MGSvar <- GINS2_MGSvar[GINS2_MGSvar$Sums !=0, ]
# 42 obs
# also need to remove NA values
GINS2_MGSvar <- GINS2_MGSvar[!is.na(GINS2_MGSvar$Sums),]
# 21 obs

GINS3_MGSvar <- GINS3_MGSvar[GINS3_MGSvar$Sums !=0, ]
# 22 obs
# also need to remove NA values
GINS3_MGSvar <- GINS3_MGSvar[!is.na(GINS3_MGSvar$Sums),]
# 21 obs

unique(GINS2_MGSvar[,8])
# 1 
unique(GINS3_MGSvar[,8])
# 1 

# make a list of participant IDs w relevant SNPs
GINS2_MGSvar_IDs <- GINS2_MGSvar[, 1]
GINS3_MGSvar_IDs <- GINS3_MGSvar[, 1]

GINS2_MGSvar_IDs <- as.character(GINS2_MGSvar_IDs)
GINS3_MGSvar_IDs <- as.character(GINS3_MGSvar_IDs)

# save
write(GINS2_MGSvar_IDs, "GINS2_MGSvar_IDs.txt")
write(GINS3_MGSvar_IDs, "GINS3_MGSvar_IDs.txt")

# now to save this to UKB project
system('dx upload GINS2_MGSvar_IDs.txt --path Emily-folder/MGSvar/GINS2_MGSvar_IDs.txt')
system('dx upload GINS3_MGSvar_IDs.txt --path Emily-folder/MGSvar/GINS3_MGSvar_IDs.txt')


# both only have one variant, so instead of using the method I'd used for DONSON,
# I can instead just paste the variant name into a VAR column!
# GINS2_MGSvar
GINS2_MGSvar$VAR <- "16:85678631:C:A_A"
# GINS3_MGSvar
GINS3_MGSvar$VAR <- "16:58403156:G:A_A"


# condense it down to just the ID number and the variant information
GINS2_MGSvar_test2 <- GINS2_MGSvar[,c(1, 8, 9)]
GINS3_MGSvar_test2 <- GINS3_MGSvar[,c(1, 8, 9)]
# can save this now to annotate the entire cohort dataframe

#save it to R
write.csv(GINS2_MGSvar_test2, "GINS2_MGSvar_IDs_VARs.csv")
write.csv(GINS3_MGSvar_test2, "GINS3_MGSvar_IDs_VARs.csv")

# now to save this to UKB project
system('dx upload GINS2_MGSvar_IDs_VARs.csv --path Emily-folder/MGSvar/GINS2_MGSvar_IDs_VARs.csv')
system('dx upload GINS3_MGSvar_IDs_VARs.csv --path Emily-folder/MGSvar/GINS3_MGSvar_IDs_VARs.csv')














# MCM5 and CDC45 IDs ------------------------------------------------------




# Load PLINK format SNP names for MGSvar SNPs
mcm5cdc45MGSsnps = readLines('/mnt/project/plink-snp-data-20230405/mcm5-cdc45-MGSvar-snp-ids.txt')
# 15 total
head(mcm5cdc45MGSsnps)
# need to remove "chr:pos:ref:alt"
mcm5cdc45MGSsnps <- mcm5cdc45MGSsnps[-1]

# Use data.table (faster) to load the raw data from PLINK
mcm5cdc45MGSraw = fread('/mnt/project/plink-snp-data-20230405/mcm5-cdc45-MGSvar-snp-data/mcm5-cdc45-MGSvar-snps.raw', data.table=FALSE)

class(mcm5cdc45MGSraw)
# "data.frame"
dim(mcm5cdc45MGSraw)
# 469835     20


# We only want the MGS SNPs
# First need to get the SNP id formats matching
# Create new variable with extra "_X" added:
mcm5cdc45MGSsnp_ids = strsplit(mcm5cdc45MGSsnps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(mcm5cdc45MGSsnp_ids)

# Which columns of the raw data match SNPs of interest?
match_ids_mcm5cdc45 = match(mcm5cdc45MGSsnp_ids, colnames(mcm5cdc45MGSraw))
match_ids_mcm5cdc45
# NA  7 18 17 13 20 NA  8  9 15 NA 16 10 19


#positions of extra variants 
all_names_mcm5cdc45 <- colnames(mcm5cdc45MGSraw)
all_names_mcm5cdc45
# "FID"                                 "IID"                                 "PAT"                                
# [4] "MAT"                                 "SEX"                                 "PHENOTYPE"                          
# [7] "22:19479969:A:C_C"                   "22:19482818:C:T_T"                   "22:19483983:A:G_G"                  
# [10] "22:19494594:G:A_A"                   "22:19494594:G:T_T"                   "22:19494594:GTCTCCCAGGTACATGCCA:G_G"
# [13] "22:19505448:C:A_A"                   "22:19507454:C:A_A"                   "22:19507770:C:A_A"                  
# [16] "22:19507830:C:T_T"                   "22:19514801:C:T_T"                   "22:19514996:C:T_T"                  
# [19] "22:19515024:C:T_T"                   "22:19518867:C:T_T"

setdiff(all_names_mcm5cdc45, mcm5cdc45MGSsnp_ids)
# "FID"                                 "IID"                                 "PAT"                                
# [4] "MAT"                                 "SEX"                                 "PHENOTYPE"                          
# [7] "22:19494594:G:T_T"                   "22:19494594:GTCTCCCAGGTACATGCCA:G_G" "22:19507454:C:A_A"

extra_variants_mcm5cdc45 <- setdiff(all_names_mcm5cdc45, mcm5cdc45MGSsnp_ids)
match(extra_variants_mcm5cdc45, colnames(mcm5cdc45MGSraw))
# 1  2  3  4  5  6 11 12 14
# 1:6 are meant to be there

mcm5cdc45MGSvar <- mcm5cdc45MGSraw[, -c(11, 12, 14)]
match(extra_variants_mcm5cdc45, colnames(mcm5cdc45MGSvar))
# 1  2  3  4  5  6 NA NA NA

# removing the first 6 columns for the next part of the analysis
mcm5cdc45genotypes <- mcm5cdc45MGSvar[, -c(1:6)]

# How many 0's, 1's and 2's are there across alls SNPs and samples?
# minor allele homoygotes
sum(mcm5cdc45genotypes == 2, na.rm=TRUE)
# 0

# heterozygotes
sum(mcm5cdc45genotypes == 1, na.rm=TRUE)
# 1629

# major allele homoygotes
sum(mcm5cdc45genotypes == 0, na.rm=TRUE)
# 5166401

# Get counts of 1's and 2's (i.e., any minor allele) per SNP
mcm5cdc45_allele_counts = colSums(mcm5cdc45genotypes!=0, na.rm=TRUE)
mcm5cdc45_allele_counts
# 22:19479969:A:C_C 22:19482818:C:T_T 22:19483983:A:G_G 22:19494594:G:A_A 22:19505448:C:A_A 22:19507770:C:A_A 22:19507830:C:T_T 
# 1                 7              1452                 7                 3                 2                99 
# 22:19514801:C:T_T 22:19514996:C:T_T 22:19515024:C:T_T 22:19518867:C:T_T 
# 10                19                15                14 




#
# separating by gene
#install.packages("readr")
#library(readr)
system('dx download Emily-folder/MGSvar/gene_MGSvar_csv/MCM5_CDC45_MGSvar.csv')
MCM5_CDC45_MGSvar <- read.csv("MCM5_CDC45_MGSvar.csv")
View(MCM5_CDC45_MGSvar)

# first making the variant IDs match
MCM5_CDC45_MGSvar$variants <- paste(MCM5_CDC45_MGSvar$chr, 
                                    MCM5_CDC45_MGSvar$pos,
                                    MCM5_CDC45_MGSvar$ref,
                                    MCM5_CDC45_MGSvar$alt,
                                    sep = ":")

MCM5_CDC45_MGSvar$variants <- paste(MCM5_CDC45_MGSvar$variants,
                                    MCM5_CDC45_MGSvar$alt,
                                    sep = "_")
MCM5_CDC45_MGSvar$variants
# looks good!

MCM5_variants <- MCM5_CDC45_MGSvar$variants[MCM5_CDC45_MGSvar$gene == "MCM5"]
MCM5_variants

CDC45_variants <- MCM5_CDC45_MGSvar$variants[MCM5_CDC45_MGSvar$gene == "CDC45"]
CDC45_variants

#df <- df[!(names(df) %in% rem)]
MCM5_MGSvar <- mcm5cdc45MGSvar[!(names(mcm5cdc45MGSvar) %in% CDC45_variants)]

CDC45_MGSvar <- mcm5cdc45MGSvar[!(names(mcm5cdc45MGSvar) %in% MCM5_variants)]

# to check:
colnames(MCM5_MGSvar)
# "FID"       "IID"       "PAT"       "MAT"       "SEX"       "PHENOTYPE"
MCM5_variants
# "22:35416388:C:T_T"

colnames(CDC45_MGSvar)
# "FID"               "IID"               "PAT"               "MAT"               "SEX"               "PHENOTYPE"        
# [7] "22:19479969:A:C_C" "22:19482818:C:T_T" "22:19483983:A:G_G" "22:19494594:G:A_A" "22:19505448:C:A_A" "22:19507770:C:A_A"
# [13] "22:19507830:C:T_T" "22:19514801:C:T_T" "22:19514996:C:T_T" "22:19515024:C:T_T" "22:19518867:C:T_T"
CDC45_variants
# "22:19479969:A:C_C" "22:19514996:C:T_T" "22:19514801:C:T_T" "22:19505448:C:A_A" "22:19518867:C:T_T" "22:19481044:A:G_G"
# [7] "22:19482818:C:T_T" "22:19483983:A:G_G" "22:19507770:C:A_A" "22:19507454:C:T_T" "22:19507830:C:T_T" "22:19494594:G:A_A"
# [13] "22:19515024:C:T_T"

# so there are no MCM5 variants in the biobank

###
# extracting list of patient IDs
# add a new column that contains the sum of SNP columns for each row(participant)
CDC45_MGSvar$Sums <- rowSums(CDC45_MGSvar[,7:17], na.rm = TRUE)

# remove the rows that have $Sum = 0
CDC45_MGSvar <- CDC45_MGSvar[CDC45_MGSvar$Sums !=0, ]
# 1629 obs

unique(CDC45_MGSvar[,18])
# 1

# make a list of participant IDs w relevant SNPs
CDC45_MGSvar_IDs <- CDC45_MGSvar[, 1]

CDC45_MGSvar_IDs <- as.character(CDC45_MGSvar_IDs)

# save
write(CDC45_MGSvar_IDs, "CDC45_MGSvar_IDs.txt")

# now to save this to UKB project
system('dx upload CDC45_MGSvar_IDs.txt --path Emily-folder/MGSvar/CDC45_MGSvar_IDs.txt')










### making MGSvar_IDs_VARs.csv
# no individuals have more than one of the varaiants
# but this may not be the case for other times I want to use this code

# make a copy to mess around with
CDC45_MGSvar_test1 <- CDC45_MGSvar

CDC45_MGSvar_test1$VAR <- names(CDC45_MGSvar_test1[,c(7:17)][max.col(CDC45_MGSvar_test1[,c(7:17)])])
# it works!!!
# but it adds a number to the end of some of the variants, not ideal
# this code ammends this:
CDC45_MGSvar_test1$VAR <- sub('\\..*', '', CDC45_MGSvar_test1$VAR)

# condense it down to just the ID number and the variant information
CDC45_MGSvar_test2 <- CDC45_MGSvar_test1[,c(1, 18, 19)]
# can save this now to annotate the entire cohort dataframe

#save it to R
write.csv(CDC45_MGSvar_test2, "CDC45_MGSvar_IDs_VARs.csv")
# now to save this to UKB project
system('dx upload CDC45_MGSvar_IDs_VARs.csv --path Emily-folder/MGSvar/CDC45_MGSvar_IDs_VARs.csv')
















# MCM7 IDs ----------------------------------------------------------------


# Load PLINK format SNP names for MGSvar SNPs
mcm7MGSsnps = readLines('/mnt/project/plink-snp-data-20230405/mcm7-MGSvar-snp-ids.txt')
head(mcm7MGSsnps)
# "chr:pos:ref:alt" "7:100099190:G:A" "7:100095450:T:C"

# Get rid of first entry ("chr:pos:ref:alt")
mcm7MGSsnps = mcm7MGSsnps[-1]
head(mcm7MGSsnps)

# Use data.table (faster) to load the raw data from PLINK
mcm7MGSraw = fread('/mnt/project/plink-snp-data-20230405/mcm7-MGSvar-snp-data/mcm7-MGSvar-snps.raw', data.table=FALSE)

#checking
class(mcm7MGSraw)
# "data.frame"
dim(mcm7MGSraw)
# 469835      8

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
mcm7MGSraw[1:5,1:8]

# We only want the MGS SNPs
# First need to get the SNP id formats matching
head(mcm7MGSsnps)
# "7:100099190:G:A" "7:100095450:T:C"
colnames(mcm7MGSraw)[7:10]
# "7:100095450:T:C_C" "7:100099190:G:A_A" NA                  NA 

# Create new variable with extra "_X" added:
mcm7MGSsnp_ids = strsplit(mcm7MGSsnps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(mcm7MGSsnp_ids)
# "7:100099190:G:A_A" "7:100095450:T:C_C"




## checking frequency data 


# removing the first 6 columns for the next part of the analysis
mcm7genotypes <- mcm7MGSraw[, -c(1:6)]

# How many 0's, 1's and 2's are there across alls SNPs and samples?
# minor allele homoygotes
sum(mcm7genotypes == 2, na.rm=TRUE)
# 0 (perfect)

# heterozygotes
sum(mcm7genotypes == 1, na.rm=TRUE)
# 54 (looks good!)

# major allele homoygotes
sum(mcm7genotypes == 0, na.rm=TRUE)
# 939615

# Get counts of 1's and 2's (i.e., any minor allele) per SNP
mcm7_allele_counts = colSums(mcm7genotypes!=0, na.rm=TRUE)
mcm7_allele_counts
# 7:100095450:T:C_C 7:100099190:G:A_A 
# 5                49
# looks perfect! perfectly matches the preliminary look!

# i don't think I need to check against plink - set of variants is small enough that i can see that everything is good :)









## removing variants that I don't care about
# because plink extracts all recorded variants at each location i'm interested in!

# want to duplicate this dataframe to mess around with it
mcm7MGSvar <- mcm7MGSraw

# Which columns of the raw data match SNPs of interest?
match_ids_mcm7 = match(mcm7MGSsnp_ids, colnames(mcm7MGSraw))
match_ids_mcm7
# 8 7

#positions of extra variants 
all_names_mcm7 <- colnames(mcm7MGSraw)
all_names_mcm7

setdiff(all_names_mcm7, mcm7MGSsnp_ids)
# "FID"       "IID"       "PAT"       "MAT"       "SEX"       "PHENOTYPE"

# no unwanted variants :)



###
# extracting list of patient IDs
# add a new column that contains the sum of SNP columns for each row(participant)
mcm7MGSvar$Sums <- rowSums(mcm7MGSvar[ , 7:8], na.rm = TRUE)
# remove the rows that have $Sum = 0
mcm7MGSvar <- mcm7MGSvar[mcm7MGSvar$Sums !=0, ]
# 54 obs. (perfect!)
unique(mcm7MGSvar[,9])
# 1 (each individual only has the 1 variant)

# make a list of participant IDs w relevant SNPs
MCM7_MGSvar_IDs <- mcm7MGSvar[, 1]

### figured out that its essential to change to character - or will save as chunks of entries
MCM7_MGSvar_IDs <- as.character(MCM7_MGSvar_IDs)
###

#save 'ORC1_MGSvar_IDs'
write(MCM7_MGSvar_IDs, "MCM7_MGSvar_IDs.txt")

# now to save this to UKB project
system('dx upload MCM7_MGSvar_IDs.txt --path Emily-folder/MGSvar/MCM7_MGSvar_IDs.txt')








### making MGSvar_IDs_VARs.csv
# working with donsonMGSvar
# no individuals have more than one of the varaiants
# but this may not be the case for other times I want to use this code

# make a copy to mess around with
mcm7MGSvar_test1 <- mcm7MGSvar

mcm7MGSvar_test1$VAR <- names(mcm7MGSvar_test1[,c(7:8)][max.col(mcm7MGSvar_test1[,c(7:8)])])
# it works!!!
# but it adds a number to the end of some of the variants, not ideal
# this code ammends this:
mcm7MGSvar_test1$VAR <- sub('\\..*', '', mcm7MGSvar_test1$VAR)

# condense it down to just the ID number and the variant information
mcm7MGSvar_test2 <- mcm7MGSvar_test1[,c(1, 9, 10)]
# can save this now to annotate the entire cohort dataframe

#save it to R
write.csv(mcm7MGSvar_test2, "MCM7_MGSvar_IDs_VARs.csv")
# now to save this to UKB project
system('dx upload MCM7_MGSvar_IDs_VARs.csv --path Emily-folder/MGSvar/MCM7_MGSvar_IDs_VARs.csv')



















# ORC4 IDs ----------------------------------------------------------------


# Load PLINK format SNP names for MGSvar SNPs
orc4MGSsnps = readLines('/mnt/project/plink-snp-data-20230405/orc4-MGSvar-snp-ids.txt')
head(orc4MGSsnps)
# "chr:pos:ref:alt" "2:147952440:T:C" "2:147948190:G:C" "2:147939142:T:C"

# Get rid of first entry ("chr:pos:ref:alt")
orc4MGSsnps = orc4MGSsnps[-1]
head(orc4MGSsnps)

# Use data.table (faster) to load the raw data from PLINK
orc4MGSraw = fread('/mnt/project/plink-snp-data-20230405/orc4-MGSvar-snp-data/orc4-MGSvar-snps.raw', data.table=FALSE)

#checking
class(orc4MGSraw)
# "data.frame"
dim(orc4MGSraw)
# 469835      7

# Take a look at the first 5 rows, and first 7 columns
# First 6 columns are sample info
orc4MGSraw[1:5,1:7]

# We only want the MGS SNPs
# First need to get the SNP id formats matching
head(orc4MGSsnps)
# "2:147952440:T:C" "2:147948190:G:C" "2:147939142:T:C"
colnames(orc4MGSraw)[7:10]
# "2:147948190:G:C_C" NA                  NA                  NA

# Create new variable with extra "_X" added:
orc4MGSsnp_ids = strsplit(orc4MGSsnps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(orc4MGSsnp_ids)




## checking frequency data 


# removing the first 6 columns for the next part of the analysis
orc4genotypes <- orc4MGSraw[, -c(1:6)]

# How many 0's, 1's and 2's are there across alls SNPs and samples?
# minor allele homoygotes
sum(orc4genotypes == 2, na.rm=TRUE)
# 0 (perfect)

# heterozygotes
sum(orc4genotypes == 1, na.rm=TRUE)
# 3 (looks good!)

# major allele homoygotes
sum(orc4genotypes == 0, na.rm=TRUE)
# 469813


# i don't think I need to check against plink - set of variants is small enough that i can see that everything is good :)









## removing variants that I don't care about
# because plink extracts all recorded variants at each location i'm interested in!

# want to duplicate this dataframe to mess around with it
orc4MGSvar <- orc4MGSraw

# Which columns of the raw data match SNPs of interest?
match_ids_orc4 = match(orc4MGSsnp_ids, colnames(orc4MGSraw))
match_ids_orc4
# NA  7 NA

#positions of extra variants 
all_names_orc4 <- colnames(orc4MGSraw)
all_names_orc4
# "FID"               "IID"               "PAT"               "MAT"               "SEX"               "PHENOTYPE"        
# [7] "2:147948190:G:C_C"

setdiff(all_names_orc4, orc4MGSsnp_ids)
# "FID"       "IID"       "PAT"       "MAT"       "SEX"       "PHENOTYPE"

# no unwanted variants :)



###
# extracting list of patient IDs

# the method I had previously used relied on summing multiple columns in a row
# but there is only one variant in each of these dfs
# so can just duplicate the column and rename
orc4MGSvar$Sums <- orc4MGSvar$`2:147948190:G:C_C`


# remove the rows that have $Sum = 0
orc4MGSvar <- orc4MGSvar[orc4MGSvar$Sums !=0, ]
# 22 obs
# also need to remove NA values
orc4MGSvar <- orc4MGSvar[!is.na(orc4MGSvar$Sums),]
# 3 obs



# make a list of participant IDs w relevant SNPs
ORC4_MGSvar_IDs <- orc4MGSvar[, 1]

### figured out that its essential to change to character - or will save as chunks of entries
ORC4_MGSvar_IDs <- as.character(ORC4_MGSvar_IDs)
###

#save 'ORC1_MGSvar_IDs'
write(ORC4_MGSvar_IDs, "ORC4_MGSvar_IDs.txt")

# now to save this to UKB project
system('dx upload ORC4_MGSvar_IDs.txt --path Emily-folder/MGSvar/ORC4_MGSvar_IDs.txt')






# only has one variant, so instead of using the method I'd used for DONSON,
# I can instead just paste the variant name into a VAR column!

orc4MGSvar$VAR <- "2:147948190:G:C_C"


# condense it down to just the ID number and the variant information
orc4MGSvar_test2 <- orc4MGSvar[,c(1, 8, 9)]
# can save this now to annotate the entire cohort dataframe


#save it to R
write.csv(orc4MGSvar_test2, "ORC4_MGSvar_IDs_VARs.csv")
# now to save this to UKB project
system('dx upload ORC4_MGSvar_IDs_VARs.csv --path Emily-folder/MGSvar/ORC4_MGSvar_IDs_VARs.csv')

















## ORC1 / CDT1 / ORC6 were in a different script, need to add them here as well
# use MGSvar_script_END.R as a reference



# ORC1 IDs ----------------------------------------------------------------


# Load PLINK format SNP names for MGSvar SNPs
orc1MGSsnps = readLines('/mnt/project/plink-MGSvar-snp-data-20230217/orc1-MGSvar-snp-ids.txt')
head(orc1MGSsnps)
# "chr:pos:ref:alt" "1:52397773:C:T"  "1:52385264:T:C"  "1:52383437:G:A"  "1:52397707:T:C"  "1:52397821:A:G" 

# Get rid of first entry ("chr:pos:ref:alt")
orc1MGSsnps = orc1MGSsnps[-1]
head(orc1MGSsnps)
#"1:52397773:C:T" "1:52385264:T:C" "1:52383437:G:A" "1:52397707:T:C" "1:52397821:A:G" "1:52375574:C:T"

# Use data.table (faster) to load the raw data from PLINK
orc1MGSraw = fread('/mnt/project/plink-MGSvar-snp-data-20230217/orc1-MGSvar-snp-data/orc1-MGSvar-snps.raw', data.table=FALSE)

class(orc1MGSraw)
# "data.frame"
dim(orc1MGSraw)
# 469835     15

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
orc1MGSraw[1:5,1:8]

# We only want the MGS SNPs
# First need to get the SNP id formats matching
head(orc1MGSsnps)
# "1:52397773:C:T" "1:52385264:T:C" "1:52383437:G:A" "1:52397707:T:C" "1:52397821:A:G" "1:52375574:C:T"
colnames(orc1MGSraw)[7:10]
# "1:52375441:G:A_A" "1:52375574:C:T_T" "1:52383437:G:A_A" "1:52384584:G:A_A"

# Create new variable with extra "_X" added:
orc1MGSsnp_ids = strsplit(orc1MGSsnps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(orc1MGSsnp_ids)
# "1:52397773:C:T_T" "1:52385264:T:C_C" "1:52383437:G:A_A" "1:52397707:T:C_C" "1:52397821:A:G_G" "1:52375574:C:T_T"

# Which columns of the raw data match SNPs of interest?
match_ids_orc1 = match(orc1MGSsnp_ids, colnames(orc1MGSraw))
match_ids_orc1
# 12 11  9 NA 13  8 10  7 15 NA

#positions of extra variants 
all_names_orc1 <- colnames(orc1MGSraw)
all_names_orc1

setdiff(all_names_orc1, orc1MGSsnp_ids)
# "FID"              "IID"              "PAT"              "MAT"              "SEX"              "PHENOTYPE"       
# "1:52401368:C:G_G"

extra_variants_orc1 <- setdiff(all_names_orc1, orc1MGSsnp_ids)
match(extra_variants_orc1, colnames(orc1MGSraw))
# 1  2  3  4  5  6 14
# 1:6 are all meant to be there!

orc1MGSraw <- orc1MGSraw[, -14]
match(extra_variants_orc1, colnames(orc1MGSraw))





# removing the first 6 columns for the next part of the analysis
orc1genotypes <- orc1MGSraw[, -c(1:6)]

# How many 0's, 1's and 2's are there across alls SNPs and samples?
# minor allele homoygotes
sum(orc1genotypes == 2, na.rm=TRUE)
# 0 (perfect)

# heterozygotes
sum(orc1genotypes == 1, na.rm=TRUE)
# 654 (looks good!)

# major allele homoygotes
sum(orc1genotypes == 0, na.rm=TRUE)
# 3757872

# Get counts of 1's and 2's (i.e., any minor allele) per SNP
orc1_allele_counts = colSums(orc1genotypes!=0, na.rm=TRUE)
orc1_allele_counts
# 1:52375441:G:A_A 1:52375574:C:T_T 1:52383437:G:A_A 1:52384584:G:A_A 1:52385264:T:C_C 1:52397773:C:T_T 1:52397821:A:G_G 1:52401368:C:T_T 
# 14              195               16                9               13              399                1                7 
# looks perfect! perfectly matches the preliminary look!

# i don't think I need to check against plink - set of variants is small enough that i can see that everything is good :)


###
orc1MGSvar <- orc1MGSraw
# extracting list of patient IDs
# add a new column that contains the sum of SNP columns for each row(participant)
orc1MGSvar$Sums <- rowSums(orc1MGSvar[ , 7:14], na.rm = TRUE)
# remove the rows that have $Sum = 0
orc1MGSvar <- orc1MGSvar[orc1MGSvar$Sums !=0, ]
# 654 obs. (perfect!)
unique(orc1MGSvar[,15])
# 1 (each individual only has the 1 variant)

# make a list of participant IDs w relevant SNPs
ORC1_MGSvar_IDs <- orc1MGSvar[, 1]

### figured out that its essential to change to character - or will save as chunks of entries
ORC1_MGSvar_IDs <- as.character(ORC1_MGSvar_IDs)
###

#save 'ORC1_MGSvar_IDs'
write(ORC1_MGSvar_IDs, "ORC1_MGSvar_IDs.txt")

# now to save this to UKB project
system('dx upload ORC1_MGSvar_IDs.txt --path Emily-folder/MGSvar/ORC1_MGSvar_IDs.txt')











### making MGSvar_IDs_VARs.csv
# no individuals have more than one of the varaiants
# but this may not be the case for other times I want to use this code

# make a copy to mess around with
orc1MGSvar_test1 <- orc1MGSvar

# moving the var names into a new column when a participant has the variant
orc1MGSvar_test1$VAR <- names(orc1MGSvar_test1[,c(7:14)][max.col(orc1MGSvar_test1[,c(7:14)])])
# for some reason this doesn't work, but it does work for all the other genes
# I can't tell any difference between this and the others!
summary(orc1MGSvar_test1)
# there is an NA value in the 1:52385264:T:C_C column
# replace NA with 0
d[is.na(d)] <- 0
orc1MGSvar_test1$`1:52385264:T:C_C`[is.na(orc1MGSvar_test1$`1:52385264:T:C_C`)] <- 0
orc1MGSvar_test1$`1:52385264:T:C_C` <- as.integer(orc1MGSvar_test1$`1:52385264:T:C_C`)
# now it should work
orc1MGSvar_test1$VAR <- names(orc1MGSvar_test1[,c(7:14)][max.col(orc1MGSvar_test1[,c(7:14)])])
# it works!!!
# but it adds a number to the end of some of the variants, not ideal
# this code ammends this:
orc1MGSvar_test1$VAR <- sub('\\..*', '', orc1MGSvar_test1$VAR)

# condense it down to just the ID number and the variant information
orc1MGSvar_test2 <- orc1MGSvar_test1[,c(1, 15, 16)]
# can save this now to annotate the entire cohort dataframe

#save it to R
write.csv(orc1MGSvar_test2, "ORC1_MGSvar_IDs_VARs.csv")
# now to save this to UKB project
system('dx upload ORC1_MGSvar_IDs_VARs.csv --path Emily-folder/MGSvar/ORC1_MGSvar_IDs_VARs.csv')















# CDT1 and ORC6 IDs -------------------------------------------------------


# Load PLINK format SNP names for MGSvar SNPs
cdt1orc6MGSsnps = readLines('/mnt/project/plink-MGSvar-snp-data-20230217b/cdt1-orc6-MGSvar-snp-ids.txt')
# 16 total
head(cdt1orc6MGSsnps)
# need to remove "chr:pos:ref:alt"
cdt1orc6MGSsnps <- cdt1orc6MGSsnps[-1]

# Use data.table (faster) to load the raw data from PLINK
cdt1orc6MGSraw = fread('/mnt/project/plink-MGSvar-snp-data-20230217b/cdt1-orc6-MGSvar-snp-data/cdt1-orc6-MGSvar-snps.raw', data.table=FALSE)

class(cdt1orc6MGSraw)
# "data.frame"
dim(cdt1orc6MGSraw)
# 469835     23



# We only want the MGS SNPs
# First need to get the SNP id formats matching
# Create new variable with extra "_X" added:
cdt1orc6MGSsnp_ids = strsplit(cdt1orc6MGSsnps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(cdt1orc6MGSsnp_ids)

# Which columns of the raw data match SNPs of interest?
match_ids_cdt1orc6 = match(cdt1orc6MGSsnp_ids, colnames(cdt1orc6MGSraw))
match_ids_cdt1orc6
# NA  7 11  9 10 NA 19 22 12 14 NA 18 20 15 NA


#positions of extra variants 
all_names_cdt1orc6 <- colnames(cdt1orc6MGSraw)
all_names_cdt1orc6
# [1] "FID"               "IID"               "PAT"               "MAT"               "SEX"               "PHENOTYPE"        
# [7] "16:46689707:T:C_C" "16:46689707:T:G_G" "16:46690992:A:G_G" "16:46690996:C:T_T" "16:46693187:G:A_A" "16:88804027:G:A_A"
# [13] "16:88804027:G:C_C" "16:88804667:G:C_C" "16:88805869:G:T_T" "16:88807362:C:A_A" "16:88807362:C:G_G" "16:88807362:C:T_T"
# [19] "16:88807390:G:A_A" "16:88807407:G:A_A" "16:88807407:G:T_T" "16:88808197:C:A_A" "16:88808197:C:T_T"

setdiff(all_names_cdt1orc6, cdt1orc6MGSsnp_ids)
# [1] "FID"               "IID"               "PAT"               "MAT"               "SEX"               "PHENOTYPE"        
# [7] "16:46689707:T:G_G" "16:88804027:G:C_C" "16:88807362:C:A_A" "16:88807362:C:G_G" "16:88807407:G:T_T" "16:88808197:C:T_T"

extra_variants_cdt1orc6 <- setdiff(all_names_cdt1orc6, cdt1orc6MGSsnp_ids)
match(extra_variants_cdt1orc6, colnames(cdt1orc6MGSraw))
# 1  2  3  4  5  6  8 13 16 17 21 23
# 1:6 are meant to be there

cdt1orc6MGSraw <- cdt1orc6MGSraw[, -c(8, 13, 16, 17, 21, 23)]
match(extra_variants_cdt1orc6, colnames(cdt1orc6MGSraw))
# 1  2  3  4  5  6 NA NA NA NA NA NA




# removing the first 6 columns for the next part of the analysis
cdt1orc6genotypes <- cdt1orc6MGSraw[, -c(1:6)]

# How many 0's, 1's and 2's are there across alls SNPs and samples?
# minor allele homoygotes
sum(cdt1orc6genotypes == 2, na.rm=TRUE)
# 4 !!!

# heterozygotes
sum(cdt1orc6genotypes == 1, na.rm=TRUE)
# 1520 (looks good!)

# major allele homoygotes
sum(cdt1orc6genotypes == 0, na.rm=TRUE)
# 5166454

# Get counts of 1's and 2's (i.e., any minor allele) per SNP
cdt1orc6_allele_counts = colSums(cdt1orc6genotypes!=0, na.rm=TRUE)
cdt1orc6_allele_counts
# 16:46689707:T:C_C 16:46690992:A:G_G 16:46690996:C:T_T 16:46693187:G:A_A 16:88804027:G:A_A 16:88804667:G:C_C 16:88805869:G:T_T 16:88807362:C:T_T 
# 660                 2                 1               207                20               107                26               247 
# 16:88807390:G:A_A 16:88807407:G:A_A 16:88808197:C:A_A 
# 145                87                22 
# great






# separating by gene
#install.packages("readr")
library(readr)
system('dx download Emily-folder/MGSvar/gene_MGSvar_csv/CDT1_ORC6_MGSvar.csv')
CDT1_ORC6_MGSvar <- read.csv("CDT1_ORC6_MGSvar.csv")
View(CDT1_ORC6_MGSvar)

# first making the variant IDs match
CDT1_ORC6_MGSvar$variants <- paste(CDT1_ORC6_MGSvar$chr, 
                                   CDT1_ORC6_MGSvar$pos,
                                   CDT1_ORC6_MGSvar$ref,
                                   CDT1_ORC6_MGSvar$alt,
                                   sep = ":")

CDT1_ORC6_MGSvar$variants <- paste(CDT1_ORC6_MGSvar$variants,
                                   CDT1_ORC6_MGSvar$alt,
                                   sep = "_")
CDT1_ORC6_MGSvar$variants
# looks good!

ORC6_variants <- CDT1_ORC6_MGSvar$variants[CDT1_ORC6_MGSvar$gene == "ORC6"]
ORC6_variants

CDT1_variants <- CDT1_ORC6_MGSvar$variants[CDT1_ORC6_MGSvar$gene == "CDT1"]
CDT1_variants

#df <- df[!(names(df) %in% rem)]
ORC6_MGSvar <- cdt1orc6MGSraw[!(names(cdt1orc6MGSraw) %in% CDT1_variants)]

CDT1_MGSvar <- cdt1orc6MGSraw[!(names(cdt1orc6MGSraw) %in% ORC6_variants)]

# to check:
colnames(CDT1_MGSvar)
#"16:88804027:G:A_A" "16:88804667:G:C_C" "16:88805869:G:T_T" "16:88807362:C:T_T" 
#"16:88807390:G:A_A" "16:88807407:G:A_A" "16:88808197:C:A_A"
CDT1_variants
#"16:88807390:G:A_A" "16:88808197:C:A_A" "16:88804027:G:A_A" "16:88804667:G:C_C" 
#"16:88806633:C:T_T" "16:88807362:C:T_T" "16:88807407:G:A_A" "16:88805869:G:T_T" 
#"16:88807257:A:G_G"
# looks good!





###
# extracting list of patient IDs
# add a new column that contains the sum of SNP columns for each row(participant)
CDT1_MGSvar$Sums <- rowSums(CDT1_MGSvar[,7:13], na.rm = TRUE)
ORC6_MGSvar$Sums <- rowSums(ORC6_MGSvar[,7:10], na.rm = TRUE)

# remove the rows that have $Sum = 0
CDT1_MGSvar <- CDT1_MGSvar[CDT1_MGSvar$Sums !=0, ]
# 653obs
ORC6_MGSvar <- ORC6_MGSvar[ORC6_MGSvar$Sums !=0, ]
# 870obs
unique(CDT1_MGSvar[,14])
# 1 2 (they either have 1 or 2 vars)
unique(ORC6_MGSvar[,11])
# 1 (no one has more than 1 ORC6 var)

# make a list of participant IDs w relevant SNPs
CDT1_MGSvar_IDs <- CDT1_MGSvar[, 1]
ORC6_MGSvar_IDs <- ORC6_MGSvar[, 1]

CDT1_MGSvar_IDs <- as.character(CDT1_MGSvar_IDs)
ORC6_MGSvar_IDs <- as.character(ORC6_MGSvar_IDs)

# save
write(CDT1_MGSvar_IDs, "CDT1_MGSvar_IDs.txt")
write(ORC6_MGSvar_IDs, "ORC6_MGSvar_IDs.txt")

# now to save this to UKB project
system('dx upload CDT1_MGSvar_IDs.txt --path Emily-folder/MGSvar/CDT1_MGSvar_IDs.txt')
system('dx upload ORC6_MGSvar_IDs.txt --path Emily-folder/MGSvar/ORC6_MGSvar_IDs.txt')







### making MGSvar_IDs_VARs.csv
# this may be different to the other times, because at least one individual has more than one CDT1 variant
# unsure currently whether this is only homozygous or heterozygous as well!

# starting with ORC6
# make a copy to mess around with
ORC6_MGSvar_test1 <- ORC6_MGSvar

ORC6_MGSvar_test1$VAR <- names(ORC6_MGSvar_test1[,c(7:10)][max.col(ORC6_MGSvar_test1[,c(7:10)])])
# it works!!!
# but it adds a number to the end of some of the variants, not ideal
# this code ammends this:
ORC6_MGSvar_test1$VAR <- sub('\\..*', '', ORC6_MGSvar_test1$VAR)

# condense it down to just the ID number and the variant information
ORC6_MGSvar_test2 <- ORC6_MGSvar_test1[,c(1, 11, 12)]
# can save this now to annotate the entire cohort dataframe

#save it to R
write.csv(ORC6_MGSvar_test2, "ORC6_MGSvar_IDs_VARs.csv")
# now to save this to UKB project
system('dx upload ORC6_MGSvar_IDs_VARs.csv --path Emily-folder/MGSvar/ORC6_MGSvar_IDs_VARs.csv')






# now to do CDT1 (which has some 2 values in Sum)
# have an investigate of these columns
CDT1_MGSvar[CDT1_MGSvar$Sums >1,]
# 4 individuals have homozygous 16:88807362:C:T_T
# 1 individual has comp het 16:88807362:C:T_T and 16:88808197:C:A_A

# this makes the next bit of code a little more tricky
# someone has two variants associated! which means they need 2 VAR values
# I could do this manually but that won't work when the datasets are much larger
# for now I'm going to separate by ';'


# make a copy to mess around with
CDT1_MGSvar_test1 <- CDT1_MGSvar

CDT1_MGSvar_test1$VAR <- names(CDT1_MGSvar_test1[,c(7:10)][max.col(CDT1_MGSvar_test1[,c(7:10)])])
# checking the 2 Sum people
CDT1_MGSvar_test1[CDT1_MGSvar_test1$Sums >1,]
# for the comp het person they only have the 16:88807362:C:T_T valiant here

# it works!!!
# but it adds a number to the end of some of the variants, not ideal
# this code ammends this:
CDT1_MGSvar_test1$VAR <- sub('\\..*', '', CDT1_MGSvar_test1$VAR)


## now to add the extra variant for the comp het person
CDT1_MGSvar_test1$VAR[CDT1_MGSvar_test1$`16:88808197:C:A_A`==1 & CDT1_MGSvar_test1$`16:88807362:C:T_T`==1] <- "16:88807362:C:T_T;16:88808197:C:A_A"
# need to figure out a more automated way of doing this for larger datasets!


# condense it down to just the ID number and the variant information
CDT1_MGSvar_test2 <- CDT1_MGSvar_test1[,c(1, 14, 15)]
# can save this now to annotate the entire cohort dataframe

#save it to R
write.csv(CDT1_MGSvar_test2, "CDT1_MGSvar_IDs_VARs.csv")
# now to save this to UKB project
system('dx upload CDT1_MGSvar_IDs_VARs.csv --path Emily-folder/MGSvar/CDT1_MGSvar_IDs_VARs.csv')







system('dx upload MGSvar_script3_END.R --path Emily-folder/MGSvar/MGSvar_script3_END.R')