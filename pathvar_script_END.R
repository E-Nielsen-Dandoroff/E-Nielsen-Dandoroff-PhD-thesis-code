# 2024 02 20

# particiapant ID extraction for path data.

# using Emily-folder/20230110-O1C1O6-BayesDel/CDT1_BayesDel_script_END.R as a reference

# instance type: mem1_hdd1_v2_x16 = 31.3GB memory, 1600GB storage, 16 cores
# lots of memory required for loading in the large files!

# Load packages
library(data.table)
library(dplyr)

# set things up so that I only need to execute this code the one time







# chr 1 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for Chr1 snps
chr1snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr1-snp-ids.txt')
head(chr1snps)
# output as expected: "chr:pos:ref:alt" "1:109621036:T:C" "1:109621102:G:T" "1:109621233:C:G" "1:109621233:C:G" "1:109625333:T:A"

# Get rid of first entry ("chr:pos:ref:alt")
chr1snps = chr1snps[-1]
head(chr1snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr1raw_head = fread('/mnt/project/plink-snp-data-20240108-path/chr1-snp-data/chr1-snps.raw', data.table=FALSE)
# this gives the following warning:
# Warning message:
#   In fread("/mnt/project/plink-snp-data-20240108-path/chr1-snp-data/chr1-snps.raw",  :
#              Stopped early on line 9502. Expected 7074 fields but found 7870. Consider fill=TRUE and comment.char=. First discarded non-empty line: <<3979939 3979939 0 0 1 -9 0 0...

#the following site shows that others are having the same issue with no obvious solution
# https://github.com/Rdatatable/data.table/issues/2727

# loading in the rest of the data - except for the one problem line
#9500 7074
chr1raw_tail = fread('/mnt/project/plink-snp-data-20240108-path/chr1-snp-data/chr1-snps.raw', skip = 9502, data.table=FALSE)

dim(chr1raw_tail)
# 460334   7074



## merge chr1raw_head with chr1raw_tail
# rbind to add new obeservations
# requires common column names between the two datasets
colnames(chr1raw_tail) <- colnames(chr1raw_head)
# check usin identical
identical(names(chr1raw_tail), names(chr1raw_head))
# TRUE

chr1raw <- rbind(chr1raw_head, chr1raw_tail)

# process this in the same way you processed the MGSvar ones:

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr1raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr1snps)
# output: "1:109621036:T:C" "1:109621102:G:T" "1:109621233:C:G" "1:109621233:C:G" "1:109625333:T:A" "1:109625344:C:T"
colnames(chr1raw)[7:10]
# output: "1:1232291:C:T_T" "1:1232292:G:C_C" "1:1232301:G:A_A" "1:1232309:C:T_T"

# Create new variable with extra "_X" added:
chr1snp_ids <- strsplit(chr1snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr1snp_ids)
# output: "1:109621036:T:C_C" "1:109621102:G:T_T" "1:109621233:C:G_G" "1:109621233:C:G_G" "1:109625333:T:A_A" "1:109625344:C:T_T"

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr1 <- match(chr1snp_ids, colnames(chr1raw))
# 6619, so there are 455 that I need to remove

# Keep just those columns
#chr1genotypes = chr1raw[, match_ids_chr1]
#Error in `[.data.frame`(chr1raw, , match_ids_chr1) : 
#undefined columns selected


chr1colnames <- colnames(chr1raw)

# this still works - do this instead
extra_variants_chr1 <- setdiff(chr1colnames, chr1snp_ids)
indexchr1 <- match(extra_variants_chr1, colnames(chr1raw))
indexchr1 <- indexchr1[-c(1:6)]

# retain variants of interest by index
chr1raw <- chr1raw[, -indexchr1]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr1raw$Sums <- rowSums(chr1raw[, 7:6535], na.rm = TRUE)
# check
unique(chr1raw$Sums)
# 0 1 2 3 4

# remove the rows that have $Sum = 0
chr1raw <- chr1raw[chr1raw$Sums !=0, ]
# obs went from 469834 to 60203
unique(chr1raw[,6536])
# 1 2 3 4


# make two dfs to merge together later
chr1vars <- chr1raw[, c(7:6535)]
chr1ids <- chr1raw[, c(1, 6536)]

# 6531 columns now - the last is still 'Sums'

# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr1vars != 0, arr.ind = TRUE)
chr1vars[w] <- names(chr1vars)[w[,"col"]]

unique(chr1vars$`1:1232301:G:A_A`)
unique(chr1vars$`1:1232292:G:C_C`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr1vars[chr1vars == "0"] <- NA
chr1vars <- unite(chr1vars, "VAR", c(1:6529), sep = ";", na.rm = TRUE)

#check
unique(chr1vars)

# now to cbind with the ids df
chr1_final <- cbind(chr1ids, chr1vars)

# check to make sure things line up
chr1_final$VAR[10]
# "1:21884839:A:C_C"
chr1raw$`1:21884839:A:C_C`[10]
# 1
# matches!!!

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr1_final, "chr1_final.csv")
# now to save this to UKB project
system('dx upload chr1_final.csv --path Emily-folder/pathvar/chr1_final.csv')







# chr 2 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr2snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr2-snp-ids.txt')
head(chr2snps)
# output as expected: "chr:pos:ref:alt" "2:177392895:C:T" "2:177392904:C:T" "2:177392916:G:A" "2:177392923:T:C" "2:177392970:C:T"

# Get rid of first entry ("chr:pos:ref:alt")
chr2snps = chr2snps[-1]
head(chr2snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr2raw = fread('/mnt/project/plink-snp-data-20240108-path/chr2-snp-data/chr2-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr2raw)
# 469835   5095

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr2raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr2snps)
# output: "2:177392895:C:T" "2:177392904:C:T" "2:177392916:G:A" "2:177392923:T:C" "2:177392970:C:T" "2:177393040:C:T"
colnames(chr2raw)[7:10]
# output: "2:5692765:C:A_A" "2:5692768:G:A_A" "2:5692770:G:C_C" "2:5692771:A:G_G"

# Create new variable with extra "_X" added:
chr2snp_ids <- strsplit(chr2snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr2snp_ids)
# output: "2:177392895:C:T_T" "2:177392904:C:T_T" "2:177392916:G:A_A" "2:177392923:T:C_C" "2:177392970:C:T_T" "2:177393040:C:T_T"

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr2 <- match(chr2snp_ids, colnames(chr2raw))
# 4768, so there are 327 that I need to remove



chr2colnames <- colnames(chr2raw)

# this still works - do this instead
extra_variants_chr2 <- setdiff(chr2colnames, chr2snp_ids)
indexchr2 <- match(extra_variants_chr2, colnames(chr2raw))
indexchr2 <- indexchr2[-c(1:6)]

# retain variants of interest by index
chr2raw <- chr2raw[, -indexchr2]
#4712 now







### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr2raw$Sums <- rowSums(chr2raw[, 7:4712], na.rm = TRUE)
# check
unique(chr2raw$Sums)
# 0 1 2 3 4 5

# remove the rows that have $Sum = 0
chr2raw <- chr2raw[chr2raw$Sums !=0, ]
# obs went from 
unique(chr2raw[, 4713])
# 1 2 3 4 5


# make two dfs to merge together later
chr2vars <- chr2raw[, c(7:4712)]
chr2ids <- chr2raw[, c(1, 4713)]

# ### columns now - the last is still 'Sums'

# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr2vars != 0, arr.ind = TRUE)
chr2vars[w] <- names(chr2vars)[w[,"col"]]

unique(chr2vars$`2:5692765:C:A_A`)
unique(chr2vars$`2:5692786:A:C_C`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr2vars[chr2vars == "0"] <- NA
chr2vars <- unite(chr2vars, "VAR", c(1:4706), sep = ";", na.rm = TRUE)

#check
unique(chr2vars)

# now to cbind with the ids df
chr2_final <- cbind(chr2ids, chr2vars)

# check to make sure things line up
chr2_final$VAR[10]
# "2:215380924:T:A_A"
chr2raw$`2:215380924:T:A_A`[10]
# 1
# looks good

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr2_final, "chr2_final.csv")
# now to save this to UKB project
system('dx upload chr2_final.csv --path Emily-folder/pathvar/chr2_final.csv')
















# chr 3 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr3snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr3-snp-ids.txt')
head(chr3snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr3snps = chr3snps[-1]
head(chr3snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr3raw = fread('/mnt/project/plink-snp-data-20240108-path/chr3-snp-data/chr3-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr3raw)
# 469835   4593


# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr3raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr3snps)
# output: 
colnames(chr3raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr3snp_ids <- strsplit(chr3snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr3snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr3 <- match(chr3snp_ids, colnames(chr3raw))
# 4375, so there are 218 that I need to remove



chr3colnames <- colnames(chr3raw)

# this still works - do this instead
extra_variants_chr3 <- setdiff(chr3colnames, chr3snp_ids)
indexchr3 <- match(extra_variants_chr3, colnames(chr3raw))
indexchr3 <- indexchr3[-c(1:6)]

# retain variants of interest by index
chr3raw <- chr3raw[, -indexchr3]
# 4284 now






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr3raw$Sums <- rowSums(chr3raw[, 7:4284], na.rm = TRUE)
# check
unique(chr3raw$Sums)
# 0 1 2 3 4

# remove the rows that have $Sum = 0
chr3raw <- chr3raw[chr3raw$Sums !=0, ]
# obs went from 469835 to 41258
unique(chr3raw[, 4285])
# 1 2 3 4


# make two dfs to merge together later
chr3vars <- chr3raw[, c(7: 4284)]
chr3ids <- chr3raw[, c(1, 4285)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr3vars != 0, arr.ind = TRUE)
chr3vars[w] <- names(chr3vars)[w[,"col"]]

unique(chr3vars$`3:12489811:C:A_A`)
unique(chr3vars$`3:12489862:C:T_T`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr3vars[chr3vars == "0"] <- NA
chr3vars <- unite(chr3vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr3vars)

# now to cbind with the ids df
chr3_final <- cbind(chr3ids, chr3vars)

# check to make sure things line up
chr3_final$VAR[10]
# "3:58124392:G:C_C"
chr3raw$`3:58124392:G:C_C`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr3_final, "chr3_final.csv")
# now to save this to UKB project
system('dx upload chr3_final.csv --path Emily-folder/pathvar/chr3_final.csv')













# chr 4 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr4snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr4-snp-ids.txt')
head(chr4snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr4snps = chr4snps[-1]
head(chr4snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr4raw = fread('/mnt/project/plink-snp-data-20240108-path/chr4-snp-data/chr4-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr4raw)
# 469835   3179

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr4raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr4snps)
# output: 
colnames(chr4raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr4snp_ids <- strsplit(chr4snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr4snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr4 <- match(chr4snp_ids, colnames(chr4raw))




chr4colnames <- colnames(chr4raw)

# this still works - do this instead
extra_variants_chr4 <- setdiff(chr4colnames, chr4snp_ids)
indexchr4 <- match(extra_variants_chr4, colnames(chr4raw))
indexchr4 <- indexchr4[-c(1:6)]

# retain variants of interest by index
chr4raw <- chr4raw[, -indexchr4]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr4raw$Sums <- rowSums(chr4raw[, 7:2941], na.rm = TRUE)
# check
unique(chr4raw$Sums)
# 0 2 1 3 4

# remove the rows that have $Sum = 0
chr4raw <- chr4raw[chr4raw$Sums !=0, ]
# obs went from 469835 to 26257
unique(chr4raw[, 2942])
# 2 1 3 4


# make two dfs to merge together later
chr4vars <- chr4raw[, c(7:2941)]
chr4ids <- chr4raw[, c(1, 2942)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr4vars != 0, arr.ind = TRUE)
chr4vars[w] <- names(chr4vars)[w[,"col"]]

unique(chr4vars$`4:1804392:G:A_A`)
unique(chr4vars$`4:5562952:A:G_G`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr4vars[chr4vars == "0"] <- NA
chr4vars <- unite(chr4vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr4vars)

# now to cbind with the ids df
chr4_final <- cbind(chr4ids, chr4vars)

# check to make sure things line up
chr4_final$VAR[10]
# "4:95125072:T:C_C"
chr4raw$`4:95125072:T:C_C`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr4_final, "chr4_final.csv")
# now to save this to UKB project
system('dx upload chr4_final.csv --path Emily-folder/pathvar/chr4_final.csv')












# chr 5 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr5snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr5-snp-ids.txt')
head(chr5snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr5snps = chr5snps[-1]
head(chr5snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr5raw = fread('/mnt/project/plink-snp-data-20240108-path/chr5-snp-data/chr5-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr5raw)
# 469835   2736

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr5raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr5snps)
# output: 
colnames(chr5raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr5snp_ids <- strsplit(chr5snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr5snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr5 <- match(chr5snp_ids, colnames(chr5raw))
# 2536, so there are 200 that I need to remove



chr5colnames <- colnames(chr5raw)

# this still works - do this instead
extra_variants_chr5 <- setdiff(chr5colnames, chr5snp_ids)
indexchr5 <- match(extra_variants_chr5, colnames(chr5raw))
indexchr5 <- indexchr5[-c(1:6)]

# retain variants of interest by index
chr5raw <- chr5raw[, -indexchr5]
# 2519 now





### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr5raw$Sums <- rowSums(chr5raw[, 7:2519], na.rm = TRUE)
# check
unique(chr5raw$Sums)
# 1 0 2 3

# remove the rows that have $Sum = 0
chr5raw <- chr5raw[chr5raw$Sums !=0, ]
# obs went from 469835 to 28306
unique(chr5raw[, 2520])
# 1 2 3


# make two dfs to merge together later
chr5vars <- chr5raw[, c(7:2519)]
chr5ids <- chr5raw[, c(1, 2520)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr5vars != 0, arr.ind = TRUE)
chr5vars[w] <- names(chr5vars)[w[,"col"]]

unique(chr5vars$`5:36953706:G:C_C`)
unique(chr5vars$`5:36955532:T:G_G`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr5vars[chr5vars == "0"] <- NA
chr5vars <- unite(chr5vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr5vars)

# now to cbind with the ids df
chr5_final <- cbind(chr5ids, chr5vars)

# check to make sure things line up
chr5_final$VAR[10]
#  "5:141957473:C:G_G"
chr5raw$`5:141957473:C:G_G`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr5_final, "chr5_final.csv")
# now to save this to UKB project
system('dx upload chr5_final.csv --path Emily-folder/pathvar/chr5_final.csv')














# chr 6 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr6snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr6-snp-ids.txt')
head(chr6snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr6snps = chr6snps[-1]
head(chr6snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr6raw = fread('/mnt/project/plink-snp-data-20240108-path/chr6-snp-data/chr6-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr6raw)
# 469835   1984

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr6raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr6snps)
# output: 
colnames(chr6raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr6snp_ids <- strsplit(chr6snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr6snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr6 <- match(chr6snp_ids, colnames(chr6raw))
# 1857, so there are 127 that I need to remove



chr6colnames <- colnames(chr6raw)

# this still works - do this instead
extra_variants_chr6 <- setdiff(chr6colnames, chr6snp_ids)
indexchr6 <- match(extra_variants_chr6, colnames(chr6raw))
indexchr6 <- indexchr6[-c(1:6)]

# retain variants of interest by index
chr6raw <- chr6raw[, -indexchr6]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr6raw$Sums <- rowSums(chr6raw[, 7:1854], na.rm = TRUE)
# check
unique(chr6raw$Sums)
# 0 1 2 3

# remove the rows that have $Sum = 0
chr6raw <- chr6raw[chr6raw$Sums !=0, ]
# obs went from 469835 to 
unique(chr6raw[, 1855])
# 1 2 3


# make two dfs to merge together later
chr6vars <- chr6raw[, c(7:1854)]
chr6ids <- chr6raw[, c(1, 1855)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr6vars != 0, arr.ind = TRUE)
chr6vars[w] <- names(chr6vars)[w[,"col"]]

unique(chr6vars$`6:3224959:A:G_G`)
unique(chr6vars$`6:3225163:C:T_T`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr6vars[chr6vars == "0"] <- NA
chr6vars <- unite(chr6vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr6vars)

# now to cbind with the ids df
chr6_final <- cbind(chr6ids, chr6vars)

# check to make sure things line up
chr6_final$VAR[10]
# "6:33178948:C:T_T"
chr6raw$`6:33178948:C:T_T`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr6_final, "chr6_final.csv")
# now to save this to UKB project
system('dx upload chr6_final.csv --path Emily-folder/pathvar/chr6_final.csv')
















# chr 7 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr7snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr7-snp-ids.txt')
head(chr7snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr7snps = chr7snps[-1]
head(chr7snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr7raw = fread('/mnt/project/plink-snp-data-20240108-path/chr7-snp-data/chr7-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr7raw)
# 469835   2582

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr7raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr7snps)
# output: 
colnames(chr7raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr7snp_ids <- strsplit(chr7snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr7snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr7 <- match(chr7snp_ids, colnames(chr7raw))
# 2452, so there are 130 that I need to remove



chr7colnames <- colnames(chr7raw)

# this still works - do this instead
extra_variants_chr7 <- setdiff(chr7colnames, chr7snp_ids)
indexchr7 <- match(extra_variants_chr7, colnames(chr7raw))
indexchr7 <- indexchr7[-c(1:6)]

# retain variants of interest by index
chr7raw <- chr7raw[, -indexchr7]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr7raw$Sums <- rowSums(chr7raw[, 7:2367], na.rm = TRUE)
# check
unique(chr7raw$Sums)
# 0 1 2 3

# remove the rows that have $Sum = 0
chr7raw <- chr7raw[chr7raw$Sums !=0, ]
# obs went from 469835 to 16083
unique(chr7raw[,  2368])
# 1 2 3 


# make two dfs to merge together later
chr7vars <- chr7raw[, c(7:2367)]
chr7ids <- chr7raw[, c(1, 2368)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr7vars != 0, arr.ind = TRUE)
chr7vars[w] <- names(chr7vars)[w[,"col"]]

unique(chr7vars$`7:193335:T:A_A`)
unique(chr7vars$`7:193471:T:C_C`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr7vars[chr7vars == "0"] <- NA
chr7vars <- unite(chr7vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr7vars)

# now to cbind with the ids df
chr7_final <- cbind(chr7ids, chr7vars)

# check to make sure things line up
chr7_final$VAR[10]
#  "7:83134313:G:C_C"
chr7raw$`7:83134313:G:C_C`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr7_final, "chr7_final.csv")
# now to save this to UKB project
system('dx upload chr7_final.csv --path Emily-folder/pathvar/chr7_final.csv')














# chr 8 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr8snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr8-snp-ids.txt')
head(chr8snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr8snps = chr8snps[-1]
head(chr8snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr8raw = fread('/mnt/project/plink-snp-data-20240108-path/chr8-snp-data/chr8-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr8raw)
# 469835   2167

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr8raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr8snps)
# output: 
colnames(chr8raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr8snp_ids <- strsplit(chr8snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr8snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr8 <- match(chr8snp_ids, colnames(chr8raw))




chr8colnames <- colnames(chr8raw)

# this still works - do this instead
extra_variants_chr8 <- setdiff(chr8colnames, chr8snp_ids)
indexchr8 <- match(extra_variants_chr8, colnames(chr8raw))
indexchr8 <- indexchr8[-c(1:6)]

# retain variants of interest by index
chr8raw <- chr8raw[, -indexchr8]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr8raw$Sums <- rowSums(chr8raw[, 7:1994], na.rm = TRUE)
# check
unique(chr8raw$Sums)
# 0 1 2 3 4

# remove the rows that have $Sum = 0
chr8raw <- chr8raw[chr8raw$Sums !=0, ]
# obs went from 469835 to 16119
unique(chr8raw[, 1995])
# 1 2 3 4


# make two dfs to merge together later
chr8vars <- chr8raw[, c(7:1994)]
chr8ids <- chr8raw[, c(1, 1995)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr8vars != 0, arr.ind = TRUE)
chr8vars[w] <- names(chr8vars)[w[,"col"]]

unique(chr8vars$`8:6409278:G:A_A`)
unique(chr8vars$`8:6409324:A:G_G`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr8vars[chr8vars == "0"] <- NA
chr8vars <- unite(chr8vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr8vars)

# now to cbind with the ids df
chr8_final <- cbind(chr8ids, chr8vars)

# check to make sure things line up
chr8_final$VAR[10]
#  "8:22201829:G:A_A"
chr8raw$`8:22201829:G:A_A`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr8_final, "chr8_final.csv")
# now to save this to UKB project
system('dx upload chr8_final.csv --path Emily-folder/pathvar/chr8_final.csv')















# chr 9 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr9snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr9-snp-ids.txt')
head(chr9snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr9snps = chr9snps[-1]
head(chr9snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr9raw = fread('/mnt/project/plink-snp-data-20240108-path/chr9-snp-data/chr9-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr9raw)
# 469835   2231

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr9raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr9snps)
# output: 
colnames(chr9raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr9snp_ids <- strsplit(chr9snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr9snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr9 <- match(chr9snp_ids, colnames(chr9raw))
# 2065, so there are 166 that I need to remove



chr9colnames <- colnames(chr9raw)

# this still works - do this instead
extra_variants_chr9 <- setdiff(chr9colnames, chr9snp_ids)
indexchr9 <- match(extra_variants_chr9, colnames(chr9raw))
indexchr9 <- indexchr9[-c(1:6)]

# retain variants of interest by index
chr9raw <- chr9raw[, -indexchr9]
# 2047






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr9raw$Sums <- rowSums(chr9raw[, 7:2047], na.rm = TRUE)
# check
unique(chr9raw$Sums)
# 0 1 2 3

# remove the rows that have $Sum = 0
chr9raw <- chr9raw[chr9raw$Sums !=0, ]
# obs went from 469835 to 
unique(chr9raw[,  2048])
# 1 2 3


# make two dfs to merge together later
chr9vars <- chr9raw[, c(7:2047)]
chr9ids <- chr9raw[, c(1, 2048)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr9vars != 0, arr.ind = TRUE)
chr9vars[w] <- names(chr9vars)[w[,"col"]]

unique(chr9vars$`9:35792485:C:T_T`)
unique(chr9vars$`9:35792964:G:A_A`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr9vars[chr9vars == "0"] <- NA
chr9vars <- unite(chr9vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr9vars)

# now to cbind with the ids df
chr9_final <- cbind(chr9ids, chr9vars)

# check to make sure things line up
chr9_final$VAR[10]
# "9:114306644:G:A_A"
chr9raw$`9:114306644:G:A_A`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr9_final, "chr9_final.csv")
# now to save this to UKB project
system('dx upload chr9_final.csv --path Emily-folder/pathvar/chr9_final.csv')














# chr 10 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr10snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr10-snp-ids.txt')
head(chr10snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr10snps = chr10snps[-1]
head(chr10snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr10raw = fread('/mnt/project/plink-snp-data-20240108-path/chr10-snp-data/chr10-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr10raw)
# 469835   2198

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr10raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr10snps)
# output: 
colnames(chr10raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr10snp_ids <- strsplit(chr10snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr10snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr10 <- match(chr10snp_ids, colnames(chr10raw))
# 2063, so there are 135 that I need to remove



chr10colnames <- colnames(chr10raw)

# this still works - do this instead
extra_variants_chr10 <- setdiff(chr10colnames, chr10snp_ids)
indexchr10 <- match(extra_variants_chr10, colnames(chr10raw))
indexchr10 <- indexchr10[-c(1:6)]

# retain variants of interest by index
chr10raw <- chr10raw[, -indexchr10]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr10raw$Sums <- rowSums(chr10raw[, 7:2056], na.rm = TRUE)
# check
unique(chr10raw$Sums)
# 0 1 2 3

# remove the rows that have $Sum = 0
chr10raw <- chr10raw[chr10raw$Sums !=0, ]
# obs went from 469835 to 24039
unique(chr10raw[,  2057])
# 1 2 3


# make two dfs to merge together later
chr10vars <- chr10raw[, c(7:2056)]
chr10ids <- chr10raw[, c(1, 2057)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr10vars != 0, arr.ind = TRUE)
chr10vars[w] <- names(chr10vars)[w[,"col"]]

unique(chr10vars$`10:13651913:C:T_T`)
unique(chr10vars$`10:13654485:G:C_C`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr10vars[chr10vars == "0"] <- NA
chr10vars <- unite(chr10vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr10vars)

# now to cbind with the ids df
chr10_final <- cbind(chr10ids, chr10vars)

# check to make sure things line up
chr10_final$VAR[10]
#  "10:49470310:A:C_C"
chr10raw$`10:49470310:A:C_C`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr10_final, "chr10_final.csv")
# now to save this to UKB project
system('dx upload chr10_final.csv --path Emily-folder/pathvar/chr10_final.csv')
















# chr 11 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr11snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr11-snp-ids.txt')
head(chr11snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr11snps = chr11snps[-1]
head(chr11snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr11raw = fread('/mnt/project/plink-snp-data-20240108-path/chr11-snp-data/chr11-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr11raw)
# 469835   3760

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr11raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr11snps)
# output: 
colnames(chr11raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr11snp_ids <- strsplit(chr11snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr11snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr11 <- match(chr11snp_ids, colnames(chr11raw))




chr11colnames <- colnames(chr11raw)

# this still works - do this instead
extra_variants_chr11 <- setdiff(chr11colnames, chr11snp_ids)
indexchr11 <- match(extra_variants_chr11, colnames(chr11raw))
indexchr11 <- indexchr11[-c(1:6)]

# retain variants of interest by index
chr11raw <- chr11raw[, -indexchr11]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr11raw$Sums <- rowSums(chr11raw[, 7:3479], na.rm = TRUE)
# check
unique(chr11raw$Sums)
# 0 1 2 3 4

# remove the rows that have $Sum = 0
chr11raw <- chr11raw[chr11raw$Sums !=0, ]
# obs went from 469835 to 
unique(chr11raw[,  3480])
# 1 2 3 4


# make two dfs to merge together later
chr11vars <- chr11raw[, c(7:3479)]
chr11ids <- chr11raw[, c(1, 3480)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr11vars != 0, arr.ind = TRUE)
chr11vars[w] <- names(chr11vars)[w[,"col"]]

unique(chr11vars$`11:747503:C:T_T`)
unique(chr11vars$`11:759033:G:A_A`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr11vars[chr11vars == "0"] <- NA
chr11vars <- unite(chr11vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr11vars)

# now to cbind with the ids df
chr11_final <- cbind(chr11ids, chr11vars)

# check to make sure things line up
chr11_final$VAR[10]
#  "11:71435840:C:G_G"
chr11raw$`11:71435840:C:G_G`[10]
# 1 

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr11_final, "chr11_final.csv")
# now to save this to UKB project
system('dx upload chr11_final.csv --path Emily-folder/pathvar/chr11_final.csv')
















# chr 12 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr12snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr12-snp-ids.txt')
head(chr12snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr12snps = chr12snps[-1]
head(chr12snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr12raw = fread('/mnt/project/plink-snp-data-20240108-path/chr12-snp-data/chr12-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr12raw)
# 469835   4060

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr12raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr12snps)
# output: 
colnames(chr12raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr12snp_ids <- strsplit(chr12snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr12snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr12 <- match(chr12snp_ids, colnames(chr12raw))




chr12colnames <- colnames(chr12raw)

# this still works - do this instead
extra_variants_chr12 <- setdiff(chr12colnames, chr12snp_ids)
indexchr12 <- match(extra_variants_chr12, colnames(chr12raw))
indexchr12 <- indexchr12[-c(1:6)]

# retain variants of interest by index
chr12raw <- chr12raw[, -indexchr12]





### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr12raw$Sums <- rowSums(chr12raw[, 7:3761], na.rm = TRUE)
# check
unique(chr12raw$Sums)
# 0 1 3 2

# remove the rows that have $Sum = 0
chr12raw <- chr12raw[chr12raw$Sums !=0, ]
# obs went from 469835 to 3762
unique(chr12raw[, 3762])
# 1 3 2


# make two dfs to merge together later
chr12vars <- chr12raw[, c(7:3761)]
chr12ids <- chr12raw[, c(1, 3762)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr12vars != 0, arr.ind = TRUE)
chr12vars[w] <- names(chr12vars)[w[,"col"]]

unique(chr12vars$`12:6495122:C:G_G`)
unique(chr12vars$`12:6510719:A:G_G`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr12vars[chr12vars == "0"] <- NA
chr12vars <- unite(chr12vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr12vars)

# now to cbind with the ids df
chr12_final <- cbind(chr12ids, chr12vars)

# check to make sure things line up
chr12_final$VAR[10]
#  "12:101768082:C:A_A"
chr12raw$`12:101768082:C:A_A`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr12_final, "chr12_final.csv")
# now to save this to UKB project
system('dx upload chr12_final.csv --path Emily-folder/pathvar/chr12_final.csv')















# chr 13 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr13snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr13-snp-ids.txt')
head(chr13snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr13snps = chr13snps[-1]
head(chr13snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr13raw = fread('/mnt/project/plink-snp-data-20240108-path/chr13-snp-data/chr13-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr13raw)
# 469835    752

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr13raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr13snps)
# output: 
colnames(chr13raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr13snp_ids <- strsplit(chr13snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr13snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr13 <- match(chr13snp_ids, colnames(chr13raw))




chr13colnames <- colnames(chr13raw)

# this still works - do this instead
extra_variants_chr13 <- setdiff(chr13colnames, chr13snp_ids)
indexchr13 <- match(extra_variants_chr13, colnames(chr13raw))
indexchr13 <- indexchr13[-c(1:6)]

# retain variants of interest by index
chr13raw <- chr13raw[, -indexchr13]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr13raw$Sums <- rowSums(chr13raw[, 7:691], na.rm = TRUE)
# check
unique(chr13raw$Sums)
# 0 1 2 

# remove the rows that have $Sum = 0
chr13raw <- chr13raw[chr13raw$Sums !=0, ]
# obs went from 469835 to 7062
unique(chr13raw[, 692])
# 1 2


# make two dfs to merge together later
chr13vars <- chr13raw[, c(7:691)]
chr13ids <- chr13raw[, c(1, 692)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr13vars != 0, arr.ind = TRUE)
chr13vars[w] <- names(chr13vars)[w[,"col"]]

unique(chr13vars$`13:24883193:A:G_G`)
unique(chr13vars$`13:24883306:C:G_G`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr13vars[chr13vars == "0"] <- NA
chr13vars <- unite(chr13vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr13vars)

# now to cbind with the ids df
chr13_final <- cbind(chr13ids, chr13vars)

# check to make sure things line up
chr13_final$VAR[10]
# "13:31269278:G:A_A"
chr13raw$`13:31269278:G:A_A`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr13_final, "chr13_final.csv")
# now to save this to UKB project
system('dx upload chr13_final.csv --path Emily-folder/pathvar/chr13_final.csv')
















# chr 14 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr14snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr14-snp-ids.txt')
head(chr14snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr14snps = chr14snps[-1]
head(chr14snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr14raw = fread('/mnt/project/plink-snp-data-20240108-path/chr14-snp-data/chr14-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr14raw)
# 469835   1259

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr14raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr14snps)
# output: 
colnames(chr14raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr14snp_ids <- strsplit(chr14snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr14snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr14 <- match(chr14snp_ids, colnames(chr14raw))




chr14colnames <- colnames(chr14raw)

# this still works - do this instead
extra_variants_chr14 <- setdiff(chr14colnames, chr14snp_ids)
indexchr14 <- match(extra_variants_chr14, colnames(chr14raw))
indexchr14 <- indexchr14[-c(1:6)]

# retain variants of interest by index
chr14raw <- chr14raw[, -indexchr14]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr14raw$Sums <- rowSums(chr14raw[, 7:1154], na.rm = TRUE)
# check
unique(chr14raw$Sums)
# 0 1 2 

# remove the rows that have $Sum = 0
chr14raw <- chr14raw[chr14raw$Sums !=0, ]
# obs went from 469835 to 
unique(chr14raw[,  1155])
# 1 2 


# make two dfs to merge together later
chr14vars <- chr14raw[, c(7:1154)]
chr14ids <- chr14raw[, c(1, 1155)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr14vars != 0, arr.ind = TRUE)
chr14vars[w] <- names(chr14vars)[w[,"col"]]

unique(chr14vars$`14:20447274:C:G_G`)
unique(chr14vars$`14:20447646:G:A_A`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr14vars[chr14vars == "0"] <- NA
chr14vars <- unite(chr14vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr14vars)

# now to cbind with the ids df
chr14_final <- cbind(chr14ids, chr14vars)

# check to make sure things line up
chr14_final$VAR[10]
#  "14:56804334:G:T_T"
chr14raw$`14:56804334:G:T_T`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr14_final, "chr14_final.csv")
# now to save this to UKB project
system('dx upload chr14_final.csv --path Emily-folder/pathvar/chr14_final.csv')

















# chr 15 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr15snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr15-snp-ids.txt')
head(chr15snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr15snps = chr15snps[-1]
head(chr15snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr15raw = fread('/mnt/project/plink-snp-data-20240108-path/chr15-snp-data/chr15-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr15raw)
# 469835   2149

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr15raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr15snps)
# output: 
colnames(chr15raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr15snp_ids <- strsplit(chr15snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr15snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr15 <- match(chr15snp_ids, colnames(chr15raw))




chr15colnames <- colnames(chr15raw)

# this still works - do this instead
extra_variants_chr15 <- setdiff(chr15colnames, chr15snp_ids)
indexchr15 <- match(extra_variants_chr15, colnames(chr15raw))
indexchr15 <- indexchr15[-c(1:6)]

# retain variants of interest by index
chr15raw <- chr15raw[, -indexchr15]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr15raw$Sums <- rowSums(chr15raw[, 7:1993], na.rm = TRUE)
# check
unique(chr15raw$Sums)
# 0 1 2 3

# remove the rows that have $Sum = 0
chr15raw <- chr15raw[chr15raw$Sums !=0, ]
# obs went from 469835 to 
unique(chr15raw[,  1994])
# 1 2 3 


# make two dfs to merge together later
chr15vars <- chr15raw[, c(7:1993)]
chr15ids <- chr15raw[, c(1, 1994)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr15vars != 0, arr.ind = TRUE)
chr15vars[w] <- names(chr15vars)[w[,"col"]]

unique(chr15vars$`15:40165078:G:A_A`)
unique(chr15vars$`15:40170117:G:T_T`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr15vars[chr15vars == "0"] <- NA
chr15vars <- unite(chr15vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr15vars)

# now to cbind with the ids df
chr15_final <- cbind(chr15ids, chr15vars)

# check to make sure things line up
chr15_final$VAR[10]
#  "15:79889204:T:C_C"
chr15raw$`15:79889204:T:C_C`[10]
# 1 

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr15_final, "chr15_final.csv")
# now to save this to UKB project
system('dx upload chr15_final.csv --path Emily-folder/pathvar/chr15_final.csv')

















# chr 16 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr16snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr16-snp-ids.txt')
head(chr16snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr16snps = chr16snps[-1]
head(chr16snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr16raw = fread('/mnt/project/plink-snp-data-20240108-path/chr16-snp-data/chr16-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr16raw)
# 469835   3280

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr16raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr16snps)
# output: 
colnames(chr16raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr16snp_ids <- strsplit(chr16snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr16snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr16 <- match(chr16snp_ids, colnames(chr16raw))




chr16colnames <- colnames(chr16raw)

# this still works - do this instead
extra_variants_chr16 <- setdiff(chr16colnames, chr16snp_ids)
indexchr16 <- match(extra_variants_chr16, colnames(chr16raw))
indexchr16 <- indexchr16[-c(1:6)]

# retain variants of interest by index
chr16raw <- chr16raw[, -indexchr16]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr16raw$Sums <- rowSums(chr16raw[, 7:2993], na.rm = TRUE)
# check
unique(chr16raw$Sums)
# 0 1 2 3

# remove the rows that have $Sum = 0
chr16raw <- chr16raw[chr16raw$Sums !=0, ]
# obs went from 469835 to 
unique(chr16raw[,  2994])
# 1 2 3 


# make two dfs to merge together later
chr16vars <- chr16raw[, c(7:2993)]
chr16ids <- chr16raw[, c(1, 2994)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr16vars != 0, arr.ind = TRUE)
chr16vars[w] <- names(chr16vars)[w[,"col"]]

unique(chr16vars$`16:1790634:C:T_T`)
unique(chr16vars$`16:1791032:G:C_C`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr16vars[chr16vars == "0"] <- NA
chr16vars <- unite(chr16vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr16vars)

# now to cbind with the ids df
chr16_final <- cbind(chr16ids, chr16vars)

# check to make sure things line up
chr16_final$VAR[10]
#  "16:88835251:G:A_A"
chr16raw$`16:88835251:G:A_A`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr16_final, "chr16_final.csv")
# now to save this to UKB project
system('dx upload chr16_final.csv --path Emily-folder/pathvar/chr16_final.csv')
















# chr 17 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr17snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr17-snp-ids.txt')
head(chr17snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr17snps = chr17snps[-1]
head(chr17snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr17raw = fread('/mnt/project/plink-snp-data-20240108-path/chr17-snp-data/chr17-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr17raw)
# 469835   2566

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr17raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr17snps)
# output: 
colnames(chr17raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr17snp_ids <- strsplit(chr17snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr17snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr17 <- match(chr17snp_ids, colnames(chr17raw))




chr17colnames <- colnames(chr17raw)

# this still works - do this instead
extra_variants_chr17 <- setdiff(chr17colnames, chr17snp_ids)
indexchr17 <- match(extra_variants_chr17, colnames(chr17raw))
indexchr17 <- indexchr17[-c(1:6)]

# retain variants of interest by index
chr17raw <- chr17raw[, -indexchr17]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr17raw$Sums <- rowSums(chr17raw[, 7:2399], na.rm = TRUE)
# check
unique(chr17raw$Sums)
# 0 1 2 3

# remove the rows that have $Sum = 0
chr17raw <- chr17raw[chr17raw$Sums !=0, ]
# obs went from 469835 to 
unique(chr17raw[, 2400])
# 1 2 3


# make two dfs to merge together later
chr17vars <- chr17raw[, c(7:2399)]
chr17ids <- chr17raw[, c(1, 2400)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr17vars != 0, arr.ind = TRUE)
chr17vars[w] <- names(chr17vars)[w[,"col"]]

unique(chr17vars$`17:800972:T:G_G`)
unique(chr17vars$`17:803758:G:A_A`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr17vars[chr17vars == "0"] <- NA
chr17vars <- unite(chr17vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr17vars)

# now to cbind with the ids df
chr17_final <- cbind(chr17ids, chr17vars)

# check to make sure things line up
chr17_final$VAR[10]
#  "17:2039805:C:G_G"
chr17raw$`17:2039805:C:G_G`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr17_final, "chr17_final.csv")
# now to save this to UKB project
system('dx upload chr17_final.csv --path Emily-folder/pathvar/chr17_final.csv')

















# chr 18 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr18snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr18-snp-ids.txt')
head(chr18snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr18snps = chr18snps[-1]
head(chr18snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr18raw = fread('/mnt/project/plink-snp-data-20240108-path/chr18-snp-data/chr18-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr18raw)
# 469835    917

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr18raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr18snps)
# output: 
colnames(chr18raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr18snp_ids <- strsplit(chr18snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr18snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr18 <- match(chr18snp_ids, colnames(chr18raw))




chr18colnames <- colnames(chr18raw)

# this still works - do this instead
extra_variants_chr18 <- setdiff(chr18colnames, chr18snp_ids)
indexchr18 <- match(extra_variants_chr18, colnames(chr18raw))
indexchr18 <- indexchr18[-c(1:6)]

# retain variants of interest by index
chr18raw <- chr18raw[, -indexchr18]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr18raw$Sums <- rowSums(chr18raw[, 7:859], na.rm = TRUE)
# check
unique(chr18raw$Sums)
# 0 1 2 

# remove the rows that have $Sum = 0
chr18raw <- chr18raw[chr18raw$Sums !=0, ]
# obs went from 469835 to 
unique(chr18raw[, 860])
# 1 2


# make two dfs to merge together later
chr18vars <- chr18raw[, c(7:859)]
chr18ids <- chr18raw[, c(1, 860)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr18vars != 0, arr.ind = TRUE)
chr18vars[w] <- names(chr18vars)[w[,"col"]]

unique(chr18vars$`18:22936874:G:T_T`)
unique(chr18vars$`18:22949697:A:G_G`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr18vars[chr18vars == "0"] <- NA
chr18vars <- unite(chr18vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr18vars)

# now to cbind with the ids df
chr18_final <- cbind(chr18ids, chr18vars)

# check to make sure things line up
chr18_final$VAR[10]
#  "18:23022237:G:C_C"
chr18raw$`18:23022237:G:C_C`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr18_final, "chr18_final.csv")
# now to save this to UKB project
system('dx upload chr18_final.csv --path Emily-folder/pathvar/chr18_final.csv')



















# chr 19 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr19snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr19-snp-ids.txt')
head(chr19snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr19snps = chr19snps[-1]
head(chr19snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr19raw = fread('/mnt/project/plink-snp-data-20240108-path/chr19-snp-data/chr19-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr19raw)
# 469835   2190

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr19raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr19snps)
# output: 
colnames(chr19raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr19snp_ids <- strsplit(chr19snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr19snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr19 <- match(chr19snp_ids, colnames(chr19raw))




chr19colnames <- colnames(chr19raw)

# this still works - do this instead
extra_variants_chr19 <- setdiff(chr19colnames, chr19snp_ids)
indexchr19 <- match(extra_variants_chr19, colnames(chr19raw))
indexchr19 <- indexchr19[-c(1:6)]

# retain variants of interest by index
chr19raw <- chr19raw[, -indexchr19]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr19raw$Sums <- rowSums(chr19raw[, 7:1989], na.rm = TRUE)
# check
unique(chr19raw$Sums)
# 0 1 2 3

# remove the rows that have $Sum = 0
chr19raw <- chr19raw[chr19raw$Sums !=0, ]
# obs went from 469835 to 
unique(chr19raw[, 1990])
# 1 2 3  


# make two dfs to merge together later
chr19vars <- chr19raw[, c(7:1989)]
chr19ids <- chr19raw[, c(1, 1990)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr19vars != 0, arr.ind = TRUE)
chr19vars[w] <- names(chr19vars)[w[,"col"]]

unique(chr19vars$`19:1104045:T:A_A`)
unique(chr19vars$`19:1105461:A:T_T`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr19vars[chr19vars == "0"] <- NA
chr19vars <- unite(chr19vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr19vars)

# now to cbind with the ids df
chr19_final <- cbind(chr19ids, chr19vars)

# check to make sure things line up
chr19_final$VAR[10]
#  "19:36067377:C:G_G"
chr19raw$`19:36067377:C:G_G`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr19_final, "chr19_final.csv")
# now to save this to UKB project
system('dx upload chr19_final.csv --path Emily-folder/pathvar/chr19_final.csv')


















# chr 20 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr20snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr20-snp-ids.txt')
head(chr20snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr20snps = chr20snps[-1]
head(chr20snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr20raw = fread('/mnt/project/plink-snp-data-20240108-path/chr20-snp-data/chr20-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr20raw)
# 469835   1205

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr20raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr20snps)
# output: 
colnames(chr20raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr20snp_ids <- strsplit(chr20snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr20snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr20 <- match(chr20snp_ids, colnames(chr20raw))




chr20colnames <- colnames(chr20raw)

# this still works - do this instead
extra_variants_chr20 <- setdiff(chr20colnames, chr20snp_ids)
indexchr20 <- match(extra_variants_chr20, colnames(chr20raw))
indexchr20 <- indexchr20[-c(1:6)]

# retain variants of interest by index
chr20raw <- chr20raw[, -indexchr20]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr20raw$Sums <- rowSums(chr20raw[, 7:1109], na.rm = TRUE)
# check
unique(chr20raw$Sums)
# 0 1 2

# remove the rows that have $Sum = 0
chr20raw <- chr20raw[chr20raw$Sums !=0, ]
# obs went from  to 
unique(chr20raw[, 1110])
# 1 2 


# make two dfs to merge together later
chr20vars <- chr20raw[, c(7: 1109)]
chr20ids <- chr20raw[, c(1, 1110)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr20vars != 0, arr.ind = TRUE)
chr20vars[w] <- names(chr20vars)[w[,"col"]]

unique(chr20vars$`20:3190667:G:A_A`)
unique(chr20vars$`20:3190808:C:T_T`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr20vars[chr20vars == "0"] <- NA
chr20vars <- unite(chr20vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr20vars)

# now to cbind with the ids df
chr20_final <- cbind(chr20ids, chr20vars)

# check to make sure things line up
chr20_final$VAR[10]
#  "20:58853803:C:T_T"
chr20raw$`20:58853803:C:T_T`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr20_final, "chr20_final.csv")
# now to save this to UKB project
system('dx upload chr20_final.csv --path Emily-folder/pathvar/chr20_final.csv')
















# chr 21 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr21snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr21-snp-ids.txt')
head(chr21snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr21snps = chr21snps[-1]
head(chr21snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr21raw = fread('/mnt/project/plink-snp-data-20240108-path/chr21-snp-data/chr21-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr21raw)
# 469835    816

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr21raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr21snps)
# output: 
colnames(chr21raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr21snp_ids <- strsplit(chr21snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr21snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr21 <- match(chr21snp_ids, colnames(chr21raw))




chr21colnames <- colnames(chr21raw)

# this still works - do this instead
extra_variants_chr21 <- setdiff(chr21colnames, chr21snp_ids)
indexchr21 <- match(extra_variants_chr21, colnames(chr21raw))
indexchr21 <- indexchr21[-c(1:6)]

# retain variants of interest by index
chr21raw <- chr21raw[, -indexchr21]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr21raw$Sums <- rowSums(chr21raw[, 7:730], na.rm = TRUE)
# check
unique(chr21raw$Sums)
# 0 1 2 

# remove the rows that have $Sum = 0
chr21raw <- chr21raw[chr21raw$Sums !=0, ]
# obs went from  to 
unique(chr21raw[, 731])
# 1 2 


# make two dfs to merge together later
chr21vars <- chr21raw[, c(7:730)]
chr21ids <- chr21raw[, c(1, 731)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr21vars != 0, arr.ind = TRUE)
chr21vars[w] <- names(chr21vars)[w[,"col"]]

unique(chr21vars$`21:33578332:T:G_G`)
unique(chr21vars$`21:33579480:G:A_A`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr21vars[chr21vars == "0"] <- NA
chr21vars <- unite(chr21vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr21vars)

# now to cbind with the ids df
chr21_final <- cbind(chr21ids, chr21vars)

# check to make sure things line up
chr21_final$VAR[10]
#  "21:46398042:G:A_A"
chr21raw$`21:46398042:G:A_A`[10]
# 1 

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr21_final, "chr21_final.csv")
# now to save this to UKB project
system('dx upload chr21_final.csv --path Emily-folder/pathvar/chr21_final.csv')


















# chr 22 [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chr22snps = readLines('/mnt/project/plink-snp-data-20240108-path/chr22-snp-ids.txt')
head(chr22snps)
# output as expected: 

# Get rid of first entry ("chr:pos:ref:alt")
chr22snps = chr22snps[-1]
head(chr22snps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chr22raw = fread('/mnt/project/plink-snp-data-20240108-path/chr22-snp-data/chr22-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chr22raw)
# 469835   2781

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chr22raw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chr22snps)
# output: 
colnames(chr22raw)[7:10]
# output: 

# Create new variable with extra "_X" added:
chr22snp_ids <- strsplit(chr22snps , ":") %>% 
  lapply(., function(x) paste0(x[1], ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chr22snp_ids)
# output: 

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chr22 <- match(chr22snp_ids, colnames(chr22raw))




chr22colnames <- colnames(chr22raw)

# this still works - do this instead
extra_variants_chr22 <- setdiff(chr22colnames, chr22snp_ids)
indexchr22 <- match(extra_variants_chr22, colnames(chr22raw))
indexchr22 <- indexchr22[-c(1:6)]

# retain variants of interest by index
chr22raw <- chr22raw[, -indexchr22]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chr22raw$Sums <- rowSums(chr22raw[, 7:2554], na.rm = TRUE)
# check
unique(chr22raw$Sums)
# 0 1 2 3 

# remove the rows that have $Sum = 0
chr22raw <- chr22raw[chr22raw$Sums !=0, ]
# obs went from  to 
unique(chr22raw[,  2555])
# 1 2 3 


# make two dfs to merge together later
chr22vars <- chr22raw[, c(7:2554)]
chr22ids <- chr22raw[, c(1, 2555)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chr22vars != 0, arr.ind = TRUE)
chr22vars[w] <- names(chr22vars)[w[,"col"]]

unique(chr22vars$`22:19479971:G:A_A`)
unique(chr22vars$`22:19480984:C:T_T`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chr22vars[chr22vars == "0"] <- NA
chr22vars <- unite(chr22vars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chr22vars)

# now to cbind with the ids df
chr22_final <- cbind(chr22ids, chr22vars)

# check to make sure things line up
chr22_final$VAR[10]
#  "22:31648181:G:A_A"
chr22raw$`22:31648181:G:A_A`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chr22_final, "chr22_final.csv")
# now to save this to UKB project
system('dx upload chr22_final.csv --path Emily-folder/pathvar/chr22_final.csv')
















# chr X [DONE]-------------------------------------------------------------------


### load in the data


# Load PLINK format SNP names for snps
chrXsnps = readLines('/mnt/project/plink-snp-data-20240108-path/chrX-snp-ids.txt')
head(chrXsnps)
# output as expected: "chr:pos:ref:alt" "X:110198535:G:C" "X:110198585:G:A" "X:110198619:C:G" "X:110198629:C:T" "X:110198630:G:C"

# Get rid of first entry ("chr:pos:ref:alt")
chrXsnps = chrXsnps[-1]
head(chrXsnps)
# output as expected

# Use data.table (faster) to load the raw data from PLINK
chrXraw = fread('/mnt/project/plink-snp-data-20240108-path/chrX-snp-data/chrX-snps.raw', data.table=FALSE)
# hopefully this one doesn't have any frustrating rows!!!
dim(chrXraw)
# 469835    629

# Take a look at the first 5 rows, and first 8 columns
# First 6 columns are sample info
chrXraw[1:5,1:8]
# looks correct



### remove unwanted variants



# We only want the SNPs that we supplied
# First need to get the SNP id formats matching
head(chrXsnps)
# output: "X:110198535:G:C" "X:110198585:G:A" "X:110198619:C:G" "X:110198629:C:T" "X:110198630:G:C" "X:110200957:G:A"
colnames(chrXraw)[7:10]
# output: "23:630961:G:T_T" "23:630961:G:A_A" "23:631004:G:A_A" "23:631007:A:G_G"

# Create new variable with extra "_X" added: (also change chr from X to 23)
chrXsnp_ids <- strsplit(chrXsnps , ":") %>% 
  lapply(., function(x) paste0("23", ":", x[2], ":", x[3], ":", x[4], "_", x[4])) %>% 
  unlist()
# Now they are the same format...
head(chrXsnp_ids)
# output: "23:110198535:G:C_C" "23:110198585:G:A_A" "23:110198619:C:G_G" "23:110198629:C:T_T" "23:110198630:G:C_C"

# Which columns of the raw data match Emily's SNPs of interest?
match_ids_chrX <- match(chrXsnp_ids, colnames(chrXraw))




chrXcolnames <- colnames(chrXraw)

# this still works - do this instead
extra_variants_chrX <- setdiff(chrXcolnames, chrXsnp_ids)
indexchrX <- match(extra_variants_chrX, colnames(chrXraw))
indexchrX <- indexchrX[-c(1:6)]

# retain variants of interest by index
chrXraw <- chrXraw[, -indexchrX]






### reformat the df into something I can save and use


# add a new column that contains the sum of SNP columns for each row(participant)
chrXraw$Sums <- rowSums(chrXraw[, 7:588], na.rm = TRUE)
# check
unique(chrXraw$Sums)
# 0 1 2 4 3

# remove the rows that have $Sum = 0
chrXraw <- chrXraw[chrXraw$Sums !=0, ]
# obs went from  to 
unique(chrXraw[, 589])
# 1 2 4 3


# make two dfs to merge together later
chrXvars <- chrXraw[, c(7:588)]
chrXids <- chrXraw[, c(1, 589)]


# now to save variant information into a variable
# code obtained from here: https://stackoverflow.com/questions/52650425/how-do-i-replace-1s-in-data-frame-with-the-column-name
#install.packages('tidyverse')
#library(tidyverse)

# EMILY: check that the index numbers are correct for the size of the df!!!!!

# first: print the column name in each cell that doesn’t contain ‘0’
# m1 <- col(chr1raw[2:6530]) * chr1raw[2:6530]
# i1 <- m1 != 0
# chr1raw[2:6530][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
##### currently getting stuck here:
# Error in `[<-.data.frame`(`*tmp*`, i1, value = c("1:1232291:C:T_T", "1:1232291:C:T_T",  : 
# 'value' is the wrong length


# trying something different
w <- which(chrXvars != 0, arr.ind = TRUE)
chrXvars[w] <- names(chrXvars)[w[,"col"]]

unique(chrXvars$`23:631004:G:A_A`)
unique(chrXvars$`23:644512:C:G_G`)


library(tidyr)
# now merge all of the columns together, ignoring NA, and separating by ‘;’
# first make 0 = NA
chrXvars[chrXvars == "0"] <- NA
chrXvars <- unite(chrXvars, "VAR", sep = ";", na.rm = TRUE)

#check
unique(chrXvars)

# now to cbind with the ids df
chrX_final <- cbind(chrXids, chrXvars)

# check to make sure things line up
chrX_final$VAR[10]
#  "23:22096977:G:A_A"
chrXraw$`23:22096977:G:A_A`[10]
# 1

# save this df into the pathvar file - to combine with the others later on
# combining all together at the same time would take too much memory

#save it to R
write.csv(chrX_final, "chrX_final.csv")
# now to save this to UKB project
system('dx upload chrX_final.csv --path Emily-folder/pathvar/chrX_final.csv')













# saving participant IDs --------------------------------------------------

library(data.table)
library(dplyr)


# Use data.table (faster) to load the raw data from PLINK
chr21raw = fread('/mnt/project/plink-snp-data-20240108-path/chr21-snp-data/chr21-snps.raw', data.table=FALSE)

# Use data.table (faster) to load the raw data from PLINK
chr22raw = fread('/mnt/project/plink-snp-data-20240108-path/chr22-snp-data/chr22-snps.raw', data.table=FALSE)


setdiff(chr22raw$FID, chr21raw$FID)
setdiff(chr21raw$FID, chr22raw$FID)
# 0
# so this is probably the same for all 


#checking that this is the same for the MGSvar data:
donsonMGSraw = fread('/mnt/project/plink-snp-data-20230405/donson-MGSvar-snp-data/donson-MGSvar-snps.raw', data.table=FALSE)
setdiff(chr22raw$FID, donsonMGSraw$FID)
setdiff(donsonMGSraw$FID, chr22raw$FID)
# 0

# can make a list of genotyped participants from this!
genotyped_participants <- c(chr21raw$FID)


#save it to R
write.csv(genotyped_participants, "genotyped_participants.csv")
# now to save this to UKB project
system('dx upload genotyped_participants.csv --path Emily-folder/genotyped_participants.csv')








system('dx upload pathvar_script_END.R --path Emily-folder/pathvar/pathvar_script_END.R')