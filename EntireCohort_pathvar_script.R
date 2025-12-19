# EntireCohort_pathvar_script.R

# this script aims to annotate pathvar data - creates path_var_data.csv 
# this contains three parts (combine CHRs, annotate by GENE, and annotate by phenotypes and cancer association)

# the second part of this script uses path_var_data.csv to annotate the EntireCohort_ready.csv dataframe
# this creates EntireCohort_pathvar.csv


# combine CHRs ------------------------------------------------------------


# this first section was completed on 22/02/24 and does not need to be executed again:
# using the chr files for the pathvar data to annotate the EntireCohort_Cancer_2 df
library(dplyr)
library(readr)

chr1_final <- read_csv("pathvar/chr1_final.csv")
chr2_final <- read_csv("pathvar/chr2_final.csv")
chr3_final <- read_csv("pathvar/chr3_final.csv")
chr4_final <- read_csv("pathvar/chr4_final.csv")
chr5_final <- read_csv("pathvar/chr5_final.csv")
chr6_final <- read_csv("pathvar/chr6_final.csv")
chr7_final <- read_csv("pathvar/chr7_final.csv")
chr8_final <- read_csv("pathvar/chr8_final.csv")
chr9_final <- read_csv("pathvar/chr9_final.csv")
chr10_final <- read_csv("pathvar/chr10_final.csv")
chr11_final <- read_csv("pathvar/chr11_final.csv")
chr12_final <- read_csv("pathvar/chr12_final.csv")
chr13_final <- read_csv("pathvar/chr13_final.csv")
chr14_final <- read_csv("pathvar/chr14_final.csv")
chr15_final <- read_csv("pathvar/chr15_final.csv")
chr16_final <- read_csv("pathvar/chr16_final.csv")
chr17_final <- read_csv("pathvar/chr17_final.csv")
chr18_final <- read_csv("pathvar/chr18_final.csv")
chr19_final <- read_csv("pathvar/chr19_final.csv")
chr20_final <- read_csv("pathvar/chr20_final.csv")
chr21_final <- read_csv("pathvar/chr21_final.csv")
chr22_final <- read_csv("pathvar/chr22_final.csv")
chrX_final <- read_csv("pathvar/chrX_final.csv")

# remove the first column from each
chr1_final <- chr1_final[-1]
chr2_final <- chr2_final[-1]
chr3_final <- chr3_final[-1]
chr4_final <- chr4_final[-1]
chr5_final <- chr5_final[-1]
chr6_final <- chr6_final[-1]
chr7_final <- chr7_final[-1]
chr8_final <- chr8_final[-1]
chr9_final <- chr9_final[-1]
chr10_final <- chr10_final[-1]
chr11_final <- chr11_final[-1]
chr12_final <- chr12_final[-1]
chr13_final <- chr13_final[-1]
chr14_final <- chr14_final[-1]
chr15_final <- chr15_final[-1]
chr16_final <- chr16_final[-1]
chr17_final <- chr17_final[-1]
chr18_final <- chr18_final[-1]
chr19_final <- chr19_final[-1]
chr20_final <- chr20_final[-1]
chr21_final <- chr21_final[-1]
chr22_final <- chr22_final[-1]
chrX_final <- chrX_final[-1]




# want to combine all of these into one
# need to be careful with adding using FID as reference
# sum all of the variant Sums
# combine all of the VAR columns

## to combine into one dataframe:
# use full_join with the argument 'by=FID'
df1 <- full_join(chr1_final, chr2_final, by='FID')
df2 <- full_join(df1, chr3_final, by='FID')
df3 <- full_join(df2, chr4_final, by='FID')
df4 <- full_join(df3, chr5_final, by='FID')
df5 <- full_join(df4, chr6_final, by='FID')
df6 <- full_join(df5, chr7_final, by='FID')
df7 <- full_join(df6, chr8_final, by='FID')
df8 <- full_join(df7, chr9_final, by='FID')
df9 <- full_join(df8, chr10_final, by='FID')
df10 <- full_join(df9, chr11_final, by='FID')
df11 <- full_join(df10, chr12_final, by='FID')
df12 <- full_join(df11, chr13_final, by='FID')
df13 <- full_join(df12, chr14_final, by='FID')
df14 <- full_join(df13, chr15_final, by='FID')
df15 <- full_join(df14, chr16_final, by='FID')
df16 <- full_join(df15, chr17_final, by='FID')
df17 <- full_join(df16, chr18_final, by='FID')
df18 <- full_join(df17, chr19_final, by='FID')
df19 <- full_join(df18, chr20_final, by='FID')
df20 <- full_join(df19, chr21_final, by='FID')
df21 <- full_join(df20, chr22_final, by='FID')
all_chr_df <- full_join(df21, chrX_final, by='FID')

dim(all_chr_df)
# 327272     47

# next, sum all of the Sums columns

# changing DIF to character for now so that I don't accidentally do math on it!
all_chr_df$FID <- as.character(all_chr_df$FID)

#sum all Sums columns
all_chr_df$SUMS <- rowSums(all_chr_df[, c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46)], na.rm = TRUE)
unique(all_chr_df$SUMS)
# 1 2 5 3 4 6 7 8 9

# can remove the Sums columns now
all_chr_df <- all_chr_df[,-c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46)]


# compile all of the VAR columns
library(tidyr)
all_chr_df <- unite(all_chr_df, "Var", c(2:24), sep = ";", na.rm = TRUE)
# looks good :)

# save to R
write.csv(all_chr_df, "all_chr_df.csv")
# save this for later
system('dx upload all_chr_df.csv --path Emily-folder/pathvar/all_chr_df.csv')
















# annotate by GENE --------------------------------------------------------


library(tidyr)
library(readr)
library(dplyr)

# load in the all_chr_df.csv file and remove the first column
system('dx download Emily-folder/pathvar/all_chr_df.csv')
all_chr_df <- read_csv("all_chr_df.csv")
all_chr_df <- all_chr_df[-1]





# Can't use grep to find variant locations - as it could mess up with the chr numbers in the 10s and 20s
# so need to split the variant data into one cell per variant

unique(all_chr_df$SUMS)
# up to 9 variables in a cell
all_chr_df <- separate(data = all_chr_df, col = Var, into = c("Var1","Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9"), sep = ";")
# check
unique(all_chr_df$Var9)
# there are two variants in this column - so I've split into just the right amount! :)

# next: add gene name into each gene name column based on these variants







# # testing things out before I apply it to the code
# 
# # trying to merge/combine two dfs with possible duplicates in the reference
# # make sure variants are separated into one column each
# df_1 <- data.frame(id = c("id1", "id2","id3","id4","id5","id6","id7","id8"),
#                    Var1 = c("1:1:A:A_A","1:2:G:G_G","3:3:C:C_C","1:1:A:A_A","1:1:A:A_A","1:1:A:A_A","1:1:A:A_A","1:1:A:A_A"),
#                    Var2 = c(NA, "1:1:A:A_A", "1:1:A:A_A", "1:2:G:G_G", NA, NA, "4:4:T:T_T", NA))
# 
# # make sure I have a df that contains all variants of interest and their associated genes!
# df_2 <- data.frame(Var1 = c("1:1:A:A_A", "1:2:G:G_G", "3:3:C:C_C", "4:4:T:T_T"),
#                    GENE = c("GENE1", "GENE2", "GENE3", "GENE4"))
# 
# # can merge by each column of variants - so will end up with 9 columns of genes that I can merge manually!
# 
# # left join to combine participant df with the gene reference df by the Var1 column (as ref)
# df_3 <- left_join(df_1, df_2, by = 'Var1')
# # make sure to rename the variant column in the GENE ref df
# df_2 <- rename(df_2, Var2 = Var1)
# # repeat left join using the second variant column as reference!!!
# df_4 <- left_join(df_3, df_2, by = 'Var2')
# # repeat rename and left join for the rest of the variant columns
# df_2 <- rename(df_2, Var3 = Var2)
# df_5 <- left_join(df_4, df_2, by = 'Var3')
# #...






# so I have a plan for how to do it
# download var-genes data
system('dx download -r Emily-folder/variants-genes')

chr1_variants_genes <- read_csv("variants-genes/chr1_variants_genes.csv")
chr2_variants_genes <- read_csv("variants-genes/chr2_variants_genes.csv")
chr3_variants_genes <- read_csv("variants-genes/chr3_variants_genes.csv")
chr4_variants_genes <- read_csv("variants-genes/chr4_variants_genes.csv")
chr5_variants_genes <- read_csv("variants-genes/chr5_variants_genes.csv")
chr6_variants_genes <- read_csv("variants-genes/chr6_variants_genes.csv")
chr7_variants_genes <- read_csv("variants-genes/chr7_variants_genes.csv")
chr8_variants_genes <- read_csv("variants-genes/chr8_variants_genes.csv")
chr9_variants_genes <- read_csv("variants-genes/chr9_variants_genes.csv")
chr10_variants_genes <- read_csv("variants-genes/chr10_variants_genes.csv")
chr11_variants_genes <- read_csv("variants-genes/chr11_variants_genes.csv")
chr12_variants_genes <- read_csv("variants-genes/chr12_variants_genes.csv")
chr13_variants_genes <- read_csv("variants-genes/chr13_variants_genes.csv")
chr14_variants_genes <- read_csv("variants-genes/chr14_variants_genes.csv")
chr15_variants_genes <- read_csv("variants-genes/chr15_variants_genes.csv")
chr16_variants_genes <- read_csv("variants-genes/chr16_variants_genes.csv")
chr17_variants_genes <- read_csv("variants-genes/chr17_variants_genes.csv")
chr18_variants_genes <- read_csv("variants-genes/chr18_variants_genes.csv")
chr19_variants_genes <- read_csv("variants-genes/chr19_variants_genes.csv")
chr20_variants_genes <- read_csv("variants-genes/chr20_variants_genes.csv")
chr21_variants_genes <- read_csv("variants-genes/chr21_variants_genes.csv")
chr22_variants_genes <- read_csv("variants-genes/chr22_variants_genes.csv")
chrX_variants_genes <- read_csv("variants-genes/chrX_variants_genes.csv")

# remove first column from each
chr1_variants_genes <- chr1_variants_genes[-1]
chr2_variants_genes <- chr2_variants_genes[-1]
chr3_variants_genes <- chr3_variants_genes[-1]
chr4_variants_genes <- chr4_variants_genes[-1]
chr5_variants_genes <- chr5_variants_genes[-1]
chr6_variants_genes <- chr6_variants_genes[-1]
chr7_variants_genes <- chr7_variants_genes[-1]
chr8_variants_genes <- chr8_variants_genes[-1]
chr9_variants_genes <- chr9_variants_genes[-1]
chr10_variants_genes <- chr10_variants_genes[-1]
chr11_variants_genes <- chr11_variants_genes[-1]
chr12_variants_genes <- chr12_variants_genes[-1]
chr13_variants_genes <- chr13_variants_genes[-1]
chr14_variants_genes <- chr14_variants_genes[-1]
chr15_variants_genes <- chr15_variants_genes[-1]
chr16_variants_genes <- chr16_variants_genes[-1]
chr17_variants_genes <- chr17_variants_genes[-1]
chr18_variants_genes <- chr18_variants_genes[-1]
chr19_variants_genes <- chr19_variants_genes[-1]
chr20_variants_genes <- chr20_variants_genes[-1]
chr21_variants_genes <- chr21_variants_genes[-1]
chr22_variants_genes <- chr22_variants_genes[-1]
chrX_variants_genes <- chrX_variants_genes[-1]

# want to combine the var-gene reference docs together (instead of by chr)
# before merging them all together - need to change X chr from 'X' to '23'
chrX_variants_genes$variant <- strsplit(chrX_variants_genes$variant, ":") %>% 
  lapply(., function(x) paste0("23", ":", x[2], ":", x[3], ":", x[4])) %>% 
  unlist()
# looks good :)

# now to combine them all together!
# adding new observations, so should use rbind
var_genes_1 <- rbind(chr1_variants_genes, chr2_variants_genes)
var_genes_2 <- rbind(var_genes_1, chr2_variants_genes)
var_genes_3 <- rbind(var_genes_2, chr3_variants_genes)
var_genes_4 <- rbind(var_genes_3, chr4_variants_genes)
var_genes_5 <- rbind(var_genes_4, chr5_variants_genes)
var_genes_6 <- rbind(var_genes_5, chr6_variants_genes)
var_genes_7 <- rbind(var_genes_6, chr7_variants_genes)
var_genes_8 <- rbind(var_genes_7, chr8_variants_genes)
var_genes_9 <- rbind(var_genes_8, chr9_variants_genes)
var_genes_10 <- rbind(var_genes_9, chr10_variants_genes)
var_genes_11 <- rbind(var_genes_10, chr11_variants_genes)
var_genes_12 <- rbind(var_genes_11, chr12_variants_genes)
var_genes_13 <- rbind(var_genes_12, chr13_variants_genes)
var_genes_14 <- rbind(var_genes_13, chr14_variants_genes)
var_genes_15 <- rbind(var_genes_14, chr15_variants_genes)
var_genes_16 <- rbind(var_genes_15, chr16_variants_genes)
var_genes_17 <- rbind(var_genes_16, chr17_variants_genes)
var_genes_18 <- rbind(var_genes_17, chr18_variants_genes)
var_genes_19 <- rbind(var_genes_18, chr19_variants_genes)
var_genes_20 <- rbind(var_genes_19, chr20_variants_genes)
var_genes_21 <- rbind(var_genes_20, chr21_variants_genes)
var_genes_22 <- rbind(var_genes_21, chr22_variants_genes)
var_genes_X <- rbind(var_genes_22, chrX_variants_genes)
# total of 61249 variants



# now to use the code I wrote above to add gene names to all_chr_df 
# need var_genes_X to have the same Var name as all_chr_df 
var_genes_X <- rename(var_genes_X, Var1 = variant, GENE = GENE_NAME)


# now to annotate the dataframe!
# left join to combine participant df with the gene reference df by the Var1 column (as ref)
all_chr_df_1 <- left_join(all_chr_df, var_genes_X, by = 'Var1')

# make sure to rename the variant column in the GENE ref df
var_genes_X <- rename(var_genes_X, Var2 = Var1)

# repeat!
all_chr_df_2 <- left_join(all_chr_df_1, var_genes_X, by = 'Var2')

var_genes_X <- rename(var_genes_X, Var3 = Var2)
all_chr_df_3 <- left_join(all_chr_df_2, var_genes_X, by = 'Var3')

var_genes_X <- rename(var_genes_X, Var4 = Var3)
all_chr_df_4 <- left_join(all_chr_df_3, var_genes_X, by = 'Var4')

var_genes_X <- rename(var_genes_X, Var5 = Var4)
all_chr_df_5 <- left_join(all_chr_df_4, var_genes_X, by = 'Var5')

var_genes_X <- rename(var_genes_X, Var6 = Var5)
all_chr_df_6 <- left_join(all_chr_df_5, var_genes_X, by = 'Var6')

var_genes_X <- rename(var_genes_X, Var7 = Var6)
all_chr_df_7 <- left_join(all_chr_df_6, var_genes_X, by = 'Var7')

var_genes_X <- rename(var_genes_X, Var8 = Var7)
all_chr_df_8 <- left_join(all_chr_df_7, var_genes_X, by = 'Var8')

var_genes_X <- rename(var_genes_X, Var9 = Var8)
all_chr_df_9 <- left_join(all_chr_df_8, var_genes_X, by = 'Var9')

# has increased the number of observations
duplicated(all_chr_df_9)
# some of the observations have been duplicated!
annotated_1 <- all_chr_df_9[!duplicated(all_chr_df_9),]
# removed the duplicate rows - the correct number of observations again :)


annotated_1 <- rename(annotated_1, 
                      GENE9 = GENE, 
                      GENE1 = GENE.x,
                      GENE2 = GENE.y,
                      GENE3 = GENE.x.x,
                      GENE4 = GENE.y.y,
                      GENE5 = GENE.x.x.x,
                      GENE6 = GENE.y.y.y,
                      GENE7 = GENE.x.x.x.x,
                      GENE8 = GENE.y.y.y.y)


# check to make sure that all variants have been annotated with a gene name
annotated_1[!is.na(annotated_1$Var1) & is.na(annotated_1$GENE1), ]
annotated_1[!is.na(annotated_1$Var2) & is.na(annotated_1$GENE2), ]
annotated_1[!is.na(annotated_1$Var3) & is.na(annotated_1$GENE3), ]
annotated_1[!is.na(annotated_1$Var4) & is.na(annotated_1$GENE4), ]
annotated_1[!is.na(annotated_1$Var5) & is.na(annotated_1$GENE5), ]
annotated_1[!is.na(annotated_1$Var6) & is.na(annotated_1$GENE6), ]
annotated_1[!is.na(annotated_1$Var7) & is.na(annotated_1$GENE7), ]
annotated_1[!is.na(annotated_1$Var8) & is.na(annotated_1$GENE8), ]
annotated_1[!is.na(annotated_1$Var9) & is.na(annotated_1$GENE9), ]
# all give tibble: 0 x 20
# so there are no variables that have not been annotated with a gene name!!! :)






# leave VAR and GENE in separate columns for now, may make further annotation easier!













# annotate by phenotype and cancer association ----------------------------

# load in gene_phenotype_cancer.csv
system('dx download Emily-folder/gene_phenotype_cancer.csv')

gene_phenotype_cancer <- read_csv("gene_phenotype_cancer.csv")

# annotate the table in the same way I annotated with GENE
# using GENE as a reference this time
# when its done I can merge all of the fields together

# check all of the fields in this datafield (go thr all the column names)
unique(gene_phenotype_cancer$non_LOF_var)
# looks right to me!



# now to use the code I wrote above to add phenotype and cancer association data 
# need gene_phenotype_cancer to have the same GENE name as annotated_1 

gene_phenotype_cancer <- rename(gene_phenotype_cancer, GENE1 = GENE)

# now to annotate the dataframe!

annotated_1a <- left_join(annotated_1, gene_phenotype_cancer, by = 'GENE1')

gene_phenotype_cancer <- rename(gene_phenotype_cancer, GENE2 = GENE1)
annotated_1b <- left_join(annotated_1a, gene_phenotype_cancer, by = 'GENE2')

gene_phenotype_cancer <- rename(gene_phenotype_cancer, GENE3 = GENE2)
annotated_1c <- left_join(annotated_1b, gene_phenotype_cancer, by = 'GENE3')

gene_phenotype_cancer <- rename(gene_phenotype_cancer, GENE4 = GENE3)
annotated_1d <- left_join(annotated_1c, gene_phenotype_cancer, by = 'GENE4')

gene_phenotype_cancer <- rename(gene_phenotype_cancer, GENE5 = GENE4)
annotated_1e <- left_join(annotated_1d, gene_phenotype_cancer, by = 'GENE5')

gene_phenotype_cancer <- rename(gene_phenotype_cancer, GENE6 = GENE5)
annotated_1f <- left_join(annotated_1e, gene_phenotype_cancer, by = 'GENE6')

gene_phenotype_cancer <- rename(gene_phenotype_cancer, GENE7 = GENE6)
annotated_1g <- left_join(annotated_1f, gene_phenotype_cancer, by = 'GENE7')

gene_phenotype_cancer <- rename(gene_phenotype_cancer, GENE8 = GENE7)
annotated_1h <- left_join(annotated_1g, gene_phenotype_cancer, by = 'GENE8')

gene_phenotype_cancer <- rename(gene_phenotype_cancer, GENE9 = GENE8)
annotated_2 <- left_join(annotated_1h, gene_phenotype_cancer, by = 'GENE9')

# correct number of observations!

# need to check that all of the genes match the two dataframes 
# make a list of all the unique GENE names in annotated_2 
GENEs_in_a2 <-c(unique(annotated_2$GENE1), unique(annotated_2$GENE2), unique(annotated_2$GENE3), 
                unique(annotated_2$GENE4), unique(annotated_2$GENE5), unique(annotated_2$GENE6), 
                unique(annotated_2$GENE7), unique(annotated_2$GENE8), unique(annotated_2$GENE9))

GENEs_in_a2_unique <- unique(GENEs_in_a2)
# includes NA - 326 gene names and 1 NA

# use setdiff to determine if any GENE names are missing from gene_phenotype_cancer or don't match exactly
setdiff(GENEs_in_a2_unique, gene_phenotype_cancer$GENE9)
# NA
# cool - so all gene names are present and they've all been annotated!

setdiff(gene_phenotype_cancer$GENE9, GENEs_in_a2_unique)
# [1] "NOTCH2"   "LMNA"     "MACF1"    "CDKN1C"   "IFITM5"   "FAM111A"  "TUBA1A"   "PEX5"     "FBN1"     "KIF22"    "SRCAP"    "VPS4A"   
# [13] "SMARCE1"  "THRA"     "FZD2"     "SPOP"     "PRKAR1A"  "ACTG1"    "PIEZO2"   "LMNB2"    "RNU4ATAC" "DNMT3A"   "COL9A3"   "RAF1"    
# [25] "MRAS"     "RNF13"    "DVL3"     "PPP3CA"   "WDFY3"    "LMNB1"    "GMNN"     "RALA"     "COL1A2"   "FGFR1"    "PTDSS1"   "SMARCA2" 
# [37] "AIFM1"    "BGN"      "FLNA"     "MBTPS2"   "ARX"
# these are the genes that do not have variants of interest in the biobank


# now: can combine the fields together!
annotated_2 <- unite(annotated_2, "VAR", c(2:10), sep = ";", na.rm = TRUE)
annotated_2 <- unite(annotated_2, "GENE", c(4:12), sep = ";", na.rm = TRUE)
annotated_2 <- unite(annotated_2, "Short_Stature", c(5, 15, 25, 35, 45, 55, 65, 75, 85), sep = ";", na.rm = TRUE)
annotated_2 <- unite(annotated_2, "Microcephaly", c(6, 15, 24, 33, 42, 51, 60, 69, 78), sep = ";", na.rm = TRUE)
annotated_2 <- unite(annotated_2, "non_LOF_var", c(7, 15, 23, 31, 39, 47, 55, 63, 71), sep = ";", na.rm = TRUE)
annotated_2 <- unite(annotated_2, "GOF", c(8, 15, 22, 29, 36, 43, 50, 57, 64), sep = ";", na.rm = TRUE)
annotated_2 <- unite(annotated_2, "DomNeg", c(9, 15, 21, 27, 33, 39, 45, 51, 57), sep = ";", na.rm = TRUE)
annotated_2 <- unite(annotated_2, "Cluster_wo_GOF_DN", c(10, 15, 20, 25, 30, 35, 40, 45, 50), sep = ";", na.rm = TRUE)
annotated_2 <- unite(annotated_2, "Syndrome_associated_w_cancer", c(11, 15, 19, 23, 27, 31, 35, 39, 43), sep = ";", na.rm = TRUE)
annotated_2 <- unite(annotated_2, "gene_vars_associated_w_cancer_risk", c(12, 15, 18, 21, 24, 27, 30, 33, 36), sep = ";", na.rm = TRUE)
annotated_2 <- unite(annotated_2, "gene_associated_w_cancer", c(13, 15, 17, 19, 21, 23, 25, 27, 29), sep = ";", na.rm = TRUE)
annotated_2 <- unite(annotated_2, "Association_type", c(14:22), sep = ";", na.rm = TRUE)

colnames(annotated_2)


# now check unique values in each column!
unique(annotated_2[,10])
# there aren't any clustered but not GOF or DN gene variants in the dataframe
# should remove this field 
unique(annotated_2[,14])

path_var_data <- annotated_2[-10]



# change "" values in the non_LOF_var, GOF, and DomNeg columns to NA
# for some reason trying to change them doesn't work unless I change to factor first

## non_LOF_var
unique(path_var_data$non_LOF_var)
path_var_data$non_LOF_var <- as.factor(path_var_data$non_LOF_var)
levels(path_var_data$non_LOF_var)[levels(path_var_data$non_LOF_var)==""] <- "No"
path_var_data$non_LOF_var <- as.character(path_var_data$non_LOF_var)

## GOF
unique(path_var_data$GOF)
path_var_data$GOF <- as.factor(path_var_data$GOF)
levels(path_var_data$GOF)[levels(path_var_data$GOF)==""] <- "No"
path_var_data$GOF <- as.character(path_var_data$GOF)

## DomNeg
unique(path_var_data$DomNeg)
path_var_data$DomNeg <- as.factor(path_var_data$DomNeg)
levels(path_var_data$DomNeg)[levels(path_var_data$DomNeg)==""] <- "No"
path_var_data$DomNeg <- as.character(path_var_data$DomNeg)






# save to R
write.csv(path_var_data, "path_var_data.csv")
# save this for later
system('dx upload path_var_data.csv --path Emily-folder/pathvar/path_var_data.csv')



















# EntireCohort_pathvar.csv ------------------------------------------------

# only needs to be executed once

# load packages
library(data.table)
library(dplyr)

# load in the EntireCohort_ready.csv file
EntireCohort_ready <- fread('/mnt/project/Emily-folder/EntireCohort_ready.csv', data.table = FALSE)


# check the data!
colnames(EntireCohort_ready)
# remove col1
EntireCohort_ready <- EntireCohort_ready[-1]


# add the pathvar data, by eid
# load in the dataframe
path_var_data <- fread('/mnt/project/Emily-folder/pathvar/path_var_data.csv', data.table = FALSE)
# The dataframe gets loaded in with an extra column called 'V1' which I don't want
path_var_data <- path_var_data[-1]
# removed :)

# make sure their eid columns have the same name
names(path_var_data)[names(path_var_data) == 'FID'] <- 'eid'

# join the two dataframes together
EntireCohort_pathvar <- left_join(EntireCohort_ready, path_var_data, by="eid")
# check
colnames(EntireCohort_pathvar)
summary(EntireCohort_pathvar)
unique(EntireCohort_pathvar$GENE)
unique(EntireCohort_pathvar$Short_Stature)
unique(EntireCohort_pathvar$Microcephaly)
unique(EntireCohort_pathvar$Association_type)

EntireCohort_pathvar[EntireCohort_pathvar$Association_type == "" ,]
# NA eids - likely redacted data




# clean up the data


unique(EntireCohort_pathvar$non_LOF_var)
# "No"  NA    "Yes"
EntireCohort_pathvar$non_LOF_var[is.na(EntireCohort_pathvar$non_LOF_var)] <- "No"
# "No"  "Yes"

unique(EntireCohort_pathvar$GOF)
# "No"        NA          "Yes"       "Suspected"
EntireCohort_pathvar$GOF[is.na(EntireCohort_pathvar$GOF)] <- "No"
# "No"        "Yes"       "Suspected"

unique(EntireCohort_pathvar$DomNeg)
# "No"        NA          "Yes"       "Suspected"
EntireCohort_pathvar$DomNeg[is.na(EntireCohort_pathvar$DomNeg)] <- "No"
# "No"        "Yes"       "Suspected"



# save EntireCohort_pathvar so that I dont need to execute this code each time:
#save it to R
write.csv(EntireCohort_pathvar, "EntireCohort_pathvar.csv")
# now to save this to UKB project
system('dx upload EntireCohort_pathvar.csv --path Emily-folder/pathvar/EntireCohort_pathvar.csv')

















# Save the script!!! -----------------------------------------------------------


system('dx upload EntireCohort_pathvar_script.R --path Emily-folder/pathvar/EntireCohort_pathvar_script.R')
