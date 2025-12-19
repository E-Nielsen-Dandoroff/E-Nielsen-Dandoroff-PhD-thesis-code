## 2023.09.15
## processing variant scores




# data cleaning -----------------------------------------------------------



library(readr)


# made a list of all of the .csv files i want to import
data_files <- list.files(pattern = "*output.csv")








# trialling Mik's code to make a list of dfs:
csvList = list()
for (i in 1:311) {
  csvList[[i]] = read_delim(data_files[i], delim = "\t", escape_double = FALSE, trim_ws = TRUE)
}


## separate out all the dataframes!
library(purrr)
library(dplyr)
library(stringr)

## here is the template code - it gives all the dataframes numbers from 1:311, but I want them to have the gene names!
# imap(csvList, ~ set_names(tibble(.x), .y)) %>% 
#   set_names(str_c("DF", 1:311)) %>% 
#   list2env(.GlobalEnv)


# this code removes '.csv' from each of the data_files names
gene_names <- tools::file_path_sans_ext(data_files)


## this code extracts all of the dataframes from the list and names them by gene!
imap(csvList, ~ set_names(tibble(.x), .y)) %>% 
  set_names(gene_names) %>% 
  list2env(.GlobalEnv)

## this code doesn't retain the column names however - I'll add them back later :)
#c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")






## combine by chr
chr1_genes <- list.files(pattern = "*_1_output.csv")
chr1_genes <- tools::file_path_sans_ext(chr1_genes)
chr1_genes
chr1_list <- list(end_AMPD2_1_output, end_ARHGEF2_1_output, end_ASPM_1_output, 
end_ATAD3A_1_output, end_B3GALT6_1_output, end_CENPF_1_output, end_COL11A1_1_output, 
end_DDR2_1_output, end_DHCR24_1_output, end_GNPAT_1_output, end_GORAB_1_output, 
end_HHAT_1_output, end_HSPG2_1_output, end_KIF14_1_output, end_LBR_1_output, 
end_LHX4_1_output, end_MFSD2A_1_output, end_NUP133_1_output, end_ORC1_1_output, 
end_P3H1_1_output, end_PHGDH_1_output, end_SASS6_1_output, end_SLC35D1_1_output, 
end_STIL_1_output, end_TAF13_1_output, end_TBCE_1_output, end_TBX15_1_output, 
end_TOE1_1_output, end_TSEN15_1_output, end_YARS1_1_output, end_YRDC_1_output, 
end_ZMPSTE24_1_output)


chr1 <- data.table::rbindlist(chr1_list, use.names = FALSE)
colnames(chr1) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr1$BayesDel_addAF_score <- as.numeric(chr1$BayesDel_addAF_score)

summary(chr1)
# BayesDel has 262 NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
# far too many variants to go through and manually remove duplicates in diff transcripts

#to check that the correct genes are represented!
unique(chr1$genename)                                                                           
# none of the gene names are missing :)
# genenames that are present that shouldn't be!                                                                                  
# [83] "GGPS1;GGPS1;GGPS1;GGPS1;GGPS1"                                                       
# [84] "GGPS1;GGPS1;GGPS1;GGPS1;GGPS1;GGPS1"                                                 
# [85] "GGPS1;GGPS1;GGPS1;GGPS1"                                                                             
# [94] "B3GALNT2"                                                                            
# [98] "MUTYH;MUTYH;MUTYH;MUTYH;MUTYH"                                                       
# will need to write code to extract the rows with these genename values!                                                                               
      
chr1 <- chr1[!chr1$genename %in% c("GGPS1;GGPS1;GGPS1;GGPS1;GGPS1", "GGPS1;GGPS1;GGPS1;GGPS1;GGPS1;GGPS1", 
                                   "GGPS1;GGPS1;GGPS1;GGPS1", "B3GALNT2", "MUTYH;MUTYH;MUTYH;MUTYH;MUTYH")]


# chr1 is done (i think) :)

write.csv(chr1, "chr1.csv")












# now to repeat this for all of the other genes!

## chromosome 2
chr2_genes <- list.files(pattern = "*_2_output.csv")
chr2_genes <- tools::file_path_sans_ext(chr2_genes)
chr2_genes
chr2_list <- list(end_AGPS_2_output, end_BCS1L_2_output, end_BUB1_2_output, 
end_CKAP2L_2_output, end_CRIPT_2_output, end_DYNC1I2_2_output, end_EIF2AK3_2_output, 
end_ERCC3_2_output, end_EXOC6B_2_output, end_FN1_2_output, end_GLI2_2_output, 
end_IHH_2_output, end_MATN3_2_output, end_MYCN_2_output, end_NCAPH_2_output, 
end_NHEJ1_2_output, end_NPPC_2_output, end_OBSL1_2_output, end_ORC4_2_output, 
end_SLC1A4_2_output, end_SMARCAL1_2_output, end_SMPD4_2_output, end_SOX11_2_output, 
end_SPRED2_2_output, end_STAMBP_2_output, end_TPRKB_2_output, end_WDR35_2_output)



chr2 <- data.table::rbindlist(chr2_list, use.names = FALSE)
colnames(chr2) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr2$BayesDel_addAF_score <- as.numeric(chr2$BayesDel_addAF_score)

summary(chr2)
# BayesDel has 1 NA value

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr2)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr2$genename)             
# genenames that are present that shouldn't be!                                                                                                                                                    
# [3] "TTC30B"                                                                                                                                                                           
# [85] "MZT2B"                                                                                                                                                  
# [86] "MZT2B;MZT2B"                                                                                                                                            
# [87] "MZT2B;MZT2B;MZT2B"                                                         
# will need to write code to extract the rows with these genename values!                                                                               

chr2 <- chr2[!chr2$genename %in% c("TTC30B", "MZT2B", "MZT2B;MZT2B", "MZT2B;MZT2B;MZT2B")]


# chr2 is done (i think) :)

write.csv(chr2, "chr2.csv")













## chromosome 3
chr3_genes <- list.files(pattern = "*_3_output.csv")
chr3_genes <- tools::file_path_sans_ext(chr3_genes)
chr3_genes
chr3_list <- list(end_ATR_3_output, end_ATRIP_3_output, end_CEP63_3_output, end_COPB2_3_output, 
end_CRTAP_3_output, end_FLNB_3_output, end_GHSR_3_output, end_GLB1_3_output, end_HESX1_3_output, 
end_IFT122_3_output, end_NEPRO_3_output, end_PCYT1A_3_output, end_PDCD6IP_3_output, 
end_POC1A_3_output, end_POU1F1_3_output, end_PTH1R_3_output, end_RASA2_3_output, 
end_TBC1D23_3_output, end_TRAIP_3_output, end_TSEN2_3_output, end_WNT5A_3_output)



chr3 <- data.table::rbindlist(chr3_list, use.names = FALSE)
colnames(chr3) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr3$BayesDel_addAF_score <- as.numeric(chr3$BayesDel_addAF_score)

summary(chr3)
# BayesDel has no NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr3)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr3$genename)             
# genenames that are present that shouldn't be!                                                      
# [7] "TREX1"                                                                
# [8] "TREX1;TREX1"                                                          
# [9] "TREX1;TREX1;TREX1;TREX1;TREX1;TREX1"                                                    
# [16] "MRPS22;MRPS22;MRPS22;MRPS22"                                                 
# [38] "TMPPE;TMPPE"                                                          
# [39] "TMPPE"                                                                 
# [55] "SLC51A"                                                               
# [56] "SLC51A;SLC51A"                                                                                                 
# will need to write code to extract the rows with these genename values!                                                                               

chr3 <- chr3[!chr3$genename %in% c("TREX1", "TREX1;TREX1", "TREX1;TREX1;TREX1;TREX1;TREX1;TREX1", 
                                   "MRPS22;MRPS22;MRPS22;MRPS22", "TMPPE;TMPPE", 
                                   "TMPPE", "SLC51A", "SLC51A;SLC51A")]


# chr3 is done (i think) :)

write.csv(chr3, "chr3.csv")










## chromosome 4
chr4_genes <- list.files(pattern = "*_4_output.csv")
chr4_genes <- tools::file_path_sans_ext(chr4_genes)
chr4_genes

# need to change end_NKX3-2_4_output > end_NKX3_2_4_output
end_NKX3_2_4_output <- `end_NKX3-2_4_output`

chr4_list <- list(end_BMPR1B_4_output, end_CENPE_4_output, end_CEP135_4_output, 
end_EVC_4_output, end_EVC2_4_output, end_LARP7_4_output, end_MSMO1_4_output, 
end_NEK1_4_output, end_NKX3_2_4_output, end_PLK4_4_output, end_PRKG2_4_output, 
end_RAB33B_4_output, end_SEPSECS_4_output, end_SLC10A7_4_output, end_SLC39A8_4_output, 
end_TMEM165_4_output, end_TRMT10A_4_output, end_UFSP2_4_output, end_WDR19_4_output)



chr4 <- data.table::rbindlist(chr4_list, use.names = FALSE)
colnames(chr4) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr4$BayesDel_addAF_score <- as.numeric(chr4$BayesDel_addAF_score)

summary(chr4)
# BayesDel has no NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr4)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr4$genename)             
# genenames that are present that shouldn't be!                                      
# [50] "CLOCK;CLOCK;CLOCK"                                                                                                 
# will need to write code to extract the rows with these genename values!                                                                               

chr4 <- chr4[!chr4$genename %in% "CLOCK;CLOCK;CLOCK"]


# chr4 is done (i think) :)

write.csv(chr4, "chr4.csv")










## chromosome 5
chr5_genes <- list.files(pattern = "*_5_output.csv")
chr5_genes <- tools::file_path_sans_ext(chr5_genes)
chr5_genes
chr5_list <- list(end_ADAMTS2_5_output, end_ARSB_5_output, end_B4GALT7_5_output, 
end_ERCC8_5_output, end_GHR_5_output, end_LIFR_5_output, end_NIPBL_5_output, 
end_OCLN_5_output, end_PCDH12_5_output, end_PIK3R1_5_output, end_PROP1_5_output, 
end_RAD50_5_output, end_SLC26A2_5_output, end_THG1L_5_output, end_XRCC4_5_output)



chr5 <- data.table::rbindlist(chr5_list, use.names = FALSE)
colnames(chr5) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr5$BayesDel_addAF_score <- as.numeric(chr5$BayesDel_addAF_score)

summary(chr5)
# BayesDel has 46 NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr5)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr5$genename)             
# genenames that are present that shouldn't be!                                                                              
# [37] "IL5"                                                                                                                                       
# will need to write code to extract the rows with these genename values!                                                                               

chr5 <- chr5[!chr5$genename %in% "IL5"]


# chr5 is done (i think) :)

write.csv(chr5, "chr5.csv")










## chromosome 6
chr6_genes <- list.files(pattern = "*_6_output.csv")
chr6_genes <- tools::file_path_sans_ext(chr6_genes)
chr6_genes
chr6_list <- list(end_CCN6_6_output, end_CDC40_6_output, end_COL10A1_6_output, 
end_COL11A2_6_output, end_CUL7_6_output, end_ICK_6_output, end_MCM3_6_output, 
end_PEX7_6_output, end_PPIL1_6_output, end_TUBB2B_6_output)



chr6 <- data.table::rbindlist(chr6_list, use.names = FALSE)
colnames(chr6) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr6$BayesDel_addAF_score <- as.numeric(chr6$BayesDel_addAF_score)

summary(chr6)
# BayesDel has 1 NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr6)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr6$genename)             
# genenames that are present that shouldn't be!    
# none :)
# will not need to write code to extract the rows with these genename values!                                                                               

# chr6 <- chr6[!chr6$genename %in% ]


# chr6 is done (i think) :)

write.csv(chr6, "chr6.csv")












## chromosome 7
chr7_genes <- list.files(pattern = "*_7_output.csv")
chr7_genes <- tools::file_path_sans_ext(chr7_genes)
chr7_genes
chr7_list <- list(end_ACTB_7_output, end_ASNS_7_output, end_CDK6_7_output, 
end_FAM20C_7_output, end_GHRHR_7_output, end_GLI3_7_output, end_GUSB_7_output, 
end_KDELR2_7_output, end_MCM7_7_output, end_PCLO_7_output, end_SBDS_7_output, 
end_SHH_7_output, end_TRAPPC14_7_output)



chr7 <- data.table::rbindlist(chr7_list, use.names = FALSE)
colnames(chr7) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr7$BayesDel_addAF_score <- as.numeric(chr7$BayesDel_addAF_score)

summary(chr7)
# BayesDel has no NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr7)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr7$genename)             
# genenames that are present that shouldn't be!                                                   
# none :)
# will not need to write code to extract the rows with these genename values!                                                                               

#chr7 <- chr7[!chr7$genename %in% ]


# chr7 is done (i think) :)

write.csv(chr7, "chr7.csv")











## chromosome 8
chr8_genes <- list.files(pattern = "*_8_output.csv")
chr8_genes <- tools::file_path_sans_ext(chr8_genes)
chr8_genes
chr8_list <- list(end_BMP1_8_output, end_BPNT2_8_output, end_ESCO2_8_output, 
end_MCPH1_8_output, end_NSMCE2_8_output, end_POP1_8_output, end_PUF60_8_output, 
end_RECQL4_8_output, end_TNFRSF11B_8_output, end_TONSL_8_output, end_TRPS1_8_output, 
end_RAD21_output)
# RAD21 is the only gene I forgot to add the gene number to the name!


chr8 <- data.table::rbindlist(chr8_list, use.names = FALSE)
colnames(chr8) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr8$BayesDel_addAF_score <- as.numeric(chr8$BayesDel_addAF_score)

summary(chr8)
# BayesDel has no NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr8)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr8$genename)             
# genenames that are present that shouldn't be!                                                       
# [43] "UTP23"                         
# will not need to write code to extract the rows with these genename values!                                                                               

chr8 <- chr8[!chr8$genename %in% "UTP23" ]


# chr8 is done (i think) :)

write.csv(chr8, "chr8.csv")













## chromosome 9
chr9_genes <- list.files(pattern = "*_9_output.csv")
chr9_genes <- tools::file_path_sans_ext(chr9_genes)
chr9_genes
chr9_list <- list(end_ADAMSTL2_9_output, end_CDK5RAP2_9_output, end_COL27A1_9_output, 
end_EXOSC3_9_output, end_LHX3_9_output, end_NANS_9_output, end_NPR2_9_output, 
end_PSAT1_9_output, end_ROR2_9_output, end_TMEM38B_9_output)



chr9 <- data.table::rbindlist(chr9_list, use.names = FALSE)
colnames(chr9) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr9$BayesDel_addAF_score <- as.numeric(chr9$BayesDel_addAF_score)

summary(chr9)
# BayesDel has no NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr9)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr9$genename)             
# genenames that are present that shouldn't be!                                          
# [18] "SPAG8"                                      
# will not need to write code to extract the rows with these genename values!                                                                               

chr9 <- chr9[!chr9$genename %in% "SPAG8" ]


# chr9 is done (i think) :)

write.csv(chr9, "chr9.csv")












## chromosome 10
chr10_genes <- list.files(pattern = "*_10_output.csv")
chr10_genes <- tools::file_path_sans_ext(chr10_genes)
chr10_genes
chr10_list <- list(end_ALDH18A1_10_output, end_CHST3_10_output, end_DNA2_10_output, 
end_ERCC6_10_output, end_FRMD4A_10_output, end_KIF11_10_output, end_KIFBP_10_output, 
end_MINPP1_10_output, end_PAPSS2_10_output, end_SMC3_10_output)



chr10 <- data.table::rbindlist(chr10_list, use.names = FALSE)
colnames(chr10) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr10$BayesDel_addAF_score <- as.numeric(chr10$BayesDel_addAF_score)

summary(chr10)
# BayesDel has 7 NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr10)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr10$genename)             
# genenames that are present that shouldn't be!                                                              
# none                                  
# will not need to write code to extract the rows with these genename values!                                                                               

#chr10 <- chr10[!chr10$genename %in%  ]


# chr10 is done (i think) :)

write.csv(chr10, "chr10.csv")












## chromosome 11
chr11_genes <- list.files(pattern = "*_11_output.csv")
chr11_genes <- tools::file_path_sans_ext(chr11_genes)
chr11_genes
chr11_list <- list(end_ARCN1_11_output, end_B3GAT3_11_output, end_CLP1_11_output, 
end_CREB3L1_11_output, end_DHCR7_11_output, end_DYNC2H1_11_output, end_INPPL1_11_output, 
end_LTBP3_11_output, end_MMP13_11_output, end_NCAPD3_11_output, end_SERPINH1_11_output, 
end_SIK3_11_output, end_SLC35C1_11_output, end_SLC39A13_11_output, end_TALDO1_11_output)



chr11 <- data.table::rbindlist(chr11_list, use.names = FALSE)
colnames(chr11) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr11$BayesDel_addAF_score <- as.numeric(chr11$BayesDel_addAF_score)

summary(chr11)
# BayesDel has 5 NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr11)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr11$genename)             
# genenames that are present that shouldn't be!                                                                                                          
# [47] "VPS26B;VPS26B"                                                 
# will not need to write code to extract the rows with these genename values!                                                                               

chr11 <- chr11[!chr11$genename %in% "VPS26B;VPS26B" ]


# chr11 is done (i think) :)

write.csv(chr11, "chr11.csv")














## chromosome 12
chr12_genes <- list.files(pattern = "*_12_output.csv")
chr12_genes <- tools::file_path_sans_ext(chr12_genes)
chr12_genes
chr12_list <- list(end_ANKLE2_12_output, end_CIT_12_output, end_COL2A1_12_output, 
end_GLI1_12_output, end_GNPTAB_12_output, end_HMGA2_12_output, end_IGF1_12_output, 
end_LEMD3_12_output, end_NCAPD2_12_output, end_NUP37_12_output, end_PHC1_12_output, 
end_POLE_12_output, end_PPFIBP1_12_output, end_PRIM1_12_output, end_PTHLH_12_output, 
end_SP7_12_output, end_TBX3_12_output, end_WNT1_12_output)



chr12 <- data.table::rbindlist(chr12_list, use.names = FALSE)
colnames(chr12) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr12$BayesDel_addAF_score <- as.numeric(chr12$BayesDel_addAF_score)

summary(chr12)
# BayesDel has no NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr12)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr12$genename)             
# genenames that are present that shouldn't be!                                            
# none                                            
# will not need to write code to extract the rows with these genename values!                                                                               

#chr12 <- chr12[!chr12$genename %in%  ]


# chr12 is done (i think) :)

write.csv(chr12, "chr12.csv")












## chromosome 13
chr13_genes <- list.files(pattern = "*_13_output.csv")
chr13_genes <- tools::file_path_sans_ext(chr13_genes)
chr13_genes
chr13_list <- list(end_B3GLCT_13_output, end_CENPJ_13_output, end_GPC6_13_output, end_LIG4_13_output)



chr13 <- data.table::rbindlist(chr13_list, use.names = FALSE)
colnames(chr13) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr13$BayesDel_addAF_score <- as.numeric(chr13$BayesDel_addAF_score)

summary(chr13)
# BayesDel has no NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr13)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr13$genename)             
# genenames that are present that shouldn't be!                                            
#none                                                    
# will not need to write code to extract the rows with these genename values!                                                                               

#chr13 <- chr13[!chr13$genename %in%  ]


# chr13 is done (i think) :)

write.csv(chr13, "chr13.csv")















## chromosome 14
chr14_genes <- list.files(pattern = "*_14_output.csv")
chr14_genes <- tools::file_path_sans_ext(chr14_genes)
chr14_genes
chr14_list <- list(end_FOXG1_14_output, end_GON7_14_output, end_GSC_14_output, 
end_NIN_14_output, end_OSGEP_14_output, end_OTX2_14_output, end_TRAPPC6B_14_output, 
end_TRIP11_14_output, end_VRK1_14_output)



chr14 <- data.table::rbindlist(chr14_list, use.names = FALSE)
colnames(chr14) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr14$BayesDel_addAF_score <- as.numeric(chr14$BayesDel_addAF_score)

summary(chr14)
# BayesDel has no NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr14)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr14$genename)             
# genenames that are present that shouldn't be!                                  
#none                                               
# will not need to write code to extract the rows with these genename values!                                                                               

#chr14 <- chr14[!chr14$genename %in%  ]


# chr14 is done (i think) :)

write.csv(chr14, "chr14.csv")













## chromosome 15
chr15_genes <- list.files(pattern = "*_15_output.csv")
chr15_genes <- tools::file_path_sans_ext(chr15_genes)
chr15_genes
chr15_list <- list(end_ACAN_15_output, end_BUB1B_15_output, end_CEP152_15_output, 
end_IFG1R_15_output, end_KNL1_15_output, end_MTHFS_15_output, end_RECQL3_15_output, 
end_TUBGCP4_15_output, end_UBR1_15_output, end_WDR73_15_output)



chr15 <- data.table::rbindlist(chr15_list, use.names = FALSE)
colnames(chr15) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr15$BayesDel_addAF_score <- as.numeric(chr15$BayesDel_addAF_score)

summary(chr15)
# BayesDel has no NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr15)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr15$genename)             
# genenames that are present that shouldn't be!                                                    
# [30] "TP53BP1;TP53BP1;TP53BP1;TP53BP1"                                                 
# will not need to write code to extract the rows with these genename values!                                                                               

chr15 <- chr15[!chr15$genename %in% "TP53BP1;TP53BP1;TP53BP1;TP53BP1" ]


# chr15 is done (i think) :)

write.csv(chr15, "chr15.csv")













## chromosome 16
chr16_genes <- list.files(pattern = "*_16_output.csv")
chr16_genes <- tools::file_path_sans_ext(chr16_genes)
chr16_genes
chr16_list <- list(end_AARS1_16_output, end_CDT1_16_output, end_CENPT_16_output, 
end_CHMP1A_16_output, end_CREBBP_16_output, end_CTU2_16_output, end_FANCA_16_output, 
end_GALNS_16_output, end_GINS2_16_output, end_GINS3_16_output, end_IGFALS_16_output, 
end_KATNB1_16_output, end_NDE1_16_output, end_ORC6_16_output, end_PAM16_16_output, 
end_PRMT7_16_output, end_RPL13_16_output, end_RSPRY1_16_output, end_THOC6_16_output)



chr16 <- data.table::rbindlist(chr16_list, use.names = FALSE)
colnames(chr16) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr16$BayesDel_addAF_score <- as.numeric(chr16$BayesDel_addAF_score)

summary(chr16)
# BayesDel has 2 NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr16)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr16$genename)             
# genenames that are present that shouldn't be!                                               
# [14] "SPATA33"                                         
# [22] "ZNF276;ZNF276;ZNF276"                           
# [23] "ZNF276;ZNF276"                                  
# [24] "ZNF276"                                                                                      
# [51] "MYH11;MYH11"                                    
# [52] "MYH11;MYH11;MYH11;MYH11"                                                        
# will not need to write code to extract the rows with these genename values!                                                                               

chr16 <- chr16[!chr16$genename %in% c("SPATA33", "ZNF276;ZNF276;ZNF276", "ZNF276;ZNF276", 
                                      "ZNF276", "MYH11;MYH11", "MYH11;MYH11;MYH11;MYH11") ]


# chr16 is done (i think) :)

write.csv(chr16, "chr16.csv")











## chromosome 17
chr17_genes <- list.files(pattern = "*_17_output.csv")
chr17_genes <- tools::file_path_sans_ext(chr17_genes)
chr17_genes
chr17_list <- list(end_CANT1_17_output, end_CDC6_17_output, end_COASY_17_output, 
end_COL1A1_17_output, end_DPH1_17_output, end_EIF5A_17_output, end_FKBP10_17_output, 
end_GH1_17_output, end_NXN_17_output, end_PYCR1_17_output, end_SERPINF1_17_output, 
end_SLC25A19_17_output, end_SOX9_17_output, end_STAT5B_17_output, end_TOP3A_17_output, 
end_TRIM37_17_output, end_TSEN54_17_output)



chr17 <- data.table::rbindlist(chr17_list, use.names = FALSE)
colnames(chr17) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr17$BayesDel_addAF_score <- as.numeric(chr17$BayesDel_addAF_score)

summary(chr17)
# BayesDel has 2 NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr17)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr17$genename)             
# genenames that are present that shouldn't be!                                                                                                        
# [25] "OVCA2"                                                                               
# will not need to write code to extract the rows with these genename values!                                                                               

chr17 <- chr17[!chr17$genename %in% "OVCA2" ]


# chr17 is done (i think) :)

write.csv(chr17, "chr17.csv")













## chromosome 18
chr18_genes <- list.files(pattern = "*_18_output.csv")
chr18_genes <- tools::file_path_sans_ext(chr18_genes)
chr18_genes
chr18_list <- list(end_DYM_18_output, end_IER3IP1_18_output, end_RBBP8_18_output, end_RTTN_18_output)



chr18 <- data.table::rbindlist(chr18_list, use.names = FALSE)
colnames(chr18) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr18$BayesDel_addAF_score <- as.numeric(chr18$BayesDel_addAF_score)

summary(chr18)
# BayesDel has no NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr18)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr18$genename)             
# genenames that are present that shouldn't be!                                                      
#none                                                                              
# will not need to write code to extract the rows with these genename values!                                                                               

#chr18 <- chr18[!chr18$genename %in%  ]


# chr18 is done (i think) :)

write.csv(chr18, "chr18.csv")










## chromosome 19
chr19_genes <- list.files(pattern = "*_19_output.csv")
chr19_genes <- tools::file_path_sans_ext(chr19_genes)
chr19_genes
chr19_list <- list(end_ACP5_19_output, end_CCDC8_19_output, end_COMP_19_output, 
end_ERCC2_19_output, end_GPX4_19_output, end_INSR_19_output, end_PNKP_19_output, 
end_PNPLA6_19_output, end_TRMT1_19_output, end_WDR62_19_output)



chr19 <- data.table::rbindlist(chr19_list, use.names = FALSE)
colnames(chr19) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr19$BayesDel_addAF_score <- as.numeric(chr19$BayesDel_addAF_score)

summary(chr19)
# BayesDel has 2 NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr19)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr19$genename)             
# genenames that are present that shouldn't be!                                                                    
# [12] "KLC3;KLC3;KLC3;KLC3"                                                  
# [13] "KLC3;KLC3;KLC3"                                                                                                       
# will need to write code to extract the rows with these genename values!                                                                               

chr19 <- chr19[!chr19$genename %in% c("KLC3;KLC3;KLC3;KLC3", "KLC3;KLC3;KLC3") ]


# chr19 is done (i think) :)

write.csv(chr19, "chr19.csv")












## chromosome 20
chr20_genes <- list.files(pattern = "*_20_output.csv")
chr20_genes <- tools::file_path_sans_ext(chr20_genes)
chr20_genes
chr20_list <- list(end_BMP2_20_output, end_CTSA_20_output, end_DDRGK1_20_output, 
end_GDF5_20_output, end_GHRH_20_output, end_GNAS_20_output, end_TP53RK_20_output, 
end_ZNF335_20_output)



chr20 <- data.table::rbindlist(chr20_list, use.names = FALSE)
colnames(chr20) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr20$BayesDel_addAF_score <- as.numeric(chr20$BayesDel_addAF_score)

summary(chr20)
# BayesDel has 1 NA value

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr20)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr20$genename)             
# genenames that are present that shouldn't be!    
#none
# will not need to write code to extract the rows with these genename values!                                                                               

#chr20 <- chr20[!chr20$genename %in%  ]


# chr20 is done (i think) :)

write.csv(chr20, "chr20.csv")












## chromosome 21
chr21_genes <- list.files(pattern = "*_21_output.csv")
chr21_genes <- tools::file_path_sans_ext(chr21_genes)
chr21_genes
chr21_list <- list(end_CFAP410_21_output, end_DONSON_21_output, end_DYRK1A_21_output, 
                   end_PCNT_21_output, end_WDR4_21_output)



chr21 <- data.table::rbindlist(chr21_list, use.names = FALSE)
colnames(chr21) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr21$BayesDel_addAF_score <- as.numeric(chr21$BayesDel_addAF_score)

summary(chr21)
# BayesDel has 9 NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr21)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr21$genename)             
# genenames that are present that shouldn't be!    
#none    
# will not need to write code to extract the rows with these genename values!                                                                               

#chr21 <- chr21[!chr21$genename %in%  ]


# chr21 is done (i think) :)

write.csv(chr21, "chr21.csv")













## chromosome 22
chr22_genes <- list.files(pattern = "*_22_output.csv")
chr22_genes <- tools::file_path_sans_ext(chr22_genes)
chr22_genes
chr22_list <- list(end_CDC45_22_output, end_LZTR1_22_output, end_MCM5_22_output, 
end_PISD_22_output, end_RRP7A_22_output, end_SBF1_22_output, end_TUBGCP6_22_output, 
end_EP300_X_output)
# made a mistake initially naming EP300 - make sure not included in X chr group


chr22 <- data.table::rbindlist(chr22_list, use.names = FALSE)
colnames(chr22) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chr22$BayesDel_addAF_score <- as.numeric(chr22$BayesDel_addAF_score)

summary(chr22)
# BayesDel has no NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chr22)
# far too many variants to go through and manually remove duplicates in diff transcripts



#to check that the correct genes are represented!
unique(chr22$genename)             
# genenames that are present that shouldn't be!                               
#none
# will not need to write code to extract the rows with these genename values!                                                                               

#chr22 <- chr22[!chr22$genename %in%  ]


# chr22 is done (i think) :)

write.csv(chr22, "chr22.csv")











## chromosome X
chrX_genes <- list.files(pattern = "*_X_output.csv")
chrX_genes <- tools::file_path_sans_ext(chrX_genes)
chrX_genes
chrX_list <- list(end_AMMERC1_X_output, end_ARSL_X_output, end_BTK_X_output, 
end_CASK_X_output, end_EBP_X_output, end_EIF2S3_X_output, end_FGD1_X_output, 
end_HDAC8_X_output, end_HMGB3_X_output, end_KDM6A_X_output, end_LAGE3_X_output, 
end_PDHA1_X_output, end_PHEX_X_output, end_PQBP1_X_output, end_SHOX_X_output)
# made a mistake initially naming EP300 - make sure not included in X chr group


chrX <- data.table::rbindlist(chrX_list, use.names = FALSE)
colnames(chrX) <- c("chr", "pos", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "BayesDel_addAF_score", "CADD_phred")
chrX$BayesDel_addAF_score <- as.numeric(chrX$BayesDel_addAF_score)

summary(chrX)
# BayesDel has 54 NA values

# there are some duplicate variants that have different aapos values presumably dependent on transcript
# used the 'unique()' function to check for duplicate rows (there are none)
unique(chrX)
# far too many variants to go through and manually remove duplicates in diff transcripts
unique(chrX$chr)
# just X, as it should be :)


#to check that the correct genes are represented!
unique(chrX$genename)             
# genenames that are present that shouldn't be!                                                                                                                             
# [17] "GPR34;GPR34;GPR34"                                                                                                                                                     
# [18] "GPR34;GPR34;GPR34;GPR34"                                                                                                                                               
# [19] "GPR82"  
# [58] "MAP3K15"   
# will need to write code to extract the rows with these genename values!                                                                               

##### ARX was left out from list of genes (all of them, completely not there)

chrX <- chrX[!chrX$genename %in% c("GPR34;GPR34;GPR34", "GPR34;GPR34;GPR34;GPR34", 
                                   "GPR82", "MAP3K15")  ]


# chrX is done (i think) :)

write.csv(chrX, "chrX.csv")




## have checked that all the genes are present :)













# filtering ---------------------------------------------------------------
# 14/11/23

# load each chr back in, don't want to re-execute the above script.

library(readr)

# using the more conservative model:
# model_2 <- subset(Variant_data_copy, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)




# chr1
chr1 <- read_csv("chr1.csv")
chr1[1] <- NULL               
chr1_path <- subset(chr1, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (6548 / 21219)*100 = 30.85914%
# using the less conservative method:
# (12734/21219)*100 = 60.01225%
# subsetting out the variants that have BayesDel NA values
chr1_CADD_only <- chr1[is.na(chr1$BayesDel_addAF_score),]
chr1_CADD_only_path <- subset(chr1_CADD_only, CADD_phred >= 24.6818)
# 0 obs, so none to add back into the chr1_path df
write.csv(chr1_path, "chr1_path.csv")



# chr2
chr2 <- read_csv("chr2.csv")
chr2[1] <- NULL               
chr2_path <- subset(chr2, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (4735/13747)*100 = 34.44388%
# subsetting out the variants that have BayesDel NA values
chr2_CADD_only <- chr2[is.na(chr2$BayesDel_addAF_score),]
chr2_CADD_only_path <- subset(chr2_CADD_only, CADD_phred >= 24.6818)
# 1 obs, need to add this to the chr2_path df
chr2_path_2 <- rbind(chr2_path, chr2_CADD_only_path)
write.csv(chr2_path_2, "chr2_path_2.csv")


# chr3
chr3 <- read_csv("chr3.csv")
chr3[1] <- NULL               
chr3_path <- subset(chr3, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (4302/10719)*100 = 40.13434%
# subsetting out the variants that have BayesDel NA values
chr3_CADD_only <- chr3[is.na(chr3$BayesDel_addAF_score),]
# 0 obs, don't need to change df
write.csv(chr3_path, "chr3_path.csv")


# chr4
chr4 <- read_csv("chr4.csv")
chr4[1] <- NULL               
chr4_path <- subset(chr4, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (2946/10072)*100 = 29.2494%
# subsetting out the variants that have BayesDel NA values
chr4_CADD_only <- chr4[is.na(chr4$BayesDel_addAF_score),]
# 0 obs, don't need to change df
write.csv(chr4_path, "chr4_path.csv")


# chr5
chr5 <- read_csv("chr5.csv")
chr5[1] <- NULL               
chr5_path <- subset(chr5, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (2530/7825)*100 = 32.33227%
# subsetting out the variants that have BayesDel NA values
chr5_CADD_only <- chr5[is.na(chr5$BayesDel_addAF_score),]
chr5_CADD_only_path <- subset(chr5_CADD_only, CADD_phred >= 24.6818)
# 0 obs, don't need to change df
write.csv(chr5_path, "chr5_path.csv")


# chr6
chr6 <- read_csv("chr6.csv")
chr6[1] <- NULL               
chr6_path <- subset(chr6, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (1854/4880)*100 = 37.9918%
# subsetting out the variants that have BayesDel NA values
chr6_CADD_only <- chr6[is.na(chr6$BayesDel_addAF_score),]
chr6_CADD_only_path <- subset(chr6_CADD_only, CADD_phred >= 24.6818)
# 0 obs, don't need to change df
write.csv(chr6_path, "chr6_path.csv")


# chr7
chr7 <- read_csv("chr7.csv")
chr7[1] <- NULL               
chr7_path <- subset(chr7, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (2401/8329)*100 = 28.82699%
# subsetting out the variants that have BayesDel NA values
chr7_CADD_only <- chr7[is.na(chr7$BayesDel_addAF_score),]
# 0 obs, don't need to change df
write.csv(chr7_path, "chr7_path.csv")


# chr8
chr8 <- read_csv("chr8.csv")
chr8[1] <- NULL               
chr8_path <- subset(chr8, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (2005/7458)*100 = 26.88388%
# subsetting out the variants that have BayesDel NA values
chr8_CADD_only <- chr8[is.na(chr8$BayesDel_addAF_score),]
# 0 obs, don't need to change df
write.csv(chr8_path, "chr8_path.csv")


# chr9
chr9 <- read_csv("chr9.csv")
chr9[1] <- NULL               
chr9_path <- subset(chr9, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (2044/6222)*100 = 32.85117%
# subsetting out the variants that have BayesDel NA values
chr9_CADD_only <- chr9[is.na(chr9$BayesDel_addAF_score),]
# 0 obs, don't need to change df
write.csv(chr9_path, "chr9_path.csv")


# chr10
chr10 <- read_csv("chr10.csv")
chr10[1] <- NULL               
chr10_path <- subset(chr10, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (2061/5785)*100 = 35.62662%
# subsetting out the variants that have BayesDel NA values
chr10_CADD_only <- chr10[is.na(chr10$BayesDel_addAF_score),]
chr10_CADD_only_path <- subset(chr10_CADD_only, CADD_phred >= 24.6818)
# 0 obs, don't need to change df
write.csv(chr10_path, "chr10_path.csv")


# chr11
chr11 <- read_csv("chr11.csv")
chr11[1] <- NULL               
chr11_path <- subset(chr11, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (3555/9926)*100 = 35.81503%
# subsetting out the variants that have BayesDel NA values
chr11_CADD_only <- chr11[is.na(chr11$BayesDel_addAF_score),]
chr11_CADD_only_path <- subset(chr11_CADD_only, CADD_phred >= 24.6818)
# 0 obs, don't need to alter df
write.csv(chr11_path, "chr11_path.csv")


# chr12
chr12 <- read_csv("chr12.csv")
chr12[1] <- NULL               
chr12_path <- subset(chr12, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (3731/10607)*100 =  35.17488%
# subsetting out the variants that have BayesDel NA values
chr12_CADD_only <- chr12[is.na(chr12$BayesDel_addAF_score),]
# 0 obs, don't need to alter df
write.csv(chr12_path, "chr12_path.csv")


# chr13
chr13 <- read_csv("chr13.csv")
chr13[1] <- NULL               
chr13_path <- subset(chr13, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (685/2297)*100 = 29.82151%
# subsetting out the variants that have BayesDel NA values
chr13_CADD_only <- chr13[is.na(chr13$BayesDel_addAF_score),]
# 0 obs, don't need to change df
write.csv(chr13_path, "chr13_path.csv")


# chr14
chr14 <- read_csv("chr14.csv")
chr14[1] <- NULL               
chr14_path <- subset(chr14, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (1149/3946)*100 = 29.11809%
# subsetting out the variants that have BayesDel NA values
chr14_CADD_only <- chr14[is.na(chr14$BayesDel_addAF_score),]
# 0 obs, don't need to change df
write.csv(chr14_path, "chr14_path.csv")


# chr15
chr15 <- read_csv("chr15.csv")
chr15[1] <- NULL               
chr15_path <- subset(chr15, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (1991/8606)*100 = 23.13502%
# subsetting out the variants that have BayesDel NA values
chr15_CADD_only <- chr15[is.na(chr15$BayesDel_addAF_score),]
# 0 obs, don't need to change df
write.csv(chr15_path, "chr15_path.csv")


# chr16
chr16 <- read_csv("chr16.csv")
chr16[1] <- NULL               
chr16_path <- subset(chr16, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (3019/10086)*100 = 29.93258%
# subsetting out the variants that have BayesDel NA values
chr16_CADD_only <- chr16[is.na(chr16$BayesDel_addAF_score),]
chr16_CADD_only_path <- subset(chr16_CADD_only, CADD_phred >= 24.6818)
# 0 obs, don't need to change df
write.csv(chr16_path, "chr16_path.csv")


# chr17
chr17 <- read_csv("chr17.csv")
chr17[1] <- NULL               
chr17_path <- subset(chr17, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (2411/6847)*100 = 35.2125%
# subsetting out the variants that have BayesDel NA values
chr17_CADD_only <- chr17[is.na(chr17$BayesDel_addAF_score),]
chr17_CADD_only_path <- subset(chr17_CADD_only, CADD_phred >= 24.6818)
# 0 obs, don't need to alter df
write.csv(chr17_path, "chr17_path.csv")


# chr18
chr18 <- read_csv("chr18.csv")
chr18[1] <- NULL               
chr18_path <- subset(chr18, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (877/2713)*100 = 32.32584%
# subsetting out the variants that have BayesDel NA values
chr18_CADD_only <- chr18[is.na(chr18$BayesDel_addAF_score),]
# 0 obs, don't need to change df
write.csv(chr18_path, "chr18_path.csv")


# chr19
chr19 <- read_csv("chr19.csv")
chr19[1] <- NULL               
chr19_path <- subset(chr19, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (2023/6304)*100 = 32.09074%
# subsetting out the variants that have BayesDel NA values
chr19_CADD_only <- chr19[is.na(chr19$BayesDel_addAF_score),]
chr19_CADD_only_path <- subset(chr19_CADD_only, CADD_phred >= 24.6818)
# 0 obs, don't need to change df
write.csv(chr19_path, "chr19_path.csv")


# chr20
chr20 <- read_csv("chr20.csv")
chr20[1] <- NULL               
chr20_path <- subset(chr20, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (1134/4080)*100 = 27.79412%
# subsetting out the variants that have BayesDel NA values
chr20_CADD_only <- chr20[is.na(chr20$BayesDel_addAF_score),]
chr20_CADD_only_path <- subset(chr20_CADD_only, CADD_phred >= 24.6818)
# 0 obs, don't need to change df
write.csv(chr20_path, "chr20_path.csv")


# chr21
chr21 <- read_csv("chr21.csv")
chr21[1] <- NULL               
chr21_path <- subset(chr21, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (737/4459)*100 = 16.52837%
# subsetting out the variants that have BayesDel NA values
chr21_CADD_only <- chr21[is.na(chr21$BayesDel_addAF_score),]
chr21_CADD_only_path <- subset(chr21_CADD_only, CADD_phred >= 24.6818)
# 0 obs, don't need to change df
write.csv(chr21_path, "chr21_path.csv")


# chr22
chr22 <- read_csv("chr22.csv")
chr22[1] <- NULL               
chr22_path <- subset(chr22, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (2566/7502)*100 = 34.20421%
# subsetting out the variants that have BayesDel NA values
chr22_CADD_only <- chr22[is.na(chr22$BayesDel_addAF_score),]
# 0 obs, don't need to change df
write.csv(chr22_path, "chr22_path.csv")


# chrX
chrX <- read_csv("chrX.csv")
chrX[1] <- NULL               
chrX_path <- subset(chrX, BayesDel_addAF_score >= -0.1002 & CADD_phred >= 24.6818)
# (584/2930)*100 = 19.93174%
# subsetting out the variants that have BayesDel NA values
chrX_CADD_only <- chrX[is.na(chrX$BayesDel_addAF_score),]
chrX_CADD_only_path <- subset(chrX_CADD_only, CADD_phred >= 24.6818)
# 0 obs, don't need to change df
write.csv(chrX_path, "chrX_path.csv")

## these dataframes can be saved as CSV and directly uploaded into UKB for participant ID extraction

# included:
6548 + 4736 + 4302 + 2946 + 2530 + 1854 + 2401 + 2005 + 2044 + 2061 + 3555 + 3731 + 685 + 1149 + 1991 + 3019 + 2411 + 877 + 2023 + 1134 + 737 + 2566 + 584
# 55889

# total:
21219 + 13747 + 10719 + 10072 + 7825 + 4880 + 8329 + 7458 + 6222 + 5785 + 9926 + 10607 + 2297 + 3946 + 8606 + 10086 + 6847 + 2713 + 6304 + 4080 + 4459 + 7502 + 2930
# 176559

# percentage
55889/176559
# 0.3165457 = 31.65457 %






# just for interst
# combining the variants that have CADD score but not BayesDel 

C_but_no_B_1 <- rbind(chr1_CADD_only, chr2_CADD_only)
C_but_no_B_2 <- rbind(C_but_no_B_1, chr5_CADD_only)
C_but_no_B_3 <- rbind(C_but_no_B_2, chr6_CADD_only)
C_but_no_B_4 <- rbind(C_but_no_B_3, chr10_CADD_only)
C_but_no_B_5 <- rbind(C_but_no_B_4, chr11_CADD_only)
C_but_no_B_6 <- rbind(C_but_no_B_5, chr16_CADD_only)
C_but_no_B_7 <- rbind(C_but_no_B_6, chr17_CADD_only)
C_but_no_B_8 <- rbind(C_but_no_B_7, chr19_CADD_only)
C_but_no_B_9 <- rbind(C_but_no_B_8, chr20_CADD_only)
C_but_no_B_10 <- rbind(C_but_no_B_9, chr21_CADD_only)
C_but_no_B_11 <- rbind(C_but_no_B_10, chrX_CADD_only)

summary(C_but_no_B_11)
# CADD_phred    
# Min.   : 0.013  
# 1st Qu.: 4.058  
# Median : 9.211  
# Mean   :10.506  
# 3rd Qu.:17.200  
# Max.   :25.700 
















