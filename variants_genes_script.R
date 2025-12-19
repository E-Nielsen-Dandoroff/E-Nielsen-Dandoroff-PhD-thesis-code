## 14/02/2024

# code for making gene-variant association for tagging

library(readr)




# CHR 1 -------------------------------------------------------------------


# load in the dataframe
chr1_path <- read_csv("path-scored-variants/chr1_path.csv")

# combine chr:pos:ref:alt_alt

chr1_path$variant <- paste(chr1_path$chr,":",chr1_path$pos,":",chr1_path$ref,":",chr1_path$alt,"_",chr1_path$alt, sep = "")

# just need the variant info and gene names
chr1_genes <- chr1_path[c(9,12)]

chr1_genes$GENE_NAME <- chr1_genes$genename

chr1_genes$GENE <- grepl("B3GALT6", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "B3GALT6"

chr1_genes$GENE <- grepl("ATAD3A", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "ATAD3A"

chr1_genes$GENE <- grepl("SASS6", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "SASS6"

chr1_genes$GENE <- grepl("COL11A1", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "COL11A1"

chr1_genes$GENE <- grepl("TAF13", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "TAF13"

chr1_genes$GENE <- grepl("AMPD2", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "AMPD2"

chr1_genes$GENE <- grepl("NRAS", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "NRAS"

chr1_genes$GENE <- grepl("TBX15", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "TBX15"

chr1_genes$GENE <- grepl("NOTCH2", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "NOTCH2"

chr1_genes$GENE <- grepl("PHGDH", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "PHGDH"

chr1_genes$GENE <- grepl("RIT1", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "RIT1"

chr1_genes$GENE <- grepl("ARHGEF2", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "ARHGEF2"

chr1_genes$GENE <- grepl("LMNA", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "LMNA"

chr1_genes$GENE <- grepl("DDR2", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "DDR2"

chr1_genes$GENE <- grepl("GORAB", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "GORAB"

chr1_genes$GENE <- grepl("LHX4", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "LHX4"

chr1_genes$GENE <- grepl("TSEN15", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "TSEN15"

chr1_genes$GENE <- grepl("ASPM", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "ASPM"

chr1_genes$GENE <- grepl("KIF14", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "KIF14"

chr1_genes$GENE <- grepl("HHAT", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "HHAT"

chr1_genes$GENE <- grepl("CENPF", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "CENPF"

chr1_genes$GENE <- grepl("HSPG2", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "HSPG2"

chr1_genes$GENE <- grepl("LBR", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "LBR"

chr1_genes$GENE <- grepl("NUP133", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "NUP133"

chr1_genes$GENE <- grepl("GNPAT", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "GNPAT"

chr1_genes$GENE <- grepl("TBCE", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "TBCE"

chr1_genes$GENE <- grepl("YARS", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "YARS"

chr1_genes$GENE <- grepl("YRDC", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "YRDC"

chr1_genes$GENE <- grepl("MACF1", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "MACF1"

chr1_genes$GENE <- grepl("MFSD2A", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "MFSD2A"

chr1_genes$GENE <- grepl("ZMPSTE24", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "ZMPSTE24"

chr1_genes$GENE <- grepl("P3H1", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "P3H1"

chr1_genes$GENE <- grepl("TOE1", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "TOE1"

chr1_genes$GENE <- grepl("STIL", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "STIL"

chr1_genes$GENE <- grepl("ORC1", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "ORC1"

chr1_genes$GENE <- grepl("DHCR24", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "DHCR24"

chr1_genes$GENE <- grepl("SLC35D1", chr1_genes$genename)
chr1_genes$GENE_NAME[chr1_genes$GENE == "TRUE"] <- "SLC35D1"


# check each of the genes
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "B3GALT6"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "ATAD3A"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "SASS6"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "COL11A1"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "TAF13"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "AMPD2"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "NRAS"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "TBX15"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "NOTCH2"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "PHGDH"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "RIT1"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "ARHGEF2"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "LMNA"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "DDR2"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "GORAB"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "LHX4"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "TSEN15"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "ASPM"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "KIF14"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "HHAT"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "CENPF"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "HSPG2"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "LBR"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "NUP133"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "GNPAT"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "TBCE"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "YARS"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "YRDC"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "MACF1"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "MFSD2A"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "ZMPSTE24"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "P3H1"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "TOE1"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "STIL"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "ORC1"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "DHCR24"])
unique(chr1_genes$genename[chr1_genes$GENE_NAME == "SLC35D1"])

unique(chr1_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr1_variants_genes <- chr1_genes[c(2,4)]

# save this for later
#save it to R
write.csv(chr1_variants_genes, "chr1_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr1_variants_genes.csv --path Emily-folder/variants-genes/chr1_variants_genes.csv')









# CHR 2 -------------------------------------------------------------------


# load in the dataframe
chr2_path <- read_csv("path-scored-variants/chr2_path.csv")

# combine chr:pos:ref:alt_alt

chr2_path$variant <- paste(chr2_path$chr,":",chr2_path$pos,":",chr2_path$ref,":",chr2_path$alt,"_",chr2_path$alt, sep = "")

# just need the variant info and gene names
chr2_genes <- chr2_path[c(9,12)]

chr2_genes$GENE_NAME <- chr2_genes$genename


chr2_genes$GENE <- grepl("BUB1", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "BUB1"

chr2_genes$GENE <- grepl("CKAP2L", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "CKAP2L"

chr2_genes$GENE <- grepl("GLI2", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "GLI2"

chr2_genes$GENE <- grepl("RNU4ATAC", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "RNU4ATAC"

chr2_genes$GENE <- grepl("ERCC3", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "ERCC3"

chr2_genes$GENE <- grepl("SMPD4", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "SMPD4"

chr2_genes$GENE <- grepl("ORC4", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "ORC4"

chr2_genes$GENE <- grepl("MYCN", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "MYCN"

chr2_genes$GENE <- grepl("DYNC1I2", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "DYNC1I2"

chr2_genes$GENE <- grepl("AGPS", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "AGPS"

chr2_genes$GENE <- grepl("WDR35", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "WDR35"

chr2_genes$GENE <- grepl("MATN3", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "MATN3"

chr2_genes$GENE <- grepl("FN1", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "FN1"

chr2_genes$GENE <- grepl("SMARCAL1", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "SMARCAL1"

chr2_genes$GENE <- grepl("BCS1L", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "BCS1L"

chr2_genes$GENE <- grepl("NHEJ1", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "NHEJ1"

chr2_genes$GENE <- grepl("OBSL1", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "OBSL1"

chr2_genes$GENE <- grepl("IHH", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "IHH"

chr2_genes$GENE <- grepl("NPPC", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "NPPC"

chr2_genes$GENE <- grepl("DNMT3A", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "DNMT3A"

chr2_genes$GENE <- grepl("SOS1", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "SOS1"

chr2_genes$GENE <- grepl("CRIPT", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "CRIPT"

chr2_genes$GENE <- grepl("SOX11", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "SOX11"

chr2_genes$GENE <- grepl("SLC1A4", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "SLC1A4"

chr2_genes$GENE <- grepl("SPRED2", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "SPRED2"

chr2_genes$GENE <- grepl("EXOC6B", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "EXOC6B"

chr2_genes$GENE <- grepl("STAMBP", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "STAMBP"

chr2_genes$GENE <- grepl("TPRKB", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "TPRKB"

chr2_genes$GENE <- grepl("EIF2AK3", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "EIF2AK3"

chr2_genes$GENE <- grepl("NCAPH", chr2_genes$genename)
chr2_genes$GENE_NAME[chr2_genes$GENE == "TRUE"] <- "NCAPH"


# check each of the genes
unique(chr2_genes$genename[chr2_genes$GENE_NAME == "BUB1"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "CKAP2L"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "GLI2"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "RNU4ATAC"])
# none
unique(chr2_genes$genename[chr2_genes$GENE_NAME == "ERCC3"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "SMPD4"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "ORC4"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "MYCN"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "DYNC1I2"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "AGPS"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "WDR35"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "MATN3"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "FN1"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "SMARCAL1"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "BCS1L"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "NHEJ1"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "OBSL1"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "IHH"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "NPPC"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "DNMT3A"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "SOS1"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "CRIPT"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "SOX11"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "SLC1A4"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "SPRED2"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "EXOC6B"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "STAMBP"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "TPRKB"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "EIF2AK3"])

unique(chr2_genes$genename[chr2_genes$GENE_NAME == "NCAPH"])


unique(chr2_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr2_variants_genes <- chr2_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr2_variants_genes, "chr2_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr2_variants_genes.csv --path Emily-folder/variants-genes/chr2_variants_genes.csv')









# CHR 3 -------------------------------------------------------------------

# load in the dataframe
chr3_path <- read_csv("path-scored-variants/chr3_path.csv")

# combine chr:pos:ref:alt_alt

chr3_path$variant <- paste(chr3_path$chr,":",chr3_path$pos,":",chr3_path$ref,":",chr3_path$alt,"_",chr3_path$alt, sep = "")

# just need the variant info and gene names
chr3_genes <- chr3_path[c(9,12)]

chr3_genes$GENE_NAME <- chr3_genes$genename




chr3_genes$GENE <- grepl("TBC1D23", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "TBC1D23"

chr3_genes$GENE <- grepl("NEPRO", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "NEPRO"

chr3_genes$GENE <- grepl("TSEN2", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "TSEN2"

chr3_genes$GENE <- grepl("RAF1", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "RAF1"

chr3_genes$GENE <- grepl("IFT122", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "IFT122"

chr3_genes$GENE <- grepl("CEP63", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "CEP63"

chr3_genes$GENE <- grepl("MRAS", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "MRAS"

chr3_genes$GENE <- grepl("COPB2", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "COPB2"

chr3_genes$GENE <- grepl("RASA2", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "RASA2"

chr3_genes$GENE <- grepl("ATR", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "ATR"

chr3_genes$GENE <- grepl("RNF13", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "RNF13"

chr3_genes$GENE <- grepl("GHSR", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "GHSR"

chr3_genes$GENE <- grepl("DVL3", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "DVL3"

chr3_genes$GENE <- grepl("PCYT1A", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "PCYT1A"

chr3_genes$GENE <- grepl("THRB", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "THRB"

chr3_genes$GENE <- grepl("GLB1", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "GLB1"

chr3_genes$GENE <- grepl("CRTAP", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "CRTAP"

chr3_genes$GENE <- grepl("PDCD6IP", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "PDCD6IP"

chr3_genes$GENE <- grepl("PTH1R", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "PTH1R"

chr3_genes$GENE <- grepl("ATRIP", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "ATRIP"

chr3_genes$GENE <- grepl("TRAIP", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "TRAIP"

chr3_genes$GENE <- grepl("POC1A", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "POC1A"

chr3_genes$GENE <- grepl("WNT5A", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "WNT5A"

chr3_genes$GENE <- grepl("HESX1", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "HESX1"

chr3_genes$GENE <- grepl("FLNB", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "FLNB"

chr3_genes$GENE <- grepl("POU1F1", chr3_genes$genename)
chr3_genes$GENE_NAME[chr3_genes$GENE == "TRUE"] <- "POU1F1"



# check each of the genes
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "TBC1D23"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "NEPRO"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "TSEN2"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "RAF1"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "IFT122"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "CEP63"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "MRAS"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "COPB2"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "RASA2"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "ATR"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "RNF13"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "GHSR"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "DVL3"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "PCYT1A"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "THRB"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "GLB1"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "CRTAP"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "PDCD6IP"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "PTH1R"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "ATRIP"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "TRAIP"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "POC1A"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "WNT5A"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "HESX1"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "FLNB"])
unique(chr3_genes$genename[chr3_genes$GENE_NAME == "POU1F1"])


unique(chr3_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr3_variants_genes <- chr3_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr3_variants_genes, "chr3_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr3_variants_genes.csv --path Emily-folder/variants-genes/chr3_variants_genes.csv')







# CHR 4 -------------------------------------------------------------------


# load in the dataframe
chr4_path <- read_csv("path-scored-variants/chr4_path.csv")

# combine chr:pos:ref:alt_alt

chr4_path$variant <- paste(chr4_path$chr,":",chr4_path$pos,":",chr4_path$ref,":",chr4_path$alt,"_",chr4_path$alt, sep = "")

# just need the variant info and gene names
chr4_genes <- chr4_path[c(9,12)]

chr4_genes$GENE_NAME <- chr4_genes$genename



chr4_genes$GENE <- grepl("PPP3CA", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "PPP3CA"

chr4_genes$GENE <- grepl("SLC39A8", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "SLC39A8"

chr4_genes$GENE <- grepl("CENPE", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "CENPE"

chr4_genes$GENE <- grepl("SGMS2", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "SGMS2"

chr4_genes$GENE <- grepl("LARP7", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "LARP7"

chr4_genes$GENE <- grepl("PLK4", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "PLK4"

chr4_genes$GENE <- grepl("NKX3-2", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "NKX3_2"

chr4_genes$GENE <- grepl("RAB33B", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "RAB33B"

chr4_genes$GENE <- grepl("SLC10A7", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "SLC10A7"

chr4_genes$GENE <- grepl("MSMO1", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "MSMO1"

chr4_genes$GENE <- grepl("NEK1", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "NEK1"

chr4_genes$GENE <- grepl("FGFR3", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "FGFR3"

chr4_genes$GENE <- grepl("UFSP2", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "UFSP2"

chr4_genes$GENE <- grepl("SEPSECS", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "SEPSECS"

chr4_genes$GENE <- grepl("WDR19", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "WDR19"

chr4_genes$GENE <- grepl("EVC", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "EVC"

chr4_genes$GENE <- grepl("TMEM165", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "TMEM165"

chr4_genes$GENE <- grepl("CEP135", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "CEP135"

chr4_genes$GENE <- grepl("EVC2", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "EVC2"

chr4_genes$GENE <- grepl("PRKG2", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "PRKG2"

chr4_genes$GENE <- grepl("WDFY3", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "WDFY3"

chr4_genes$GENE <- grepl("BMPR1B", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "BMPR1B"

chr4_genes$GENE <- grepl("TRMT10A", chr4_genes$genename)
chr4_genes$GENE_NAME[chr4_genes$GENE == "TRUE"] <- "TRMT10A"



# check each of the genes
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "PPP3CA"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "SLC39A8"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "CENPE"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "SGMS2"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "LARP7"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "PLK4"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "NKX3_2"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "RAB33B"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "SLC10A7"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "MSMO1"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "NEK1"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "FGFR3"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "UFSP2"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "SEPSECS"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "WDR19"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "EVC2"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "TMEM165"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "CEP135"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "EVC"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "PRKG2"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "WDFY3"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "BMPR1B"])
unique(chr4_genes$genename[chr4_genes$GENE_NAME == "TRMT10A"])


unique(chr4_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr4_variants_genes <- chr4_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr4_variants_genes, "chr4_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr4_variants_genes.csv --path Emily-folder/variants-genes/chr4_variants_genes.csv')










# CHR 5 -------------------------------------------------------------------

library(readr)

# load in the dataframe
chr5_path <- read_csv("path-scored-variants/chr5_path.csv")

# combine chr:pos:ref:alt_alt

chr5_path$variant <- paste(chr5_path$chr,":",chr5_path$pos,":",chr5_path$ref,":",chr5_path$alt,"_",chr5_path$alt, sep = "")

# just need the variant info and gene names
chr5_genes <- chr5_path[c(9,12)]

chr5_genes$GENE_NAME <- chr5_genes$genename



chr5_genes$GENE <- grepl("LMNB1", chr5_genes$genename)
chr5_genes$GENE_NAME[chr5_genes$GENE == "TRUE"] <- "LMNB1"

chr5_genes$GENE <- grepl("RAD50", chr5_genes$genename)
chr5_genes$GENE_NAME[chr5_genes$GENE == "TRUE"] <- "RAD50"

chr5_genes$GENE <- grepl("PCDH12", chr5_genes$genename)
chr5_genes$GENE_NAME[chr5_genes$GENE == "TRUE"] <- "PCDH12"

chr5_genes$GENE <- grepl("SLC26A2", chr5_genes$genename)
chr5_genes$GENE_NAME[chr5_genes$GENE == "TRUE"] <- "SLC26A2"

chr5_genes$GENE <- grepl("THG1L", chr5_genes$genename)
chr5_genes$GENE_NAME[chr5_genes$GENE == "TRUE"] <- "THG1L"

chr5_genes$GENE <- grepl("B4GALT7", chr5_genes$genename)
chr5_genes$GENE_NAME[chr5_genes$GENE == "TRUE"] <- "B4GALT7"

chr5_genes$GENE <- grepl("PROP1", chr5_genes$genename)
chr5_genes$GENE_NAME[chr5_genes$GENE == "TRUE"] <- "PROP1"

chr5_genes$GENE <- grepl("ADAMTS2", chr5_genes$genename)
chr5_genes$GENE_NAME[chr5_genes$GENE == "TRUE"] <- "ADAMTS2"

chr5_genes$GENE <- grepl("NIPBL", chr5_genes$genename)
chr5_genes$GENE_NAME[chr5_genes$GENE == "TRUE"] <- "NIPBL"

chr5_genes$GENE <- grepl("LIFR", chr5_genes$genename)
chr5_genes$GENE_NAME[chr5_genes$GENE == "TRUE"] <- "LIFR"

chr5_genes$GENE <- grepl("GHR", chr5_genes$genename)
chr5_genes$GENE_NAME[chr5_genes$GENE == "TRUE"] <- "GHR"

chr5_genes$GENE <- grepl("ERCC8", chr5_genes$genename)
chr5_genes$GENE_NAME[chr5_genes$GENE == "TRUE"] <- "ERCC8"

chr5_genes$GENE <- grepl("PIK3R1", chr5_genes$genename)
chr5_genes$GENE_NAME[chr5_genes$GENE == "TRUE"] <- "PIK3R1"

chr5_genes$GENE <- grepl("OCLN", chr5_genes$genename)
chr5_genes$GENE_NAME[chr5_genes$GENE == "TRUE"] <- "OCLN"

chr5_genes$GENE <- grepl("ARSB", chr5_genes$genename)
chr5_genes$GENE_NAME[chr5_genes$GENE == "TRUE"] <- "ARSB"

chr5_genes$GENE <- grepl("XRCC4", chr5_genes$genename)
chr5_genes$GENE_NAME[chr5_genes$GENE == "TRUE"] <- "XRCC4"






# check each of the genes
unique(chr5_genes$genename[chr5_genes$GENE_NAME == "LMNB1"])
unique(chr5_genes$genename[chr5_genes$GENE_NAME == "RAD50"])
unique(chr5_genes$genename[chr5_genes$GENE_NAME == "PCDH12"])
unique(chr5_genes$genename[chr5_genes$GENE_NAME == "SLC26A2"])
unique(chr5_genes$genename[chr5_genes$GENE_NAME == "THG1L"])
unique(chr5_genes$genename[chr5_genes$GENE_NAME == "B4GALT7"])
unique(chr5_genes$genename[chr5_genes$GENE_NAME == "PROP1"])
unique(chr5_genes$genename[chr5_genes$GENE_NAME == "ADAMTS2"])
unique(chr5_genes$genename[chr5_genes$GENE_NAME == "NIPBL"])
unique(chr5_genes$genename[chr5_genes$GENE_NAME == "LIFR"])
unique(chr5_genes$genename[chr5_genes$GENE_NAME == "GHR"])
unique(chr5_genes$genename[chr5_genes$GENE_NAME == "ERCC8"])
unique(chr5_genes$genename[chr5_genes$GENE_NAME == "PIK3R1"])
unique(chr5_genes$genename[chr5_genes$GENE_NAME == "OCLN"])
unique(chr5_genes$genename[chr5_genes$GENE_NAME == "ARSB"])
unique(chr5_genes$genename[chr5_genes$GENE_NAME == "XRCC4"])

unique(chr5_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr5_variants_genes <- chr5_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr5_variants_genes, "chr5_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr5_variants_genes.csv --path Emily-folder/variants-genes/chr5_variants_genes.csv')








# CHR 6 -------------------------------------------------------------------


# load in the dataframe
chr6_path <- read_csv("path-scored-variants/chr6_path.csv")

# combine chr:pos:ref:alt_alt

chr6_path$variant <- paste(chr6_path$chr,":",chr6_path$pos,":",chr6_path$ref,":",chr6_path$alt,"_",chr6_path$alt, sep = "")

# just need the variant info and gene names
chr6_genes <- chr6_path[c(9,12)]

chr6_genes$GENE_NAME <- chr6_genes$genename



chr6_genes$GENE <- grepl("CDC40", chr6_genes$genename)
chr6_genes$GENE_NAME[chr6_genes$GENE == "TRUE"] <- "CDC40"

chr6_genes$GENE <- grepl("WISP3", chr6_genes$genename)
chr6_genes$GENE_NAME[chr6_genes$GENE == "TRUE"] <- "CCN6"

chr6_genes$GENE <- grepl("COL10A1", chr6_genes$genename)
chr6_genes$GENE_NAME[chr6_genes$GENE == "TRUE"] <- "COL10A1"

chr6_genes$GENE <- grepl("PEX7", chr6_genes$genename)
chr6_genes$GENE_NAME[chr6_genes$GENE == "TRUE"] <- "PEX7"

chr6_genes$GENE <- grepl("GMNN", chr6_genes$genename)
chr6_genes$GENE_NAME[chr6_genes$GENE == "TRUE"] <- "GMNN"

chr6_genes$GENE <- grepl("TUBB2B", chr6_genes$genename)
chr6_genes$GENE_NAME[chr6_genes$GENE == "TRUE"] <- "TUBB2B"

chr6_genes$GENE <- grepl("COL11A2", chr6_genes$genename)
chr6_genes$GENE_NAME[chr6_genes$GENE == "TRUE"] <- "COL11A2"

chr6_genes$GENE <- grepl("PPIL1", chr6_genes$genename)
chr6_genes$GENE_NAME[chr6_genes$GENE == "TRUE"] <- "PPIL1"

chr6_genes$GENE <- grepl("CUL7", chr6_genes$genename)
chr6_genes$GENE_NAME[chr6_genes$GENE == "TRUE"] <- "CUL7"

chr6_genes$GENE <- grepl("MCM3", chr6_genes$genename)
chr6_genes$GENE_NAME[chr6_genes$GENE == "TRUE"] <- "MCM3"

chr6_genes$GENE <- grepl("ICK", chr6_genes$genename)
chr6_genes$GENE_NAME[chr6_genes$GENE == "TRUE"] <- "ICK"



# check each of the genes
unique(chr6_genes$genename[chr6_genes$GENE_NAME == "CDC40"])
unique(chr6_genes$genename[chr6_genes$GENE_NAME == "CCN6"])
unique(chr6_genes$genename[chr6_genes$GENE_NAME == "COL10A1"])
unique(chr6_genes$genename[chr6_genes$GENE_NAME == "PEX7"])
unique(chr6_genes$genename[chr6_genes$GENE_NAME == "GMNN"])
unique(chr6_genes$genename[chr6_genes$GENE_NAME == "TUBB2B"])
unique(chr6_genes$genename[chr6_genes$GENE_NAME == "COL11A2"])
unique(chr6_genes$genename[chr6_genes$GENE_NAME == "PPIL1"])
unique(chr6_genes$genename[chr6_genes$GENE_NAME == "CUL7"])
unique(chr6_genes$genename[chr6_genes$GENE_NAME == "MCM3"])
unique(chr6_genes$genename[chr6_genes$GENE_NAME == "ICK"])


unique(chr6_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr6_variants_genes <- chr6_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr6_variants_genes, "chr6_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr6_variants_genes.csv --path Emily-folder/variants-genes/chr6_variants_genes.csv')








# CHR 7 -------------------------------------------------------------------

# load in the dataframe
chr7_path <- read_csv("path-scored-variants/chr7_path.csv")

# combine chr:pos:ref:alt_alt

chr7_path$variant <- paste(chr7_path$chr,":",chr7_path$pos,":",chr7_path$ref,":",chr7_path$alt,"_",chr7_path$alt, sep = "")

# just need the variant info and gene names
chr7_genes <- chr7_path[c(9,12)]

chr7_genes$GENE_NAME <- chr7_genes$genename



chr7_genes$GENE <- grepl("MCM7", chr7_genes$genename)
chr7_genes$GENE_NAME[chr7_genes$GENE == "TRUE"] <- "MCM7"

chr7_genes$GENE <- grepl("C7orf43", chr7_genes$genename)
chr7_genes$GENE_NAME[chr7_genes$GENE == "TRUE"] <- "TRAPPC14"

chr7_genes$GENE <- grepl("BRAF", chr7_genes$genename)
chr7_genes$GENE_NAME[chr7_genes$GENE == "TRUE"] <- "BRAF"

chr7_genes$GENE <- grepl("SHH", chr7_genes$genename)
chr7_genes$GENE_NAME[chr7_genes$GENE == "TRUE"] <- "SHH"

chr7_genes$GENE <- grepl("FAM20C", chr7_genes$genename)
chr7_genes$GENE_NAME[chr7_genes$GENE == "TRUE"] <- "FAM20C"

chr7_genes$GENE <- grepl("GHRHR", chr7_genes$genename)
chr7_genes$GENE_NAME[chr7_genes$GENE == "TRUE"] <- "GHRHR"

chr7_genes$GENE <- grepl("RALA", chr7_genes$genename)
chr7_genes$GENE_NAME[chr7_genes$GENE == "TRUE"] <- "RALA"

chr7_genes$GENE <- grepl("GLI3", chr7_genes$genename)
chr7_genes$GENE_NAME[chr7_genes$GENE == "TRUE"] <- "GLI3"

chr7_genes$GENE <- grepl("ACTB", chr7_genes$genename)
chr7_genes$GENE_NAME[chr7_genes$GENE == "TRUE"] <- "ACTB"

chr7_genes$GENE <- grepl("KDELR2", chr7_genes$genename)
chr7_genes$GENE_NAME[chr7_genes$GENE == "TRUE"] <- "KDELR2"

chr7_genes$GENE <- grepl("GUSB", chr7_genes$genename)
chr7_genes$GENE_NAME[chr7_genes$GENE == "TRUE"] <- "GUSB"

chr7_genes$GENE <- grepl("SBDS", chr7_genes$genename)
chr7_genes$GENE_NAME[chr7_genes$GENE == "TRUE"] <- "SBDS"

chr7_genes$GENE <- grepl("PCLO", chr7_genes$genename)
chr7_genes$GENE_NAME[chr7_genes$GENE == "TRUE"] <- "PCLO"

chr7_genes$GENE <- grepl("CDK6", chr7_genes$genename)
chr7_genes$GENE_NAME[chr7_genes$GENE == "TRUE"] <- "CDK6"

chr7_genes$GENE <- grepl("COL1A2", chr7_genes$genename)
chr7_genes$GENE_NAME[chr7_genes$GENE == "TRUE"] <- "COL1A2"

chr7_genes$GENE <- grepl("ASNS", chr7_genes$genename)
chr7_genes$GENE_NAME[chr7_genes$GENE == "TRUE"] <- "ASNS"







# check each of the genes
unique(chr7_genes$genename[chr7_genes$GENE_NAME == "MCM7"])
unique(chr7_genes$genename[chr7_genes$GENE_NAME == "TRAPPC14"])
unique(chr7_genes$genename[chr7_genes$GENE_NAME == "BRAF"])
unique(chr7_genes$genename[chr7_genes$GENE_NAME == "SHH"])
unique(chr7_genes$genename[chr7_genes$GENE_NAME == "FAM20C"])
unique(chr7_genes$genename[chr7_genes$GENE_NAME == "GHRHR"])
unique(chr7_genes$genename[chr7_genes$GENE_NAME == "RALA"])
unique(chr7_genes$genename[chr7_genes$GENE_NAME == "GLI3"])
unique(chr7_genes$genename[chr7_genes$GENE_NAME == "ACTB"])
unique(chr7_genes$genename[chr7_genes$GENE_NAME == "KDELR2"])
unique(chr7_genes$genename[chr7_genes$GENE_NAME == "GUSB"])
unique(chr7_genes$genename[chr7_genes$GENE_NAME == "SBDS"])
unique(chr7_genes$genename[chr7_genes$GENE_NAME == "PCLO"])
unique(chr7_genes$genename[chr7_genes$GENE_NAME == "CDK6"])
unique(chr7_genes$genename[chr7_genes$GENE_NAME == "COL1A2"])
unique(chr7_genes$genename[chr7_genes$GENE_NAME == "ASNS"])


unique(chr7_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr7_variants_genes <- chr7_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr7_variants_genes, "chr7_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr7_variants_genes.csv --path Emily-folder/variants-genes/chr7_variants_genes.csv')








# CHR 8 -------------------------------------------------------------------


# load in the dataframe
chr8_path <- read_csv("path-scored-variants/chr8_path.csv")

# combine chr:pos:ref:alt_alt

chr8_path$variant <- paste(chr8_path$chr,":",chr8_path$pos,":",chr8_path$ref,":",chr8_path$alt,"_",chr8_path$alt, sep = "")

# just need the variant info and gene names
chr8_genes <- chr8_path[c(9,12)]

chr8_genes$GENE_NAME <- chr8_genes$genename



chr8_genes$GENE <- grepl("TRPS1", chr8_genes$genename)
chr8_genes$GENE_NAME[chr8_genes$GENE == "TRUE"] <- "TRPS1"

chr8_genes$GENE <- grepl("RAD21", chr8_genes$genename)
chr8_genes$GENE_NAME[chr8_genes$GENE == "TRUE"] <- "RAD21"

chr8_genes$GENE <- grepl("TNFRSF11B", chr8_genes$genename)
chr8_genes$GENE_NAME[chr8_genes$GENE == "TRUE"] <- "TNFRSF11B"

chr8_genes$GENE <- grepl("NSMCE2", chr8_genes$genename)
chr8_genes$GENE_NAME[chr8_genes$GENE == "TRUE"] <- "NSMCE2"

chr8_genes$GENE <- grepl("PUF60", chr8_genes$genename)
chr8_genes$GENE_NAME[chr8_genes$GENE == "TRUE"] <- "PUF60"

chr8_genes$GENE <- grepl("TONSL", chr8_genes$genename)
chr8_genes$GENE_NAME[chr8_genes$GENE == "TRUE"] <- "TONSL"

chr8_genes$GENE <- grepl("RECQL4", chr8_genes$genename)
chr8_genes$GENE_NAME[chr8_genes$GENE == "TRUE"] <- "RECQL4"

chr8_genes$GENE <- grepl("BMP1", chr8_genes$genename)
chr8_genes$GENE_NAME[chr8_genes$GENE == "TRUE"] <- "BMP1"

chr8_genes$GENE <- grepl("ESCO2", chr8_genes$genename)
chr8_genes$GENE_NAME[chr8_genes$GENE == "TRUE"] <- "ESCO2"

chr8_genes$GENE <- grepl("FGFR1", chr8_genes$genename)
chr8_genes$GENE_NAME[chr8_genes$GENE == "TRUE"] <- "FGFR1"

chr8_genes$GENE <- grepl("IMPAD1", chr8_genes$genename)
chr8_genes$GENE_NAME[chr8_genes$GENE == "TRUE"] <- "BPNT2"

chr8_genes$GENE <- grepl("MCPH1", chr8_genes$genename)
chr8_genes$GENE_NAME[chr8_genes$GENE == "TRUE"] <- "MCPH1"

chr8_genes$GENE <- grepl("PTDSS1", chr8_genes$genename)
chr8_genes$GENE_NAME[chr8_genes$GENE == "TRUE"] <- "PTDSS1"

chr8_genes$GENE <- grepl("POP1", chr8_genes$genename)
chr8_genes$GENE_NAME[chr8_genes$GENE == "TRUE"] <- "POP1"









# check each of the genes
unique(chr8_genes$genename[chr8_genes$GENE_NAME == "TRPS1"])
unique(chr8_genes$genename[chr8_genes$GENE_NAME == "RAD21"])
unique(chr8_genes$genename[chr8_genes$GENE_NAME == "TNFRSF11B"])
unique(chr8_genes$genename[chr8_genes$GENE_NAME == "NSMCE2"])
unique(chr8_genes$genename[chr8_genes$GENE_NAME == "PUF60"])
unique(chr8_genes$genename[chr8_genes$GENE_NAME == "TONSL"])
unique(chr8_genes$genename[chr8_genes$GENE_NAME == "RECQL4"])
unique(chr8_genes$genename[chr8_genes$GENE_NAME == "BMP1"])
unique(chr8_genes$genename[chr8_genes$GENE_NAME == "ESCO2"])
unique(chr8_genes$genename[chr8_genes$GENE_NAME == "FGFR1"])
unique(chr8_genes$genename[chr8_genes$GENE_NAME == "BPNT2"])
unique(chr8_genes$genename[chr8_genes$GENE_NAME == "MCPH1"])
unique(chr8_genes$genename[chr8_genes$GENE_NAME == "PTDSS1"])
unique(chr8_genes$genename[chr8_genes$GENE_NAME == "POP1"])


unique(chr8_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr8_variants_genes <- chr8_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr8_variants_genes, "chr8_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr8_variants_genes.csv --path Emily-folder/variants-genes/chr8_variants_genes.csv')








# CHR 9 -------------------------------------------------------------------


# load in the dataframe
chr9_path <- read_csv("path-scored-variants/chr9_path.csv")

# combine chr:pos:ref:alt_alt

chr9_path$variant <- paste(chr9_path$chr,":",chr9_path$pos,":",chr9_path$ref,":",chr9_path$alt,"_",chr9_path$alt, sep = "")

# just need the variant info and gene names
chr9_genes <- chr9_path[c(9,12)]

chr9_genes$GENE_NAME <- chr9_genes$genename



chr9_genes$GENE <- grepl("TMEM38B", chr9_genes$genename)
chr9_genes$GENE_NAME[chr9_genes$GENE == "TRUE"] <- "TMEM38B"

chr9_genes$GENE <- grepl("COL27A1", chr9_genes$genename)
chr9_genes$GENE_NAME[chr9_genes$GENE == "TRUE"] <- "COL27A1"

chr9_genes$GENE <- grepl("CDK5RAP2", chr9_genes$genename)
chr9_genes$GENE_NAME[chr9_genes$GENE == "TRUE"] <- "CDK5RAP2"

chr9_genes$GENE <- grepl("ADAMTSL2", chr9_genes$genename)
chr9_genes$GENE_NAME[chr9_genes$GENE == "TRUE"] <- "ADAMTSL2"

chr9_genes$GENE <- grepl("LHX3", chr9_genes$genename)
chr9_genes$GENE_NAME[chr9_genes$GENE == "TRUE"] <- "LHX3"

chr9_genes$GENE <- grepl("SMARCA2", chr9_genes$genename)
chr9_genes$GENE_NAME[chr9_genes$GENE == "TRUE"] <- "SMARCA2"

chr9_genes$GENE <- grepl("NPR2", chr9_genes$genename)
chr9_genes$GENE_NAME[chr9_genes$GENE == "TRUE"] <- "NPR2"

chr9_genes$GENE <- grepl("EXOSC3", chr9_genes$genename)
chr9_genes$GENE_NAME[chr9_genes$GENE == "TRUE"] <- "EXOSC3"

chr9_genes$GENE <- grepl("PSAT1", chr9_genes$genename)
chr9_genes$GENE_NAME[chr9_genes$GENE == "TRUE"] <- "PSAT1"

chr9_genes$GENE <- grepl("ROR2", chr9_genes$genename)
chr9_genes$GENE_NAME[chr9_genes$GENE == "TRUE"] <- "ROR2"

chr9_genes$GENE <- grepl("NANS", chr9_genes$genename)
chr9_genes$GENE_NAME[chr9_genes$GENE == "TRUE"] <- "NANS"





# check each of the genes
unique(chr9_genes$genename[chr9_genes$GENE_NAME == "TMEM38B"])
unique(chr9_genes$genename[chr9_genes$GENE_NAME == "COL27A1"])
unique(chr9_genes$genename[chr9_genes$GENE_NAME == "CDK5RAP2"])
unique(chr9_genes$genename[chr9_genes$GENE_NAME == "ADAMTSL2"])
unique(chr9_genes$genename[chr9_genes$GENE_NAME == "LHX3"])
unique(chr9_genes$genename[chr9_genes$GENE_NAME == "SMARCA2"])
unique(chr9_genes$genename[chr9_genes$GENE_NAME == "NPR2"])
unique(chr9_genes$genename[chr9_genes$GENE_NAME == "EXOSC3"])
unique(chr9_genes$genename[chr9_genes$GENE_NAME == "PSAT1"])
unique(chr9_genes$genename[chr9_genes$GENE_NAME == "ROR2"])
unique(chr9_genes$genename[chr9_genes$GENE_NAME == "NANS"])

unique(chr9_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr9_variants_genes <- chr9_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr9_variants_genes, "chr9_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr9_variants_genes.csv --path Emily-folder/variants-genes/chr9_variants_genes.csv')







# CHR 10 ------------------------------------------------------------------

# load in the dataframe
chr10_path <- read_csv("path-scored-variants/chr10_path.csv")

# combine chr:pos:ref:alt_alt

chr10_path$variant <- paste(chr10_path$chr,":",chr10_path$pos,":",chr10_path$ref,":",chr10_path$alt,"_",chr10_path$alt, sep = "")

# just need the variant info and gene names
chr10_genes <- chr10_path[c(9,12)]

chr10_genes$GENE_NAME <- chr10_genes$genename



chr10_genes$GENE <- grepl("SMC3", chr10_genes$genename)
chr10_genes$GENE_NAME[chr10_genes$GENE == "TRUE"] <- "SMC3"

chr10_genes$GENE <- grepl("SHOC2", chr10_genes$genename)
chr10_genes$GENE_NAME[chr10_genes$GENE == "TRUE"] <- "SHOC2"

chr10_genes$GENE <- grepl("FRMD4A", chr10_genes$genename)
chr10_genes$GENE_NAME[chr10_genes$GENE == "TRUE"] <- "FRMD4A"

chr10_genes$GENE <- grepl("ERCC6", chr10_genes$genename)
chr10_genes$GENE_NAME[chr10_genes$GENE == "TRUE"] <- "ERCC6"

chr10_genes$GENE <- grepl("KIF1BP", chr10_genes$genename)
chr10_genes$GENE_NAME[chr10_genes$GENE == "TRUE"] <- "KIFBP"

chr10_genes$GENE <- grepl("DNA2", chr10_genes$genename)
chr10_genes$GENE_NAME[chr10_genes$GENE == "TRUE"] <- "DNA2"

chr10_genes$GENE <- grepl("CHST3", chr10_genes$genename)
chr10_genes$GENE_NAME[chr10_genes$GENE == "TRUE"] <- "CHST3"

chr10_genes$GENE <- grepl("MINPP1", chr10_genes$genename)
chr10_genes$GENE_NAME[chr10_genes$GENE == "TRUE"] <- "MINPP1"

chr10_genes$GENE <- grepl("PAPSS2", chr10_genes$genename)
chr10_genes$GENE_NAME[chr10_genes$GENE == "TRUE"] <- "PAPSS2"

chr10_genes$GENE <- grepl("KIF11", chr10_genes$genename)
chr10_genes$GENE_NAME[chr10_genes$GENE == "TRUE"] <- "KIF11"

chr10_genes$GENE <- grepl("ALDH18A1", chr10_genes$genename)
chr10_genes$GENE_NAME[chr10_genes$GENE == "TRUE"] <- "ALDH18A1"



# check each of the genes
unique(chr10_genes$genename[chr10_genes$GENE_NAME == "SMC3"])
unique(chr10_genes$genename[chr10_genes$GENE_NAME == "SHOC2"])
unique(chr10_genes$genename[chr10_genes$GENE_NAME == "FRMD4A"])
unique(chr10_genes$genename[chr10_genes$GENE_NAME == "ERCC6"])
unique(chr10_genes$genename[chr10_genes$GENE_NAME == "KIFBP"])
unique(chr10_genes$genename[chr10_genes$GENE_NAME == "DNA2"])
unique(chr10_genes$genename[chr10_genes$GENE_NAME == "CHST3"])
unique(chr10_genes$genename[chr10_genes$GENE_NAME == "MINPP1"])
unique(chr10_genes$genename[chr10_genes$GENE_NAME == "PAPSS2"])
unique(chr10_genes$genename[chr10_genes$GENE_NAME == "KIF11"])
unique(chr10_genes$genename[chr10_genes$GENE_NAME == "ALDH18A1"])


unique(chr10_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr10_variants_genes <- chr10_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr10_variants_genes, "chr10_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr10_variants_genes.csv --path Emily-folder/variants-genes/chr10_variants_genes.csv')








# CHR 11 ------------------------------------------------------------------

# load in the dataframe
chr11_path <- read_csv("path-scored-variants/chr11_path.csv")

# combine chr:pos:ref:alt_alt

chr11_path$variant <- paste(chr11_path$chr,":",chr11_path$pos,":",chr11_path$ref,":",chr11_path$alt,"_",chr11_path$alt, sep = "")

# just need the variant info and gene names
chr11_genes <- chr11_path[c(9,12)]

chr11_genes$GENE_NAME <- chr11_genes$genename



chr11_genes$GENE <- grepl("MMP13", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "MMP13"

chr11_genes$GENE <- grepl("DYNC2H1", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "DYNC2H1"

chr11_genes$GENE <- grepl("SIK3", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "SIK3"

chr11_genes$GENE <- grepl("ARCN1", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "ARCN1"

chr11_genes$GENE <- grepl("CBL", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "CBL"

chr11_genes$GENE <- grepl("NCAPD3", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "NCAPD3"

chr11_genes$GENE <- grepl("RRAS2", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "RRAS2"

chr11_genes$GENE <- grepl("CDKN1C", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "CDKN1C"

chr11_genes$GENE <- grepl("IFITM5", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "IFITM5"

chr11_genes$GENE <- grepl("SLC35C1", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "SLC35C1"

chr11_genes$GENE <- grepl("CREB3L1", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "CREB3L1"

chr11_genes$GENE <- grepl("SLC39A13", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "SLC39A13"

chr11_genes$GENE <- grepl("CLP1", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "CLP1"

chr11_genes$GENE <- grepl("FAM111A", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "FAM111A"

chr11_genes$GENE <- grepl("B3GAT3", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "B3GAT3"

chr11_genes$GENE <- grepl("LTBP3", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "LTBP3"

chr11_genes$GENE <- grepl("DHCR7", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "DHCR7"

chr11_genes$GENE <- grepl("INPPL1", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "INPPL1"

chr11_genes$GENE <- grepl("TALDO1", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "TALDO1"

chr11_genes$GENE <- grepl("SERPINH1", chr11_genes$genename)
chr11_genes$GENE_NAME[chr11_genes$GENE == "TRUE"] <- "SERPINH1"





# check each of the genes
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "MMP13"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "DYNC2H1"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "SIK3"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "ARCN1"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "CBL"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "NCAPD3"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "RRAS2"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "CDKN1C"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "IFITM5"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "SLC35C1"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "CREB3L1"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "SLC39A13"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "CLP1"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "FAM111A"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "B3GAT3"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "LTBP3"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "DHCR7"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "INPPL1"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "TALDO1"])
unique(chr11_genes$genename[chr11_genes$GENE_NAME == "SERPINH1"])


unique(chr11_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr11_variants_genes <- chr11_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr11_variants_genes, "chr11_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr11_variants_genes.csv --path Emily-folder/variants-genes/chr11_variants_genes.csv')








# CHR 12 ------------------------------------------------------------------

# load in the dataframe
chr12_path <- read_csv("path-scored-variants/chr12_path.csv")

# combine chr:pos:ref:alt_alt

chr12_path$variant <- paste(chr12_path$chr,":",chr12_path$pos,":",chr12_path$ref,":",chr12_path$alt,"_",chr12_path$alt, sep = "")

# just need the variant info and gene names
chr12_genes <- chr12_path[c(9,12)]

chr12_genes$GENE_NAME <- chr12_genes$genename



chr12_genes$GENE <- grepl("GNPTAB", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "GNPTAB"

chr12_genes$GENE <- grepl("NUP37", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "NUP37"

chr12_genes$GENE <- grepl("IGF1", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "IGF1"

chr12_genes$GENE <- grepl("TRPV4", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "TRPV4"

chr12_genes$GENE <- grepl("PTPN11", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "PTPN11"

chr12_genes$GENE <- grepl("TBX3", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "TBX3"

chr12_genes$GENE <- grepl("CIT", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "CIT"

chr12_genes$GENE <- grepl("POLE", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "POLE"

chr12_genes$GENE <- grepl("ANKLE2", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "ANKLE2"

chr12_genes$GENE <- grepl("KRAS", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "KRAS"

chr12_genes$GENE <- grepl("PPFIBP1", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "PPFIBP1"

chr12_genes$GENE <- grepl("PTHLH", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "PTHLH"

chr12_genes$GENE <- grepl("COL2A1", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "COL2A1"

chr12_genes$GENE <- grepl("WNT1", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "WNT1"

chr12_genes$GENE <- grepl("TUBA1A", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "TUBA1A"

chr12_genes$GENE <- grepl("SP7", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "SP7"

chr12_genes$GENE <- grepl("PRIM1", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "PRIM1"

chr12_genes$GENE <- grepl("GLI1", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "GLI1"

chr12_genes$GENE <- grepl("NCAPD2", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "NCAPD2"

chr12_genes$GENE <- grepl("LEMD3", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "LEMD3"

chr12_genes$GENE <- grepl("HMGA2", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "HMGA2"

chr12_genes$GENE <- grepl("PEX5", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "PEX5"

chr12_genes$GENE <- grepl("PHC1", chr12_genes$genename)
chr12_genes$GENE_NAME[chr12_genes$GENE == "TRUE"] <- "PHC1"







# check each of the genes
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "GNPTAB"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "NUP37"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "IGF1"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "TRPV4"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "PTPN11"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "TBX3"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "CIT"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "POLE"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "ANKLE2"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "KRAS"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "PPFIBP1"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "PTHLH"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "COL2A1"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "WNT1"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "TUBA1A"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "SP7"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "PRIM1"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "GLI1"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "NCAPD2"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "LEMD3"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "HMGA2"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "PEX5"])
unique(chr12_genes$genename[chr12_genes$GENE_NAME == "PHC1"])

unique(chr12_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr12_variants_genes <- chr12_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr12_variants_genes, "chr12_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr12_variants_genes.csv --path Emily-folder/variants-genes/chr12_variants_genes.csv')







# CHR 13 ------------------------------------------------------------------

# load in the dataframe
chr13_path <- read_csv("path-scored-variants/chr13_path.csv")

# combine chr:pos:ref:alt_alt

chr13_path$variant <- paste(chr13_path$chr,":",chr13_path$pos,":",chr13_path$ref,":",chr13_path$alt,"_",chr13_path$alt, sep = "")

# just need the variant info and gene names
chr13_genes <- chr13_path[c(9,12)]

chr13_genes$GENE_NAME <- chr13_genes$genename



chr13_genes$GENE <- grepl("LIG4", chr13_genes$genename)
chr13_genes$GENE_NAME[chr13_genes$GENE == "TRUE"] <- "LIG4"

chr13_genes$GENE <- grepl("CENPJ", chr13_genes$genename)
chr13_genes$GENE_NAME[chr13_genes$GENE == "TRUE"] <- "CENPJ"

chr13_genes$GENE <- grepl("B3GLCT", chr13_genes$genename)
chr13_genes$GENE_NAME[chr13_genes$GENE == "TRUE"] <- "B3GLCT"

chr13_genes$GENE <- grepl("GPC6", chr13_genes$genename)
chr13_genes$GENE_NAME[chr13_genes$GENE == "TRUE"] <- "GPC6"




# check each of the genes
unique(chr13_genes$genename[chr13_genes$GENE_NAME == "LIG4"])
unique(chr13_genes$genename[chr13_genes$GENE_NAME == "CENPJ"])
unique(chr13_genes$genename[chr13_genes$GENE_NAME == "B3GLCT"])
unique(chr13_genes$genename[chr13_genes$GENE_NAME == "GPC6"])


unique(chr13_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr13_variants_genes <- chr13_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr13_variants_genes, "chr13_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr13_variants_genes.csv --path Emily-folder/variants-genes/chr13_variants_genes.csv')








# CHR 14 ------------------------------------------------------------------

# load in the dataframe
chr14_path <- read_csv("path-scored-variants/chr14_path.csv")

# combine chr:pos:ref:alt_alt

chr14_path$variant <- paste(chr14_path$chr,":",chr14_path$pos,":",chr14_path$ref,":",chr14_path$alt,"_",chr14_path$alt, sep = "")

# just need the variant info and gene names
chr14_genes <- chr14_path[c(9,12)]

chr14_genes$GENE_NAME <- chr14_genes$genename



chr14_genes$GENE <- grepl("OSGEP", chr14_genes$genename)
chr14_genes$GENE_NAME[chr14_genes$GENE == "TRUE"] <- "OSGEP"

chr14_genes$GENE <- grepl("FOXG1", chr14_genes$genename)
chr14_genes$GENE_NAME[chr14_genes$GENE == "TRUE"] <- "FOXG1"

chr14_genes$GENE <- grepl("TRAPPC6B", chr14_genes$genename)
chr14_genes$GENE_NAME[chr14_genes$GENE == "TRUE"] <- "TRAPPC6B"

chr14_genes$GENE <- grepl("NIN", chr14_genes$genename)
chr14_genes$GENE_NAME[chr14_genes$GENE == "TRUE"] <- "NIN"

chr14_genes$GENE <- grepl("OTX2", chr14_genes$genename)
chr14_genes$GENE_NAME[chr14_genes$GENE == "TRUE"] <- "OTX2"

chr14_genes$GENE <- grepl("TRIP11", chr14_genes$genename)
chr14_genes$GENE_NAME[chr14_genes$GENE == "TRUE"] <- "TRIP11"

chr14_genes$GENE <- grepl("GON7", chr14_genes$genename)
chr14_genes$GENE_NAME[chr14_genes$GENE == "TRUE"] <- "GON7"

chr14_genes$GENE <- grepl("GSC", chr14_genes$genename)
chr14_genes$GENE_NAME[chr14_genes$GENE == "TRUE"] <- "GSC"

chr14_genes$GENE <- grepl("VRK1", chr14_genes$genename)
chr14_genes$GENE_NAME[chr14_genes$GENE == "TRUE"] <- "VRK1"



# check each of the genes
unique(chr14_genes$genename[chr14_genes$GENE_NAME == "OSGEP"])
unique(chr14_genes$genename[chr14_genes$GENE_NAME == "FOXG1"])
unique(chr14_genes$genename[chr14_genes$GENE_NAME == "TRAPPC6B"])
unique(chr14_genes$genename[chr14_genes$GENE_NAME == "NIN"])
unique(chr14_genes$genename[chr14_genes$GENE_NAME == "OTX2"])
unique(chr14_genes$genename[chr14_genes$GENE_NAME == "TRIP11"])
unique(chr14_genes$genename[chr14_genes$GENE_NAME == "GON7"])
unique(chr14_genes$genename[chr14_genes$GENE_NAME == "GSC"])
unique(chr14_genes$genename[chr14_genes$GENE_NAME == "VRK1"])

unique(chr14_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr14_variants_genes <- chr14_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr14_variants_genes, "chr14_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr14_variants_genes.csv --path Emily-folder/variants-genes/chr14_variants_genes.csv')










# CHR 15 ------------------------------------------------------------------


# load in the dataframe
chr15_path <- read_csv("path-scored-variants/chr15_path.csv")

# combine chr:pos:ref:alt_alt

chr15_path$variant <- paste(chr15_path$chr,":",chr15_path$pos,":",chr15_path$ref,":",chr15_path$alt,"_",chr15_path$alt, sep = "")

# just need the variant info and gene names
chr15_genes <- chr15_path[c(9,12)]

chr15_genes$GENE_NAME <- chr15_genes$genename



chr15_genes$GENE <- grepl("BUB1B", chr15_genes$genename)
chr15_genes$GENE_NAME[chr15_genes$GENE == "TRUE"] <- "BUB1B"

chr15_genes$GENE <- grepl("KNL1", chr15_genes$genename)
chr15_genes$GENE_NAME[chr15_genes$GENE == "TRUE"] <- "KNL1"

chr15_genes$GENE <- grepl("UBR1", chr15_genes$genename)
chr15_genes$GENE_NAME[chr15_genes$GENE == "TRUE"] <- "UBR1"

chr15_genes$GENE <- grepl("TUBGCP4", chr15_genes$genename)
chr15_genes$GENE_NAME[chr15_genes$GENE == "TRUE"] <- "TUBGCP4"

chr15_genes$GENE <- grepl("CEP152", chr15_genes$genename)
chr15_genes$GENE_NAME[chr15_genes$GENE == "TRUE"] <- "CEP152"

chr15_genes$GENE <- grepl("FBN1", chr15_genes$genename)
chr15_genes$GENE_NAME[chr15_genes$GENE == "TRUE"] <- "FBN1"

chr15_genes$GENE <- grepl("MTHFS", chr15_genes$genename)
chr15_genes$GENE_NAME[chr15_genes$GENE == "TRUE"] <- "MTHFS"

chr15_genes$GENE <- grepl("WDR73", chr15_genes$genename)
chr15_genes$GENE_NAME[chr15_genes$GENE == "TRUE"] <- "WDR73"

chr15_genes$GENE <- grepl("ACAN", chr15_genes$genename)
chr15_genes$GENE_NAME[chr15_genes$GENE == "TRUE"] <- "ACAN"

chr15_genes$GENE <- grepl("BLM", chr15_genes$genename)
chr15_genes$GENE_NAME[chr15_genes$GENE == "TRUE"] <- "RECQL3"

chr15_genes$GENE <- grepl("IGF1R", chr15_genes$genename)
chr15_genes$GENE_NAME[chr15_genes$GENE == "TRUE"] <- "IGF1R"






# check each of the genes
unique(chr15_genes$genename[chr15_genes$GENE_NAME == "BUB1B"])
unique(chr15_genes$genename[chr15_genes$GENE_NAME == "KNL1"])
unique(chr15_genes$genename[chr15_genes$GENE_NAME == "UBR1"])
unique(chr15_genes$genename[chr15_genes$GENE_NAME == "TUBGCP4"])
unique(chr15_genes$genename[chr15_genes$GENE_NAME == "CEP152"])
unique(chr15_genes$genename[chr15_genes$GENE_NAME == "FBN1"])
unique(chr15_genes$genename[chr15_genes$GENE_NAME == "MTHFS"])
unique(chr15_genes$genename[chr15_genes$GENE_NAME == "WDR73"])
unique(chr15_genes$genename[chr15_genes$GENE_NAME == "ACAN"])
unique(chr15_genes$genename[chr15_genes$GENE_NAME == "RECQL3"])
unique(chr15_genes$genename[chr15_genes$GENE_NAME == "IGF1R"])


unique(chr15_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr15_variants_genes <- chr15_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr15_variants_genes, "chr15_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr15_variants_genes.csv --path Emily-folder/variants-genes/chr15_variants_genes.csv')









# CHR 16 ------------------------------------------------------------------

# load in the dataframe
chr16_path <- read_csv("path-scored-variants/chr16_path.csv")

# combine chr:pos:ref:alt_alt

chr16_path$variant <- paste(chr16_path$chr,":",chr16_path$pos,":",chr16_path$ref,":",chr16_path$alt,"_",chr16_path$alt, sep = "")

# just need the variant info and gene names
chr16_genes <- chr16_path[c(9,12)]

chr16_genes$GENE_NAME <- chr16_genes$genename



chr16_genes$GENE <- grepl("IGFALS", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "IGFALS"

chr16_genes$GENE <- grepl("NDE1", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "NDE1"

chr16_genes$GENE <- grepl("KIF22", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "KIF22"

chr16_genes$GENE <- grepl("THOC6", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "THOC6"

chr16_genes$GENE <- grepl("CREBBP", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "CREBBP"

chr16_genes$GENE <- grepl("SRCAP", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "SRCAP"

chr16_genes$GENE <- grepl("PAM16", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "PAM16"

chr16_genes$GENE <- grepl("ORC6", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "ORC6"

chr16_genes$GENE <- grepl("RSPRY1", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "RSPRY1"

chr16_genes$GENE <- grepl("KATNB1", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "KATNB1"

chr16_genes$GENE <- grepl("GINS3", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "GINS3"

chr16_genes$GENE <- grepl("CENPT", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "CENPT"

chr16_genes$GENE <- grepl("PRMT7", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "PRMT7"

chr16_genes$GENE <- grepl("VPS4A", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "VPS4A"

chr16_genes$GENE <- grepl("AARS", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "AARS1"

chr16_genes$GENE <- grepl("GINS2", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "GINS2"

chr16_genes$GENE <- grepl("CTU2", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "CTU2"

chr16_genes$GENE <- grepl("CDT1", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "CDT1"

chr16_genes$GENE <- grepl("GALNS", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "GALNS"

chr16_genes$GENE <- grepl("RPL13", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "RPL13"

chr16_genes$GENE <- grepl("CHMP1A", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "CHMP1A"

chr16_genes$GENE <- grepl("FANCA", chr16_genes$genename)
chr16_genes$GENE_NAME[chr16_genes$GENE == "TRUE"] <- "FANCA"







# check each of the genes
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "IGFALS"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "NDE1"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "KIF22"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "THOC6"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "CREBBP"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "SRCAP"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "PAM16"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "ORC6"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "RSPRY1"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "KATNB1"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "GINS3"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "CENPT"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "PRMT7"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "VPS4A"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "AARS1"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "GINS2"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "CTU2"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "CDT1"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "GALNS"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "RPL13"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "CHMP1A"])
unique(chr16_genes$genename[chr16_genes$GENE_NAME == "FANCA"])


unique(chr16_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr16_variants_genes <- chr16_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr16_variants_genes, "chr16_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr16_variants_genes.csv --path Emily-folder/variants-genes/chr16_variants_genes.csv')








# CHR 17 ------------------------------------------------------------------


# load in the dataframe
chr17_path <- read_csv("path-scored-variants/chr17_path.csv")

# combine chr:pos:ref:alt_alt

chr17_path$variant <- paste(chr17_path$chr,":",chr17_path$pos,":",chr17_path$ref,":",chr17_path$alt,"_",chr17_path$alt, sep = "")

# just need the variant info and gene names
chr17_genes <- chr17_path[c(9,12)]

chr17_genes$GENE_NAME <- chr17_genes$genename



chr17_genes$GENE <- grepl("SERPINF1", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "SERPINF1"

chr17_genes$GENE <- grepl("TOP3A", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "TOP3A"

chr17_genes$GENE <- grepl("DPH1", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "DPH1"

chr17_genes$GENE <- grepl("SMARCE1", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "SMARCE1"

chr17_genes$GENE <- grepl("THRA", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "THRA"

chr17_genes$GENE <- grepl("CDC6", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "CDC6"

chr17_genes$GENE <- grepl("FKBP10", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "FKBP10"

chr17_genes$GENE <- grepl("COASY", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "COASY"

chr17_genes$GENE <- grepl("STAT5B", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "STAT5B"

chr17_genes$GENE <- grepl("FZD2", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "FZD2"

chr17_genes$GENE <- grepl("SPOP", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "SPOP"

chr17_genes$GENE <- grepl("COL1A1", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "COL1A1"

chr17_genes$GENE <- grepl("TRIM37", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "TRIM37"

chr17_genes$GENE <- grepl("GH1", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "GH1"

chr17_genes$GENE <- grepl("PRKAR1A", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "PRKAR1A"

chr17_genes$GENE <- grepl("SOX9", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "SOX9"

chr17_genes$GENE <- grepl("EIF5A", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "EIF5A"

chr17_genes$GENE <- grepl("TSEN54", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "TSEN54"

chr17_genes$GENE <- grepl("SLC25A19", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "SLC25A19"

chr17_genes$GENE <- grepl("CANT1", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "CANT1"

chr17_genes$GENE <- grepl("NXN", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "NXN"

chr17_genes$GENE <- grepl("ACTG1", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "ACTG1"

chr17_genes$GENE <- grepl("PYCR1", chr17_genes$genename)
chr17_genes$GENE_NAME[chr17_genes$GENE == "TRUE"] <- "PYCR1"




# check each of the genes
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "SERPINF1"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "TOP3A"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "DPH1"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "SMARCE1"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "THRA"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "CDC6"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "FKBP10"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "COASY"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "STAT5B"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "FZD2"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "SPOP"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "COL1A1"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "TRIM37"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "GH1"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "PRKAR1A"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "SOX9"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "EIF5A"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "TSEN54"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "SLC25A19"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "CANT1"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "NXN"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "ACTG1"])
unique(chr17_genes$genename[chr17_genes$GENE_NAME == "PYCR1"])


unique(chr17_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr17_variants_genes <- chr17_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr17_variants_genes, "chr17_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr17_variants_genes.csv --path Emily-folder/variants-genes/chr17_variants_genes.csv')








# CHR 18 ------------------------------------------------------------------


# load in the dataframe
chr18_path <- read_csv("path-scored-variants/chr18_path.csv")

# combine chr:pos:ref:alt_alt

chr18_path$variant <- paste(chr18_path$chr,":",chr18_path$pos,":",chr18_path$ref,":",chr18_path$alt,"_",chr18_path$alt, sep = "")

# just need the variant info and gene names
chr18_genes <- chr18_path[c(9,12)]

chr18_genes$GENE_NAME <- chr18_genes$genename



chr18_genes$GENE <- grepl("PIEZO2", chr18_genes$genename)
chr18_genes$GENE_NAME[chr18_genes$GENE == "TRUE"] <- "PIEZO2"

chr18_genes$GENE <- grepl("RBBP8", chr18_genes$genename)
chr18_genes$GENE_NAME[chr18_genes$GENE == "TRUE"] <- "RBBP8"

chr18_genes$GENE <- grepl("IER3IP1", chr18_genes$genename)
chr18_genes$GENE_NAME[chr18_genes$GENE == "TRUE"] <- "IER3IP1"

chr18_genes$GENE <- grepl("DYM", chr18_genes$genename)
chr18_genes$GENE_NAME[chr18_genes$GENE == "TRUE"] <- "DYM"

chr18_genes$GENE <- grepl("SMAD4", chr18_genes$genename)
chr18_genes$GENE_NAME[chr18_genes$GENE == "TRUE"] <- "SMAD4"

chr18_genes$GENE <- grepl("RTTN", chr18_genes$genename)
chr18_genes$GENE_NAME[chr18_genes$GENE == "TRUE"] <- "RTTN"






# check each of the genes
unique(chr18_genes$genename[chr18_genes$GENE_NAME == "PIEZO2"])
unique(chr18_genes$genename[chr18_genes$GENE_NAME == "RBBP8"])
unique(chr18_genes$genename[chr18_genes$GENE_NAME == "IER3IP1"])
unique(chr18_genes$genename[chr18_genes$GENE_NAME == "DYM"])
unique(chr18_genes$genename[chr18_genes$GENE_NAME == "SMAD4"])
unique(chr18_genes$genename[chr18_genes$GENE_NAME == "RTTN"])

unique(chr18_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr18_variants_genes <- chr18_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr18_variants_genes, "chr18_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr18_variants_genes.csv --path Emily-folder/variants-genes/chr18_variants_genes.csv')







# CHR 19 ------------------------------------------------------------------

# load in the dataframe
chr19_path <- read_csv("path-scored-variants/chr19_path.csv")

# combine chr:pos:ref:alt_alt

chr19_path$variant <- paste(chr19_path$chr,":",chr19_path$pos,":",chr19_path$ref,":",chr19_path$alt,"_",chr19_path$alt, sep = "")

# just need the variant info and gene names
chr19_genes <- chr19_path[c(9,12)]

chr19_genes$GENE_NAME <- chr19_genes$genename



chr19_genes$GENE <- grepl("GPX4", chr19_genes$genename)
chr19_genes$GENE_NAME[chr19_genes$GENE == "TRUE"] <- "GPX4"

chr19_genes$GENE <- grepl("ACP5", chr19_genes$genename)
chr19_genes$GENE_NAME[chr19_genes$GENE == "TRUE"] <- "ACP5"

chr19_genes$GENE <- grepl("TRMT1", chr19_genes$genename)
chr19_genes$GENE_NAME[chr19_genes$GENE == "TRUE"] <- "TRMT1"

chr19_genes$GENE <- grepl("COMP", chr19_genes$genename)
chr19_genes$GENE_NAME[chr19_genes$GENE == "TRUE"] <- "COMP"

chr19_genes$GENE <- grepl("LMNB2", chr19_genes$genename)
chr19_genes$GENE_NAME[chr19_genes$GENE == "TRUE"] <- "LMNB2"

chr19_genes$GENE <- grepl("WDR62", chr19_genes$genename)
chr19_genes$GENE_NAME[chr19_genes$GENE == "TRUE"] <- "WDR62"

chr19_genes$GENE <- grepl("ERCC2", chr19_genes$genename)
chr19_genes$GENE_NAME[chr19_genes$GENE == "TRUE"] <- "ERCC2"

chr19_genes$GENE <- grepl("CCDC8", chr19_genes$genename)
chr19_genes$GENE_NAME[chr19_genes$GENE == "TRUE"] <- "CCDC8"

chr19_genes$GENE <- grepl("PNKP", chr19_genes$genename)
chr19_genes$GENE_NAME[chr19_genes$GENE == "TRUE"] <- "PNKP"

chr19_genes$GENE <- grepl("INSR", chr19_genes$genename)
chr19_genes$GENE_NAME[chr19_genes$GENE == "TRUE"] <- "INSR"

chr19_genes$GENE <- grepl("PNPLA6", chr19_genes$genename)
chr19_genes$GENE_NAME[chr19_genes$GENE == "TRUE"] <- "PNPLA6"




# check each of the genes
unique(chr19_genes$genename[chr19_genes$GENE_NAME == "GPX4"])
unique(chr19_genes$genename[chr19_genes$GENE_NAME == "ACP5"])
unique(chr19_genes$genename[chr19_genes$GENE_NAME == "TRMT1"])
unique(chr19_genes$genename[chr19_genes$GENE_NAME == "COMP"])
unique(chr19_genes$genename[chr19_genes$GENE_NAME == "LMNB2"])
unique(chr19_genes$genename[chr19_genes$GENE_NAME == "WDR62"])
unique(chr19_genes$genename[chr19_genes$GENE_NAME == "ERCC2"])
unique(chr19_genes$genename[chr19_genes$GENE_NAME == "CCDC8"])
unique(chr19_genes$genename[chr19_genes$GENE_NAME == "PNKP"])
unique(chr19_genes$genename[chr19_genes$GENE_NAME == "INSR"])
unique(chr19_genes$genename[chr19_genes$GENE_NAME == "PNPLA6"])


unique(chr19_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr19_variants_genes <- chr19_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr19_variants_genes, "chr19_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr19_variants_genes.csv --path Emily-folder/variants-genes/chr19_variants_genes.csv')








# CHR 20 ------------------------------------------------------------------

# load in the dataframe
chr20_path <- read_csv("path-scored-variants/chr20_path.csv")

# combine chr:pos:ref:alt_alt

chr20_path$variant <- paste(chr20_path$chr,":",chr20_path$pos,":",chr20_path$ref,":",chr20_path$alt,"_",chr20_path$alt, sep = "")

# just need the variant info and gene names
chr20_genes <- chr20_path[c(9,12)]

chr20_genes$GENE_NAME <- chr20_genes$genename



chr20_genes$GENE <- grepl("DDRGK1", chr20_genes$genename)
chr20_genes$GENE_NAME[chr20_genes$GENE == "TRUE"] <- "DDRGK1"

chr20_genes$GENE <- grepl("GDF5", chr20_genes$genename)
chr20_genes$GENE_NAME[chr20_genes$GENE == "TRUE"] <- "GDF5"

chr20_genes$GENE <- grepl("GHRH", chr20_genes$genename)
chr20_genes$GENE_NAME[chr20_genes$GENE == "TRUE"] <- "GHRH"

chr20_genes$GENE <- grepl("ZNF335", chr20_genes$genename)
chr20_genes$GENE_NAME[chr20_genes$GENE == "TRUE"] <- "ZNF335"

chr20_genes$GENE <- grepl("CTSA", chr20_genes$genename)
chr20_genes$GENE_NAME[chr20_genes$GENE == "TRUE"] <- "CTSA"

chr20_genes$GENE <- grepl("TP53RK", chr20_genes$genename)
chr20_genes$GENE_NAME[chr20_genes$GENE == "TRUE"] <- "TP53RK"

chr20_genes$GENE <- grepl("GNAS", chr20_genes$genename)
chr20_genes$GENE_NAME[chr20_genes$GENE == "TRUE"] <- "GNAS"

chr20_genes$GENE <- grepl("BMP2", chr20_genes$genename)
chr20_genes$GENE_NAME[chr20_genes$GENE == "TRUE"] <- "BMP2"

chr20_genes$GENE <- grepl("COL9A3", chr20_genes$genename)
chr20_genes$GENE_NAME[chr20_genes$GENE == "TRUE"] <- "COL9A3"



# check each of the genes
unique(chr20_genes$genename[chr20_genes$GENE_NAME == "DDRGK1"])
unique(chr20_genes$genename[chr20_genes$GENE_NAME == "GDF5"])
unique(chr20_genes$genename[chr20_genes$GENE_NAME == "GHRH"])
unique(chr20_genes$genename[chr20_genes$GENE_NAME == "ZNF335"])
unique(chr20_genes$genename[chr20_genes$GENE_NAME == "CTSA"])
unique(chr20_genes$genename[chr20_genes$GENE_NAME == "TP53RK"])
unique(chr20_genes$genename[chr20_genes$GENE_NAME == "GNAS"])
unique(chr20_genes$genename[chr20_genes$GENE_NAME == "BMP2"])
unique(chr20_genes$genename[chr20_genes$GENE_NAME == "COL9A3"])


unique(chr20_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr20_variants_genes <- chr20_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr20_variants_genes, "chr20_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr20_variants_genes.csv --path Emily-folder/variants-genes/chr20_variants_genes.csv')









# CHR 21 ------------------------------------------------------------------

# load in the dataframe
chr21_path <- read_csv("path-scored-variants/chr21_path.csv")

# combine chr:pos:ref:alt_alt

chr21_path$variant <- paste(chr21_path$chr,":",chr21_path$pos,":",chr21_path$ref,":",chr21_path$alt,"_",chr21_path$alt, sep = "")

# just need the variant info and gene names
chr21_genes <- chr21_path[c(9,12)]

chr21_genes$GENE_NAME <- chr21_genes$genename



chr21_genes$GENE <- grepl("DONSON", chr21_genes$genename)
chr21_genes$GENE_NAME[chr21_genes$GENE == "TRUE"] <- "DONSON"

chr21_genes$GENE <- grepl("DYRK1A", chr21_genes$genename)
chr21_genes$GENE_NAME[chr21_genes$GENE == "TRUE"] <- "DYRK1A"

chr21_genes$GENE <- grepl("WDR4", chr21_genes$genename)
chr21_genes$GENE_NAME[chr21_genes$GENE == "TRUE"] <- "WDR4"

chr21_genes$GENE <- grepl("CFAP410", chr21_genes$genename)
chr21_genes$GENE_NAME[chr21_genes$GENE == "TRUE"] <- "CFAP410"

chr21_genes$GENE <- grepl("PCNT", chr21_genes$genename)
chr21_genes$GENE_NAME[chr21_genes$GENE == "TRUE"] <- "PCNT"


# check each of the genes
unique(chr21_genes$genename[chr21_genes$GENE_NAME == "DONSON"])
unique(chr21_genes$genename[chr21_genes$GENE_NAME == "DYRK1A"])
unique(chr21_genes$genename[chr21_genes$GENE_NAME == "WDR4"])
unique(chr21_genes$genename[chr21_genes$GENE_NAME == "CFAP410"])
unique(chr21_genes$genename[chr21_genes$GENE_NAME == "PCNT"])



unique(chr21_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr21_variants_genes <- chr21_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr21_variants_genes, "chr21_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr21_variants_genes.csv --path Emily-folder/variants-genes/chr21_variants_genes.csv')







# CHR 22 ------------------------------------------------------------------


# load in the dataframe
chr22_path <- read_csv("path-scored-variants/chr22_path.csv")

# combine chr:pos:ref:alt_alt

chr22_path$variant <- paste(chr22_path$chr,":",chr22_path$pos,":",chr22_path$ref,":",chr22_path$alt,"_",chr22_path$alt, sep = "")

# just need the variant info and gene names
chr22_genes <- chr22_path[c(9,12)]

chr22_genes$GENE_NAME <- chr22_genes$genename



chr22_genes$GENE <- grepl("CDC45", chr22_genes$genename)
chr22_genes$GENE_NAME[chr22_genes$GENE == "TRUE"] <- "CDC45"

chr22_genes$GENE <- grepl("LZTR1", chr22_genes$genename)
chr22_genes$GENE_NAME[chr22_genes$GENE == "TRUE"] <- "LZTR1"

chr22_genes$GENE <- grepl("PISD", chr22_genes$genename)
chr22_genes$GENE_NAME[chr22_genes$GENE == "TRUE"] <- "PISD"

chr22_genes$GENE <- grepl("MCM5", chr22_genes$genename)
chr22_genes$GENE_NAME[chr22_genes$GENE == "TRUE"] <- "MCM5"

chr22_genes$GENE <- grepl("EP300", chr22_genes$genename)
chr22_genes$GENE_NAME[chr22_genes$GENE == "TRUE"] <- "EP300"

chr22_genes$GENE <- grepl("RRP7A", chr22_genes$genename)
chr22_genes$GENE_NAME[chr22_genes$GENE == "TRUE"] <- "RRP7A"

chr22_genes$GENE <- grepl("SBF1", chr22_genes$genename)
chr22_genes$GENE_NAME[chr22_genes$GENE == "TRUE"] <- "SBF1"

chr22_genes$GENE <- grepl("TUBGCP6", chr22_genes$genename)
chr22_genes$GENE_NAME[chr22_genes$GENE == "TRUE"] <- "TUBGCP6"




# check each of the genes
unique(chr22_genes$genename[chr22_genes$GENE_NAME == "CDC45"])
unique(chr22_genes$genename[chr22_genes$GENE_NAME == "LZTR1"])
unique(chr22_genes$genename[chr22_genes$GENE_NAME == "PISD"])
unique(chr22_genes$genename[chr22_genes$GENE_NAME == "MCM5"])
unique(chr22_genes$genename[chr22_genes$GENE_NAME == "EP300"])
unique(chr22_genes$genename[chr22_genes$GENE_NAME == "RRP7A"])
unique(chr22_genes$genename[chr22_genes$GENE_NAME == "SBF1"])
unique(chr22_genes$genename[chr22_genes$GENE_NAME == "TUBGCP6"])


unique(chr22_genes$GENE_NAME)
# all looks good

# just save the variant and GENE_NAME columns
chr22_variants_genes <- chr22_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chr22_variants_genes, "chr22_variants_genes.csv")
# now to save this to UKB project
system('dx upload chr22_variants_genes.csv --path Emily-folder/variants-genes/chr22_variants_genes.csv')






# CHR X -------------------------------------------------------------------


# load in the dataframe
chrX_path <- read_csv("path-scored-variants/chrX_path.csv")

# combine chr:pos:ref:alt_alt

chrX_path$variant <- paste(chrX_path$chr,":",chrX_path$pos,":",chrX_path$ref,":",chrX_path$alt,"_",chrX_path$alt, sep = "")

# just need the variant info and gene names
chrX_genes <- chrX_path[c(9,12)]

chrX_genes$GENE_NAME <- chrX_genes$genename



chrX_genes$GENE <- grepl("BTK", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "BTK"

chrX_genes$GENE <- grepl("AMMECR1", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "AMMECR1"

chrX_genes$GENE <- grepl("AIFM1", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "AIFM1"

chrX_genes$GENE <- grepl("HMGB3", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "HMGB3"

chrX_genes$GENE <- grepl("BGN", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "BGN"

chrX_genes$GENE <- grepl("FLNA", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "FLNA"

chrX_genes$GENE <- grepl("LAGE3", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "LAGE3"

chrX_genes$GENE <- grepl("PDHA1", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "PDHA1"

chrX_genes$GENE <- grepl("MBTPS2", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "MBTPS2"

chrX_genes$GENE <- grepl("PHEX", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "PHEX"

chrX_genes$GENE <- grepl("EIF2S3", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "EIF2S3"

chrX_genes$GENE <- grepl("ARX", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "ARX"

chrX_genes$GENE <- grepl("ARSE", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "ARSL"

chrX_genes$GENE <- grepl("CASK", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "CASK"

chrX_genes$GENE <- grepl("KDM6A", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "KDM6A"

chrX_genes$GENE <- grepl("EBP", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "EBP"

chrX_genes$GENE <- grepl("PQBP1", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "PQBP1"

chrX_genes$GENE <- grepl("SMC1A", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "SMC1A"

chrX_genes$GENE <- grepl("FGD1", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "FGD1"

chrX_genes$GENE <- grepl("SHOX", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "SHOX"

chrX_genes$GENE <- grepl("HDAC8", chrX_genes$genename)
chrX_genes$GENE_NAME[chrX_genes$GENE == "TRUE"] <- "HDAC8"






# check each of the genes
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "BTK"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "AMMECR1"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "AIFM1"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "HMGB3"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "BGN"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "FLNA"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "LAGE3"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "PDHA1"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "MBTPS2"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "PHEX"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "EIF2S3"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "ARX"])
#none
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "ARSL"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "CASK"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "KDM6A"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "EBP"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "PQBP1"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "SMC1A"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "FGD1"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "SHOX"])
unique(chrX_genes$genename[chrX_genes$GENE_NAME == "HDAC8"])


unique(chrX_genes$GENE_NAME)
# all looks good - still missing ARX

# just save the variant and GENE_NAME columns
chrX_variants_genes <- chrX_genes[c(2,3)]

# save this for later
#save it to R
write.csv(chrX_variants_genes, "chrX_variants_genes.csv")
# now to save this to UKB project
system('dx upload chrX_variants_genes.csv --path Emily-folder/variants-genes/chrX_variants_genes.csv')









system('dx upload variants_genes_script.R --path Emily-folder/variants-genes/variants_genes_script.R')