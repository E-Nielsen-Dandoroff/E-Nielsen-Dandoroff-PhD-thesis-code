#!/bin/bash

## execute this first part in the terminal before running the script:
# # Download the SNP information from the project folder
# dx download -r Emily-folder/path-scored-variants
# 
# # Download and extract plink software
## wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20221210.zip
## this script doesn't work in R terminal
## instead run this in console:
# download.file("https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20221210.zip", destfile = "/home/rstudio-server/plink_linux_x86_64_20221210.zip")
## the rest can be run in the terminal:
# mkdir plink_linux_x86_64_20221210
# mv plink_linux_x86_64_20221210.zip plink_linux_x86_64_20221210
# cd plink_linux_x86_64_20221210
# unzip plink_linux_x86_64_20221210.zip
# cd ..
# 
# # Set plink up so that it runs using $PLINK variable
# PLINK=$PWD/plink_linux_x86_64_20221210/plink
# 
# # Check version (just to make sure it is working)
# $PLINK --version

##########
# chr1
##########

# Create a bed file for extracting chr1 SNPs 
tail -n +2 path-scored-variants/chr1_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr1_path_distinct.bed
head chr1_path_distinct.bed

# Make a directory for all of the chr1 SNP data
mkdir chr1-snp-data

# Use plink to extract the chr1 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c1_b0_v1 --extract range chr1_path_distinct.bed --recode A --out chr1-snp-data/chr1-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr1-snp-data/chr1-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c1_b0_v1 --extract range chr1_path_distinct.bed --recode A --out chr1-snp-data/chr1-snps --freqx
# Move the allele frequency data back to the chr1 output directory
mv chr1-snps.frq chr1-snp-data/
  
  # have a look at the chr1 frequency data
  head chr1-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr1_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr1-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr1-snp-ids.txt chr1-snp-data/chr1-snps.frq | head
grep -f chr1-snp-ids.txt chr1-snp-data/chr1-snps.frqx | head





##########
# chr2
##########

# Create a bed file for extracting chr2 SNPs 
tail -n +2 path-scored-variants/chr2_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr2_path_distinct.bed
head chr2_path_distinct.bed

# Make a directory for all of the chr2 SNP data
mkdir chr2-snp-data

# Use plink to extract the chr2 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c2_b0_v1 --extract range chr2_path_distinct.bed --recode A --out chr2-snp-data/chr2-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr2-snp-data/chr2-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c2_b0_v1 --extract range chr2_path_distinct.bed --recode A --out chr2-snp-data/chr2-snps --freqx
# Move the allele frequency data back to the chr2 output directory
mv chr2-snps.frq chr2-snp-data/
  
  # have a look at the chr2 frequency data
  head chr2-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr2_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr2-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr2-snp-ids.txt chr2-snp-data/chr2-snps.frq | head
grep -f chr2-snp-ids.txt chr2-snp-data/chr2-snps.frqx | head





##########
# chr3
##########

# Create a bed file for extracting chr3 SNPs 
tail -n +2 path-scored-variants/chr3_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr3_path_distinct.bed
head chr3_path_distinct.bed

# Make a directory for all of the chr3 SNP data
mkdir chr3-snp-data

# Use plink to extract the chr3 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c3_b0_v1 --extract range chr3_path_distinct.bed --recode A --out chr3-snp-data/chr3-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr3-snp-data/chr3-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c3_b0_v1 --extract range chr3_path_distinct.bed --recode A --out chr3-snp-data/chr3-snps --freqx
# Move the allele frequency data back to the chr3 output directory
mv chr3-snps.frq chr3-snp-data/
  
  # have a look at the chr3 frequency data
  head chr3-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr3_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr3-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr3-snp-ids.txt chr3-snp-data/chr3-snps.frq | head
grep -f chr3-snp-ids.txt chr3-snp-data/chr3-snps.frqx | head







##########
# chr4
##########

# Create a bed file for extracting chr4 SNPs 
tail -n +2 path-scored-variants/chr4_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr4_path_distinct.bed
head chr4_path_distinct.bed

# Make a directory for all of the chr4 SNP data
mkdir chr4-snp-data

# Use plink to extract the chr4 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c4_b0_v1 --extract range chr4_path_distinct.bed --recode A --out chr4-snp-data/chr4-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr4-snp-data/chr4-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c4_b0_v1 --extract range chr4_path_distinct.bed --recode A --out chr4-snp-data/chr4-snps --freqx
# Move the allele frequency data back to the chr4 output directory
mv chr4-snps.frq chr4-snp-data/
  
  # have a look at the chr4 frequency data
  head chr4-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr4_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr4-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr4-snp-ids.txt chr4-snp-data/chr4-snps.frq | head
grep -f chr4-snp-ids.txt chr4-snp-data/chr4-snps.frqx | head









##########
# chr5
##########

# Create a bed file for extracting chr5 SNPs 
tail -n +2 path-scored-variants/chr5_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr5_path_distinct.bed
head chr5_path_distinct.bed

# Make a directory for all of the chr5 SNP data
mkdir chr5-snp-data

# Use plink to extract the chr5 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c5_b0_v1 --extract range chr5_path_distinct.bed --recode A --out chr5-snp-data/chr5-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr5-snp-data/chr5-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c5_b0_v1 --extract range chr5_path_distinct.bed --recode A --out chr5-snp-data/chr5-snps --freqx
# Move the allele frequency data back to the chr5 output directory
mv chr5-snps.frq chr5-snp-data/
  
  # have a look at the chr5 frequency data
  head chr5-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr5_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr5-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr5-snp-ids.txt chr5-snp-data/chr5-snps.frq | head
grep -f chr5-snp-ids.txt chr5-snp-data/chr5-snps.frqx | head







##########
# chr6
##########

# Create a bed file for extracting chr6 SNPs 
tail -n +2 path-scored-variants/chr6_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr6_path_distinct.bed
head chr6_path_distinct.bed

# Make a directory for all of the chr6 SNP data
mkdir chr6-snp-data

# Use plink to extract the chr6 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c6_b0_v1 --extract range chr6_path_distinct.bed --recode A --out chr6-snp-data/chr6-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr6-snp-data/chr6-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c6_b0_v1 --extract range chr6_path_distinct.bed --recode A --out chr6-snp-data/chr6-snps --freqx
# Move the allele frequency data back to the chr6 output directory
mv chr6-snps.frq chr6-snp-data/
  
  # have a look at the chr6 frequency data
  head chr6-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr6_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr6-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr6-snp-ids.txt chr6-snp-data/chr6-snps.frq | head
grep -f chr6-snp-ids.txt chr6-snp-data/chr6-snps.frqx | head








##########
# chr7
##########

# Create a bed file for extracting chr7 SNPs 
tail -n +2 path-scored-variants/chr7_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr7_path_distinct.bed
head chr7_path_distinct.bed

# Make a directory for all of the chr7 SNP data
mkdir chr7-snp-data

# Use plink to extract the chr7 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c7_b0_v1 --extract range chr7_path_distinct.bed --recode A --out chr7-snp-data/chr7-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr7-snp-data/chr7-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c7_b0_v1 --extract range chr7_path_distinct.bed --recode A --out chr7-snp-data/chr7-snps --freqx
# Move the allele frequency data back to the chr7 output directory
mv chr7-snps.frq chr7-snp-data/
  
  # have a look at the chr7 frequency data
  head chr7-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr7_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr7-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr7-snp-ids.txt chr7-snp-data/chr7-snps.frq | head
grep -f chr7-snp-ids.txt chr7-snp-data/chr7-snps.frqx | head








##########
# chr8
##########

# Create a bed file for extracting chr8 SNPs 
tail -n +2 path-scored-variants/chr8_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr8_path_distinct.bed
head chr8_path_distinct.bed

# Make a directory for all of the chr8 SNP data
mkdir chr8-snp-data

# Use plink to extract the chr8 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c8_b0_v1 --extract range chr8_path_distinct.bed --recode A --out chr8-snp-data/chr8-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr8-snp-data/chr8-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c8_b0_v1 --extract range chr8_path_distinct.bed --recode A --out chr8-snp-data/chr8-snps --freqx
# Move the allele frequency data back to the chr8 output directory
mv chr8-snps.frq chr8-snp-data/
  
  # have a look at the chr8 frequency data
  head chr8-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr8_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr8-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr8-snp-ids.txt chr8-snp-data/chr8-snps.frq | head
grep -f chr8-snp-ids.txt chr8-snp-data/chr8-snps.frqx | head






##########
# chr9
##########

# Create a bed file for extracting chr9 SNPs 
tail -n +2 path-scored-variants/chr9_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr9_path_distinct.bed
head chr9_path_distinct.bed

# Make a directory for all of the chr9 SNP data
mkdir chr9-snp-data

# Use plink to extract the chr9 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c9_b0_v1 --extract range chr9_path_distinct.bed --recode A --out chr9-snp-data/chr9-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr9-snp-data/chr9-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c9_b0_v1 --extract range chr9_path_distinct.bed --recode A --out chr9-snp-data/chr9-snps --freqx
# Move the allele frequency data back to the chr9 output directory
mv chr9-snps.frq chr9-snp-data/
  
  # have a look at the chr9 frequency data
  head chr9-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr9_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr9-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr9-snp-ids.txt chr9-snp-data/chr9-snps.frq | head
grep -f chr9-snp-ids.txt chr9-snp-data/chr9-snps.frqx | head







##########
# chr10
##########

# Create a bed file for extracting chr10 SNPs 
tail -n +2 path-scored-variants/chr10_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr10_path_distinct.bed
head chr10_path_distinct.bed

# Make a directory for all of the chr10 SNP data
mkdir chr10-snp-data

# Use plink to extract the chr10 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c10_b0_v1 --extract range chr10_path_distinct.bed --recode A --out chr10-snp-data/chr10-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr10-snp-data/chr10-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c10_b0_v1 --extract range chr10_path_distinct.bed --recode A --out chr10-snp-data/chr10-snps --freqx
# Move the allele frequency data back to the chr10 output directory
mv chr10-snps.frq chr10-snp-data/
  
  # have a look at the chr10 frequency data
  head chr10-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr10_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr10-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr10-snp-ids.txt chr10-snp-data/chr10-snps.frq | head
grep -f chr10-snp-ids.txt chr10-snp-data/chr10-snps.frqx | head







##########
# chr11
##########

# Create a bed file for extracting chr11 SNPs 
tail -n +2 path-scored-variants/chr11_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr11_path_distinct.bed
head chr11_path_distinct.bed

# Make a directory for all of the chr11 SNP data
mkdir chr11-snp-data

# Use plink to extract the chr11 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c11_b0_v1 --extract range chr11_path_distinct.bed --recode A --out chr11-snp-data/chr11-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr11-snp-data/chr11-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c11_b0_v1 --extract range chr11_path_distinct.bed --recode A --out chr11-snp-data/chr11-snps --freqx
# Move the allele frequency data back to the chr11 output directory
mv chr11-snps.frq chr11-snp-data/
  
  # have a look at the chr11 frequency data
  head chr11-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr11_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr11-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr11-snp-ids.txt chr11-snp-data/chr11-snps.frq | head
grep -f chr11-snp-ids.txt chr11-snp-data/chr11-snps.frqx | head






##########
# chr12
##########

# Create a bed file for extracting chr12 SNPs 
tail -n +2 path-scored-variants/chr12_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr12_path_distinct.bed
head chr12_path_distinct.bed

# Make a directory for all of the chr12 SNP data
mkdir chr12-snp-data

# Use plink to extract the chr12 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c12_b0_v1 --extract range chr12_path_distinct.bed --recode A --out chr12-snp-data/chr12-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr12-snp-data/chr12-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c12_b0_v1 --extract range chr12_path_distinct.bed --recode A --out chr12-snp-data/chr12-snps --freqx
# Move the allele frequency data back to the chr12 output directory
mv chr12-snps.frq chr12-snp-data/
  
  # have a look at the chr12 frequency data
  head chr12-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr12_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr12-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr12-snp-ids.txt chr12-snp-data/chr12-snps.frq | head
grep -f chr12-snp-ids.txt chr12-snp-data/chr12-snps.frqx | head









##########
# chr13
##########

# Create a bed file for extracting chr13 SNPs 
tail -n +2 path-scored-variants/chr13_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr13_path_distinct.bed
head chr13_path_distinct.bed

# Make a directory for all of the chr13 SNP data
mkdir chr13-snp-data

# Use plink to extract the chr13 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c13_b0_v1 --extract range chr13_path_distinct.bed --recode A --out chr13-snp-data/chr13-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr13-snp-data/chr13-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c13_b0_v1 --extract range chr13_path_distinct.bed --recode A --out chr13-snp-data/chr13-snps --freqx
# Move the allele frequency data back to the chr13 output directory
mv chr13-snps.frq chr13-snp-data/
  
  # have a look at the chr13 frequency data
  head chr13-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr13_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr13-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr13-snp-ids.txt chr13-snp-data/chr13-snps.frq | head
grep -f chr13-snp-ids.txt chr13-snp-data/chr13-snps.frqx | head







##########
# chr14
##########

# Create a bed file for extracting chr14 SNPs 
tail -n +2 path-scored-variants/chr14_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr14_path_distinct.bed
head chr14_path_distinct.bed

# Make a directory for all of the chr14 SNP data
mkdir chr14-snp-data

# Use plink to extract the chr14 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c14_b0_v1 --extract range chr14_path_distinct.bed --recode A --out chr14-snp-data/chr14-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr14-snp-data/chr14-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c14_b0_v1 --extract range chr14_path_distinct.bed --recode A --out chr14-snp-data/chr14-snps --freqx
# Move the allele frequency data back to the chr14 output directory
mv chr14-snps.frq chr14-snp-data/
  
  # have a look at the chr14 frequency data
  head chr14-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr14_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr14-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr14-snp-ids.txt chr14-snp-data/chr14-snps.frq | head
grep -f chr14-snp-ids.txt chr14-snp-data/chr14-snps.frqx | head








##########
# chr15
##########

# Create a bed file for extracting chr15 SNPs 
tail -n +2 path-scored-variants/chr15_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr15_path_distinct.bed
head chr15_path_distinct.bed

# Make a directory for all of the chr15 SNP data
mkdir chr15-snp-data

# Use plink to extract the chr15 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c15_b0_v1 --extract range chr15_path_distinct.bed --recode A --out chr15-snp-data/chr15-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr15-snp-data/chr15-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c15_b0_v1 --extract range chr15_path_distinct.bed --recode A --out chr15-snp-data/chr15-snps --freqx
# Move the allele frequency data back to the chr15 output directory
mv chr15-snps.frq chr15-snp-data/
  
  # have a look at the chr15 frequency data
  head chr15-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr15_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr15-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr15-snp-ids.txt chr15-snp-data/chr15-snps.frq | head
grep -f chr15-snp-ids.txt chr15-snp-data/chr15-snps.frqx | head








##########
# chr16
##########

# Create a bed file for extracting chr16 SNPs 
tail -n +2 path-scored-variants/chr16_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr16_path_distinct.bed
head chr16_path_distinct.bed

# Make a directory for all of the chr16 SNP data
mkdir chr16-snp-data

# Use plink to extract the chr16 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c16_b0_v1 --extract range chr16_path_distinct.bed --recode A --out chr16-snp-data/chr16-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr16-snp-data/chr16-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c16_b0_v1 --extract range chr16_path_distinct.bed --recode A --out chr16-snp-data/chr16-snps --freqx
# Move the allele frequency data back to the chr16 output directory
mv chr16-snps.frq chr16-snp-data/
  
  # have a look at the chr16 frequency data
  head chr16-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr16_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr16-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr16-snp-ids.txt chr16-snp-data/chr16-snps.frq | head
grep -f chr16-snp-ids.txt chr16-snp-data/chr16-snps.frqx | head







##########
# chr17
##########

# Create a bed file for extracting chr17 SNPs 
tail -n +2 path-scored-variants/chr17_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr17_path_distinct.bed
head chr17_path_distinct.bed

# Make a directory for all of the chr17 SNP data
mkdir chr17-snp-data

# Use plink to extract the chr17 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c17_b0_v1 --extract range chr17_path_distinct.bed --recode A --out chr17-snp-data/chr17-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr17-snp-data/chr17-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c17_b0_v1 --extract range chr17_path_distinct.bed --recode A --out chr17-snp-data/chr17-snps --freqx
# Move the allele frequency data back to the chr17 output directory
mv chr17-snps.frq chr17-snp-data/
  
  # have a look at the chr17 frequency data
  head chr17-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr17_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr17-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr17-snp-ids.txt chr17-snp-data/chr17-snps.frq | head
grep -f chr17-snp-ids.txt chr17-snp-data/chr17-snps.frqx | head






##########
# chr18
##########

# Create a bed file for extracting chr18 SNPs 
tail -n +2 path-scored-variants/chr18_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr18_path_distinct.bed
head chr18_path_distinct.bed

# Make a directory for all of the chr18 SNP data
mkdir chr18-snp-data

# Use plink to extract the chr18 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c18_b0_v1 --extract range chr18_path_distinct.bed --recode A --out chr18-snp-data/chr18-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr18-snp-data/chr18-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c18_b0_v1 --extract range chr18_path_distinct.bed --recode A --out chr18-snp-data/chr18-snps --freqx
# Move the allele frequency data back to the chr18 output directory
mv chr18-snps.frq chr18-snp-data/
  
  # have a look at the chr18 frequency data
  head chr18-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr18_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr18-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr18-snp-ids.txt chr18-snp-data/chr18-snps.frq | head
grep -f chr18-snp-ids.txt chr18-snp-data/chr18-snps.frqx | head





##########
# chr19
##########

# Create a bed file for extracting chr19 SNPs 
tail -n +2 path-scored-variants/chr19_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr19_path_distinct.bed
head chr19_path_distinct.bed

# Make a directory for all of the chr19 SNP data
mkdir chr19-snp-data

# Use plink to extract the chr19 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c19_b0_v1 --extract range chr19_path_distinct.bed --recode A --out chr19-snp-data/chr19-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr19-snp-data/chr19-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c19_b0_v1 --extract range chr19_path_distinct.bed --recode A --out chr19-snp-data/chr19-snps --freqx
# Move the allele frequency data back to the chr19 output directory
mv chr19-snps.frq chr19-snp-data/
  
  # have a look at the chr19 frequency data
  head chr19-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr19_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr19-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr19-snp-ids.txt chr19-snp-data/chr19-snps.frq | head
grep -f chr19-snp-ids.txt chr19-snp-data/chr19-snps.frqx | head






##########
# chr20
##########

# Create a bed file for extracting chr20 SNPs 
tail -n +2 path-scored-variants/chr20_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr20_path_distinct.bed
head chr20_path_distinct.bed

# Make a directory for all of the chr20 SNP data
mkdir chr20-snp-data

# Use plink to extract the chr20 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c20_b0_v1 --extract range chr20_path_distinct.bed --recode A --out chr20-snp-data/chr20-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr20-snp-data/chr20-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c20_b0_v1 --extract range chr20_path_distinct.bed --recode A --out chr20-snp-data/chr20-snps --freqx
# Move the allele frequency data back to the chr20 output directory
mv chr20-snps.frq chr20-snp-data/
  
  # have a look at the chr20 frequency data
  head chr20-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr20_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr20-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr20-snp-ids.txt chr20-snp-data/chr20-snps.frq | head
grep -f chr20-snp-ids.txt chr20-snp-data/chr20-snps.frqx | head








##########
# chr21
##########

# Create a bed file for extracting chr21 SNPs 
tail -n +2 path-scored-variants/chr21_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr21_path_distinct.bed
head chr21_path_distinct.bed

# Make a directory for all of the chr21 SNP data
mkdir chr21-snp-data

# Use plink to extract the chr21 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c21_b0_v1 --extract range chr21_path_distinct.bed --recode A --out chr21-snp-data/chr21-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr21-snp-data/chr21-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c21_b0_v1 --extract range chr21_path_distinct.bed --recode A --out chr21-snp-data/chr21-snps --freqx
# Move the allele frequency data back to the chr21 output directory
mv chr21-snps.frq chr21-snp-data/
  
  # have a look at the chr21 frequency data
  head chr21-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr21_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr21-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr21-snp-ids.txt chr21-snp-data/chr21-snps.frq | head
grep -f chr21-snp-ids.txt chr21-snp-data/chr21-snps.frqx | head







##########
# chr22
##########

# Create a bed file for extracting chr22 SNPs 
tail -n +2 path-scored-variants/chr22_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chr22_path_distinct.bed
head chr22_path_distinct.bed

# Make a directory for all of the chr22 SNP data
mkdir chr22-snp-data

# Use plink to extract the chr22 SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c22_b0_v1 --extract range chr22_path_distinct.bed --recode A --out chr22-snp-data/chr22-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chr22-snp-data/chr22-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_c22_b0_v1 --extract range chr22_path_distinct.bed --recode A --out chr22-snp-data/chr22-snps --freqx
# Move the allele frequency data back to the chr22 output directory
mv chr22-snps.frq chr22-snp-data/
  
  # have a look at the chr22 frequency data
  head chr22-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chr22_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chr22-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chr22-snp-ids.txt chr22-snp-data/chr22-snps.frq | head
grep -f chr22-snp-ids.txt chr22-snp-data/chr22-snps.frqx | head




##########
# chrX
##########

# Create a bed file for extracting chrX SNPs 
tail -n +2 path-scored-variants/chrX_path.csv | awk -F ',' '{print $1, $2, $2, $4}' > chrX_path_distinct.bed
head chrX_path_distinct.bed

# Make a directory for all of the chrX SNP data
mkdir chrX-snp-data

# Use plink to extract the chrX SNP data and calculate allele frequencies (--freq)
# Note that the exome data is already mounted in the Jupyter environment, so no need to download it
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_cX_b0_v1 --extract range chrX_path_distinct.bed --recode A --out chrX-snp-data/chrX-snps --freq
# Move the allele frequency data so that it doesn't get overwritten by the next step
mv chrX-snp-data/chrX-snps.frq .

# Re-run the above, and get count-based allele information per SNP (--freqx)
# NB: can't run --freq and --freqx at the same time, so need to run plink twice
$PLINK --bfile /mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ PLINK\ format\ -\ final\ release/ukb23158_cX_b0_v1 --extract range chrX_path_distinct.bed --recode A --out chrX-snp-data/chrX-snps --freqx
# Move the allele frequency data back to the chrX output directory
mv chrX-snps.frq chrX-snp-data/
  
  # have a look at the chrX frequency data
  head chrX-snp-data/*.f*
  
  # create plink-style SNP ids so that only the SNPs of interest can be extracted
  cat path-scored-variants/chrX_path.csv | awk -F ',' '{print $1 ":" $2 ":" $3 ":" $4}' > chrX-snp-ids.txt

# look at allele frequency info for SNPs of interest
grep -f chrX-snp-ids.txt chrX-snp-data/chrX-snps.frq | head
grep -f chrX-snp-ids.txt chrX-snp-data/chrX-snps.frqx | head





##################
# Upload data
##################

# Create date-based folder and move results there
DATE=20231118 
mkdir plink-snp-data-$DATE
mv chr* plink-snp-data-$DATE

# Upload back to the project
dx upload -r plink-snp-data-$DATE

# Put things back how they were
mv plink-snp-data-$DATE/* .
rmdir plink-snp-data-$DATE