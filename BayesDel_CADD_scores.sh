#execute this script when within the dbNSFP wd to obtain BayesDel addAF and CADD phred scores for a folder of VCF files
for f in /mnt/hcs/dsm/pathology/bicknel/Emily/SNPs/*; do java search_dbNSFP44a -v hg38 -g -w 1-7,12-13,102,130 -i $f -o ${f%.vcf}_output.vcf; done
