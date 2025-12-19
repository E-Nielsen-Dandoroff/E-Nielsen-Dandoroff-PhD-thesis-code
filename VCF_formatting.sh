## VCF format
# this script will create an output text file containing VCF formated SNP data for each text file in the wd
# make sure the wd only contains raw/unprocessed SNP files
for variable in *.txt; do sed 's/_/\t/g' $variable | awk '{print $1,$2,".",$3,$4,".",".","."}' | sed 's/ /\t/g' >d_${variable}; done

## add head
# this script will add the VCF format header to any text document in the wd with the name output*.txt
for variable in d_*.txt; do echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' | cat - $variable > n${variable}; done

## add file format
# this script will add the ## header for VCF format to all files named nd_*.txt
for variable in nd_*.txt; do echo -e '##fileformat=VCFv4.2' | cat - $variable > e${variable}; done

## txt to vcf
# use this script last to change extension from .txt to .vcf
for file in end_*.txt; do mv -- "$file" "${file%.txt}.vcf"; done