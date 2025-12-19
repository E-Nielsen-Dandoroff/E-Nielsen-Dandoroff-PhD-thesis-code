## use this script to convert all end_*.vcf files in a working directory into .csv format
for file in end_*.vcf; do mv -- "$file" "${file%.vcf}.csv"; done
