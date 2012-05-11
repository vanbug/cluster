# Shell script for filtering the unique reads from bam file
## lets call it an reliable alignment than unique as 'uniqueness' is not well defined in general cases + this is done by setting a threshold on the mapping quality. [MAQ]
# Author : Sukhdeep Singh
# Organization : Max Planck
# USAGE : sh uniqueBam.sh 
####################################################################################################################################

#!bin/sh
echo "This script will filter out the unique reads from the sorted BAM file"
files=(`cd $1;ls *.bam`)
size=${#files[@]}

for ((i=0;i<size;i++));
do
echo "Producing unique read file from ${files[i]}"
	samtools view -bq 1 ${files[i]} > $2${files[i]}.unique

done
####################################################################################################################################
