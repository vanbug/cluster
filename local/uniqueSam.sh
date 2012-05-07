###################################################################################################################################
# Shell script for filtering the unique reads from bam file using samtools
# We should call it an reliable alignment than unique as 'uniqueness' is not well defined in general cases + this is done by setting a threshold on the mapping quality. [MAQ]
# Author : Sukhdeep Singh
# Organization : Max Planck Dresden
###################################################################################################################################

#!bin/sh
echo "This script will filter out the unique reads from the sorted BAM file"
for i in `ls *.bam`
do
echo "Producing unique read file from $i"
	samtools view -bq 1 $i > /projects/globalscratch/sukhi/beijing/data/mapping/bwa/uniqueBam/$i.unique
done
