###################################################################################################################################
# Shell script for sorting BAM files
# Author : Sukhdeep Singh
# Organization : Max Planck Dresden
###################################################################################################################################

#!bin/sh
echo "This script will sort all the bam files in the user defined folder 
using samtools"
for i in `ls *.bam`
do
echo "Sorting $i"
        samtools sort $i /projects/globalscratch/sukhi/beijing/data/mapping/bwa/sortedBam/$i.sort
done
