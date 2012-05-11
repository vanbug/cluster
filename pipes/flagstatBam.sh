# Shell script for giving the mapping and other stats of the bam files
# Author : Sukhdeep Singh
# Organization : Max Planck
#############################################
#!bin/sh
echo "Some stats will be outputted such as reads mapped (forward and reverse strand) + unmapped"
for i in `ls *.bam`
do
echo "Detailing $i"
	samtools flagstat $i        
done
##########################################
