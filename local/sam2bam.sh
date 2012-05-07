###################################################################################################################################
# Shell script for sam (human readable) to bam (machine readable)
# Author : Sukhdeep Singh
# Organization : Max Planck Dresden
###################################################################################################################################

#!bin/sh
echo "This script will convert all the sam files to bam files in the user defined folder in this script using samtools"
for i in `ls *.sam`
do
echo "Converting $i to bam (machine readable & compressed)"
        samtools view -bS -o /projects/globalscratch/sukhi/beijing/data/mapping/bwa/bam/$i.bam $i
done
