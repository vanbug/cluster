# Shell Script for converting the SAM(Human Readable) to BAM (Machine Readable)
# Author : Sukhdeep Singh
# Organization : Max Planck
# Benefits : Encypted, Small Size.
# Reverting Back : samtools view file.bam > file.sam
# Usage : sh sam2bam.sh files/sam/ files/bam
#########################################################################################################################

#!/bin/sh 
#echo "This script will convert all SAM files in the user inputted folder to BAM in the user ouputted folder"
files=(`cd $1;ls *.sam`)
#outputDir=(`cd $2;pwd`)
size=${#files[@]}

for ((i=0;i<size;i++));
do
echo "Converting ${files[i]} to ${files[i]}.bam"
	samtools view -bS ${files[i]} -o $2${files[i]}.bam
done
########################################################################################################################
