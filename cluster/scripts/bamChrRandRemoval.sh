#!/bin/bash
# Author : SUkhdeep Singh
# Organization : Max Planck Soceity

# ATTENTION : PLEASE PUT ALL THE BAM FILES TO BE PROCESSED INTO A NEW FOLDER AND THEN RUN THIS SCRIPT AS IT PROCESS ALL THE FILES IN THE DIR PROVIDED AT INPUT BUT OUTPUTS THE FILES IN CURRENT DIRE
# USAGE     : sh chrRecursive.sh .
echo "Shell Script to remove random chromosomes from a directory full of BAM files or a single file";
# filename is random variable and samtools process the input bam file and thne piped to egrep
file=$1
for i in `ls $1`
do
rand=$RANDOM

# removing unwanted 4,5 column from the bed file
bamToBed -i $1 > $file.bed
cut -f1,2,3,6 $file.bed > $file.shortBed
egrep -w 'chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chrX|chrY' $file.shortBed > $file.BED
rm $file.shortBed $file.bed
echo "Done $file.removed produced for $i"
done
#################################################################################################################################################################
