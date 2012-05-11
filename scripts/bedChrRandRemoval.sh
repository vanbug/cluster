#!/bin/bash
# Author : SUkhdeep Singh
# Organization : Biotec
########################################################################################################################################################
# ATTENTION : PLEASE PUT ALL THE CLIPPED BED FILES TO BE PROCESSED INTO A NEW FOLDER & THEN RUN THIS SCRIPT AS IT PROCESS ALL THE FILES IN THE DIR PROVIDED AT INPUT BUT OUTPUTS THE FILES IN CURRENT DIR
# USAGE     : sh bedChrRandRemoval.sh *clipped
# USAGE 2   : sh bedChrRandRemoval.sh file.clipped
echo "Shell Script to remove random chromosomes from a directory full of BAM files or a single file";
# filename is random variable and samtools process the input bam file and thne piped to egrep
for i in `ls $1`
do
rand=$RANDOM

# removing unwanted 4,5 column from the bed file
egrep -w 'chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chrX|chrY' $1 > $1.nonRandomChr
echo "Done $i.nonRandomChr produced for $i"
done
#################################################################################################################################################################
