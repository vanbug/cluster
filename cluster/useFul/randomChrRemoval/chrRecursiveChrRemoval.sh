#!/bin/bash
# Author : SUkhdeep Singh
# Organization : Max Planck Soceity

# ATTENTION : PLEASE PUT ALL THE BAM FILES TO BE PROCESSED INTO A NEW FOLDER AND THEN RUN THIS SCRIPT AS IT PROCESS ALL THE FILES IN THE DIR PROVIDED AT INPUT BUT OUTPUTS THE FILES IN CURRENT DIRE
# USAGE     : sh chrRecursive.sh .
echo "Shell Script to remove random chromosomes from a directory full of BAM files";
# filename is random variable and samtools process the input bam file and thne piped to egrep

for i in `ls $1`
do
rand=$RANDOM
samtools view $i | egrep -w 'chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chrX|chrY|chrM' > "$i.$rand"
echo "Done $i.$rand produced for $i"
done
#################################################################################################################################################################
