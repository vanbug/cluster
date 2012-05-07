#!/bin/bash
# Author : SUkhdeep Singh
# Organization : Max Planck Soceity
# ATTENTION : THE OUTPUT FILE WILL BE PRODUCED IN THE INPUT FILE DIRECTORY, NAME WILL BE OUPUTTED. MANAGE ACCORDINGLY
# USAGE	    : sh chrRemove.sh /bam/1.bam

echo "Shell Script to remove random chromosomes from a BAM file";

# filename is random variable and samtools process the input bam file and thne piped to egrep
rand=$RANDOM
file=$1
samtools view $1 | egrep -w 'chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chrX|chrY|chrM' > "$file.$rand"
echo "Done $file.$rand produced"

#################################################################################################################################################################
