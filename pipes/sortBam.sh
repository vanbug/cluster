# Shell Script for sorting the bam files
#!bin/sh
# NOTE : Make sure the extension ends with .bam
# USAGE : sh sortBam.sh inputDir outputDir
###################################

echo "This script will sort all the bam files in the user defined folder using samtools"

for i in `ls $1/*.bam`
do
echo "Sorting $i"
        samtools sort $i $2/$i.sort
done

##################################
