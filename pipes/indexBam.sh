# shell script for indexing the bam files

#!bin/sh
echo "This script will index the bam files uning the samtools index and 
produces .sai files "
for i in `ls *.unique`
do
echo "Producing index of $i"
	samtools index $i
done


