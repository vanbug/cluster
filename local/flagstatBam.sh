#!bin/sh
echo "This script will output the flagstat values of the input BAM files - reads mapped (forward and reverse strand) + unmapped
for i in `ls *.bam`
do
echo "Detailing $i"
	samtools flagstat $i        
done

