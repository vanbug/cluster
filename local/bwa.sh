###################################################################################################################################
# Shell script for aligning raw fasta files to reference genome using BWA producing .sai index files
# Author : Sukhdeep Singh
# Organization : Max Planck Dresden
###################################################################################################################################

#!bin/sh
for i in `ls *.clean.fq`
do
	echo "Aligning $i to reference genome MM9 using BWA"
	nice bwa aln -t 12 /biodata/biodb/ABG/genomes/mouse/mm9/mm9 $i > /projects/globalscratch/sukhi/beijing/data/mapping/bwa/$i.sai
	echo time
done

