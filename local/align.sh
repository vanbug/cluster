###################################################################################################################################
# Shell script for aligning raw fasta files to reference genome using Bowtie
# Author : Sukhdeep Singh
# Organization : Max Planck Dresden
###################################################################################################################################

#!bin/sh 
for i in `ls *.clean.fq`
do
	echo "Aligning $i to reference genome MM9 using Bowtie"
	nice bowtie --best --strata -m 1 --seedlen 49 --sam /biodata/biodb/ABG/genomes/bowtie/mm9 /projects/globalscratch/sukhi/Beijing/data/fastq/$i /projects/globalscratch/sukhi/Beijing/data/mapping/runII/sam/$i.sam
	echo time
done
