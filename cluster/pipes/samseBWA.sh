# Shell script for the mapping raw fastq files to mm9 genome using BWA (single end reads)
# Author : Sukhdeep Singh
# Organization : Max Planck
# NOTE: pBWA is a good alternative for parallel processing (multithreading)
# USAGE : sh samseBWA.sh files/fastq/ files/sai files/output

##############################################################################

#!bin/sh
# getting raw fastq files with extension *.fq
f=(`cd $1;ls *.fq`)
files=(`ls $1/*.fq`)
# getting sai index directory
sai=(`cd $2;ls *.sai`)
size=${#files[@]}

# starting samse
for ((i=0; i<size; i++));
do
 	echo "Producing ${files[i]}.sam from ${files[i]}"
	nice bwa samse /biodata/biodb/ABG/genomes/mouse/mm9/mm9 ${sai[i]} ${files[i]} > $3${f[i]}.sam
done


###############################################################################
