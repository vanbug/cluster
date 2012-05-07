###################################################################################################################################
# Shell script for producing aligned sam files from .sai index files using BWA's samse (sam for single end reads)
# Author : Sukhdeep Singh
# Organization : Max Planck Dresden
###################################################################################################################################

#!bin/sh
files=(`ls /projects/globalscratch/sukhi/beijing/data/fastq/*.clean.fq`)
sai=(`ls *.sai`)
size=${#files[@]}

# starting samse (sam production for single end reads)
for ((i=0; i<size; i++));
do
  	echo "Producing ${files[i]}.sam from ${sai[i]}"
        nice bwa samse /biodata/biodb/ABG/genomes/mouse/mm9/mm9 ${sai[i]} ${files[i]} > /projects/globalscratch/sukhi/beijing/data/mapping/bwa/sam/${sai[i]}.sam
        #echo ${files[i]}
        #echo ${sai[i]}
done
