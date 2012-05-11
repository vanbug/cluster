# fastqc on all fastq files
for i in `ls *gz`; do fastqc $i ; done 
