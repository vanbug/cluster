# gives you number of reads in the fastq file
for i in `ls *gz`; do echo -en $i"\t" && gunzip -c $i | wc -l | awk '{print $1/4}'; done > count.log
