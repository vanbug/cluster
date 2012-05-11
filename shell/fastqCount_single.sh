# gives you number of reads in the fastq file
echo -en $1"\t" && gunzip -c $1 | wc -l |  awk '{print $1/4}'
