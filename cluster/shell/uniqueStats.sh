# Author : Sukhdeep Singh
# Source : http://www.nature.com/nprot/journal/v7/n1/full/nprot.2011.420.html
# gives the over-represented sequences and the number of unique reads in a fastq file 
# usage : ~/src/shell/uniqueStats.sh file.fastq
# usage 2 : for i in `ls *.gz`; do ~/src/shell/uniqueStats.sh $i > uniqueStats.log ; done
#for i in `ls *.fastq`; do echo -en $i "\t";cat $i | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'; done
file=$1
echo -en $file "\t"; gunzip -c $file | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'

