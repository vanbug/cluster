# script to print total reads, unique reads, percentage of unique reads, most repeated sequence, no. of times its repeated, percentage of ots repeatedness
# use 'fastx_collapser -v -i chip_dmel.fastq -o chip_dmel.fastq.collapsed' to know the number of frequency of repeated sequences in the whole dataset
cat $1 | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'
