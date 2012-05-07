# Shell Script producing bed files from unique (uniquely aligned bam files)
# NOTE : Run in the folder containing bed files
# USAGE : sh uni2bed.sh
#!/bin/bash


#file=$1
for i in `ls *unique`
 do
 # producing bed from $i 
 bamToBed -i $i > $i.bed

  # processing produced bed to remove unwanted columns
  cut -f1,2,3,6 $i.bed > $i.shortBed

  echo "Grepping out random chromosomes"
  egrep -w 'chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chrX|chrY' $i.shortBed > $i.removed
  rm $i.shortBed $i.bed
  echo "Done $i.removed produced for $i"
done
