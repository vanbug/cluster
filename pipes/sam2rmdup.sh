#!bin/sh
# sampipe will all filtering enabled (sam->bam->unique->rmdup->bed)
# Author : Sukhdeep Singh

samtools view -bS $1 | samtools sort - $1.sorted && echo "Bam sorted" && samtools flagstat $1.sorted.bam > $1.sorted.txt && samtools view -bq 1 $1.sorted.bam > $1.unique && echo "unique reads filtered" && samtools flagstat $1.unique > $1.unique.txt && samtools rmdup -s $1.unique $1.rmdup "pcr duplicates removed" && samtools flagstat $1.rmdup > $1.rmdup.txt && bamToBed -i $1.rmdup | sortBed -i stdin > $1.bed && head $1.bed && perl ~/src/perl/chIPbedExtender.pl /biodata/biodb/ABG/genomes/mouse/mm9/mm9.fa.fai $1.bed 250
~/src/pipes/bedClip.sh $1.extended.bed
egrep -v -w '.*random|chrM' $1.extended.bed.clipped > $1.extended.bed.nonRandomChr.clipped
rm $1.sorted.bam $1.unique $1.bed $1.*rep* $1.extended.bed

echo "sam->bed done"

