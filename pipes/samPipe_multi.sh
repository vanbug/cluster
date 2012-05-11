# Samtools pipe for converting SAM->BAM->sorted->unique
# Author : Sukhdeep Singh
# Organisation : Max Planck
#########################################################################

files=(`cd $1;ls *.sam`)
#size=${#files[@]}
mkdir samPipe
for i in `ls $1/*.sam`
do
echo "Sorting $i"
	samtools view -bS $i -o samPipe/$i.bam
        samtools sort samPipe/$i.bam samPipe/$i.bam.sort
        samtools view -bq 1 samPipe/$i.bam.sort.bam > samPipe/$i.bam.sort.unique
        samtools index samPipe/$i.bam.sort.unique

done

#for ;
#do
#echo "Processing $1"
 #       samtools view -bS $1 -o $1.bam
  #      samtools sort $1.bam $1.bam.sort
   #     samtools view -bq 1 $1.bam.sort.bam > $1.bam.sort.unique
#	samtools index $1.bam.sort.unique
#done

