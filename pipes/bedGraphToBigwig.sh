# date : 23/Feb/2010
# bedGraphToBigwig wrapper
# usage1 (for macs peaks) : ~/src/pipes/bedGraphToBigwig.sh peaks.bed -p
# usage2 (for normal bed): ~/src/pipes/bedGraphToBigwig.sh file.bed

#!bin/sh
if [ $# -ge 2 ]; then
	cut -f1,2,3,5 $1 > $1.cut
	~/src/useFul/ucsc/utilities/bedGraphToBigWig $1.cut ~/src/useFul/ucsc/genomeIndex/mm9.chrom $1.bw
	rm $1.cut
else
	~/src/useFul/ucsc/utilities/bedGraphToBigWig $1 ~/src/useFul/ucsc/genomeIndex/mm9.chrom $1.bw
fi
echo "Done"
