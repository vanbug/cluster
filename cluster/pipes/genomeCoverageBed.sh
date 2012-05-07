#!/bin/bash
# usage  : ~/src/pipes/genomeCoverageBed file.bed
# usage2 : ~/src/pipes/genomeCoverageBed file.bed --nobw
# calculates coverage from a bed file and converts it into bigwig

if [ $# -ge 2 ]; then
	genomeCoverageBed -bg -i $1 -g ~/src/useFul/ucsc/genomeIndex/mm9.chrom > $1.cov
else
	genomeCoverageBed -bg -i $1 -g ~/src/useFul/ucsc/genomeIndex/mm9.chrom > $1.cov
	~/src/pipes/bedGraphToBigwig.sh $1.cov
fi
echo "Done"
