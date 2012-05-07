#!/bin/bash

# extend first , then coverage
perl ~/src/chIPbedExtender.pl /biodata/biodb/ABG/genomes/mouse/mm9/mm9.fa.fai wanted.bed 250

# clipped the off-end boundaries
~/src/useFul/ucsc/utilities/bedClip wanted.extended.bed ~/src/useFul/ucsc/genomeIndex/mm9.chrom wanted.clipped

# get coverage
genomeCoverageBed -i wanted.clipped -bg -g ~/src/useFul/ucsc/genomeIndex/mm9.chrom > wanted.bg

# merge coverage
mergeBed -i wanted.bg -n > wanted.merged

# produce bigwig
~/src/useFul/ucsc/utilities/bedGraphToBigWig wanted.merged ~/src/useFul/ucsc/genomeIndex/mm9.chrom wanted.bw

#rm wanted.extended.bed wanted.clipped wanted.bg
echo "Done"
