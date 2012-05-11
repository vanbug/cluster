# Shell script for producing bed file from BAM files
# NOTE : The output bed file will be produced in the same folder as the input file. It also requires the mm9 genome index which can be downloaded from UCSC FTP

#!/bin/sh
bamFile=$1
genomeCoverageBed -bg -ibam $1 -g ../ucsc/genomeIndex/mm9.chrom > "$bamFile.bed"
echo "$bamFile.bed is the bes file for input $1"


ls 
