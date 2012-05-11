# Date : 23/Feb/2012
# bedclip script
echo "Clipping input file"
~/src/useFul/ucsc/utilities/bedClip $1 ~/src/useFul/ucsc/genomeIndex/mm9.chrom $1.clipped
echo "Done"
