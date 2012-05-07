#!/bin/bash
# true coverage pipeline bam -> bw
#samtools sort L317_Track-327_R1.fastq.sam.bam.unique L317_Track-327_R1.fastq.sam.bam.unique.sorted &
file=$1

bamToBed -i $1 | awk '$6=="+"' | genomeCoverageBed -i stdin -bg -g ~/src/useFul/ucsc/genomeIndex/mm9.chrom > plus.bg
bamToBed -i $1 | awk '$6=="-"' | genomeCoverageBed -i stdin -bg -g ~/src/useFul/ucsc/genomeIndex/mm9.chrom > minus.bg
cat plus.bg | awk '{print $0"\t+"}'>plus_strand.bg && cat minus.bg | awk '{print $0"\t-"}'>minus_strand.bg && cat plus_strand.bg minus_strand.bg > full_strand.bg
sort -k1,1 -k2,2n full_strand.bg > full_strand_sorted.bg && awk '{OFS="\t"; print $1,$2,$3,NR,$4,$5}' full_strand_sorted.bg > wanted.bed

# extending reads
perl ~/src/chIPbedExtender.pl /biodata/biodb/ABG/genomes/mouse/mm9/mm9.fa.fai wanted.bed 250

# clipping
~/src/useFul/ucsc/utilities/bedClip wanted.extended.bed ~/src/useFul/ucsc/genomeIndex/mm9.chrom wanted.clipped

# merging reads
mergeBed -i wanted.clipped -n > wanted.merged

# bedGraph -> bigwig
~/src/useFul/ucsc/utilities/bedGraphToBigWig wanted.merged ~/src/useFul/ucsc/genomeIndex/mm9.chrom wanted.bw

# removing in pipe files
rm plus.bg minus.bg plus_strand.bg minus_strand.bg full_strand.bg full_strand_sorted.bg wanted.extended.bed wanted.clipped 

echo "Done"
