# filters the top 50 peaks from a macs peaks file (bed format) - NOTE the sort is based on 5th column which has the enrichment ratio
sort -rnk5 $1 | head -n 50 > $1_top50.peaks.bed
