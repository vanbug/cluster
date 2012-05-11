# filters the top peaks from a macs peaks file (fdr threshold - bed format) - NOTE the sort is based on 7th column which has the pvalue.
if [ $# -ge 2 ]; then
	sort -rnk7 $1 | head -n $2 | sort -k1,1 -k2,2n > $1_$2.peaks.bed
else
	sort -rnk7 $1 | head -n 50 | sort -k1,1 -k2,2n > $1_top50.peaks.bed
fi
