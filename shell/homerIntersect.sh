# USAGE : ~/src/shell/homerIntersect.sh file.bed -head
# USAGE2 (without head) : ~/src/shell/homerIntersect.sh file.bed

if [ $# -ge 2 ]; then
	sed 1d $1 > $1_homerIntersections.tmp
	echo "+-2.5KB" && intersectBed -u -a $1_homerIntersections.tmp -b ~/src/geneLists/homer_genesTSS_5000.txt_uniqueGenes_noHeader.bed -f 0.5 | wc -l
	echo "+-5KB" && intersectBed -u -a $1_homerIntersections.tmp -b ~/src/geneLists/homer_genesTSS_10000.txt_uniqueGenes_noHeader.bed -f 0.5 | wc -l
	rm $1_homerIntersections.tmp
else
	# peaks intersections with custom homer (-2.5+2.5) & (-5+5) unique genelist
	echo "+-2.5KB" && intersectBed -u -a $1 -b ~/src/geneLists/homer_genesTSS_5000.txt_uniqueGenes_noHeader.bed -f 0.5 | wc -l
	echo "+-5KB" && intersectBed -u -a $1 -b ~/src/geneLists/homer_genesTSS_10000.txt_uniqueGenes_noHeader.bed -f 0.5 | wc -l
fi
