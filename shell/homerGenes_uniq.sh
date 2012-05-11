# filters unique transcripts randomly from the homer geneset, 16th column shows the geneNAME

sed 1d $1 | sort -u -k16 | sed 1d | sort -k2.5,2 -k3 >> $1_uniqueGenes.tmp
sed -n '1p' $1 > $1_header.tmp 
cat $1_header.tmp $1_uniqueGenes.tmp > $1_uniqueGenes.txt

# making a bed file for intersect bed
awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}' $1_uniqueGenes.txt > bedColumns.txt
cut -f16,18 $1_uniqueGenes.txt > geneNames.txt

paste -d" " bedColumns.txt geneNames.txt > $1_uniqueGenes.bed
echo "####################################################"
echo "$1_uniqueGenes.txt & $1_uniqueGenes.bed produced"

echo "Use .bed file for intersections"
echo "####################################################"

rm $1_header.tmp $1_uniqueGenes.tmp bedColumns.txt geneNames.txt
