# script to filter promoters for +/-2.5 KB
# removes header
sed 1d $1 > $1_fp.tmp
# 
echo "Promoter file generated" && intersectBed -u -a $1_fp.tmp -b ~/src/geneLists/homer_genesTSS_5000.txt_uniqueGenes_noHeader.bed -f 0.5 > $1_promoters.bed
echo "Number of promoters for +/-2.5 KB window" && wc -l $1_promoters.bed

# removes file
rm $1_fp.tmp
###########################################################
