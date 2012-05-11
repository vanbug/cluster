# usage ~/src/bed2macsBed.sh file.bed
awk '{OFS="\t"; print $1,$2,$3,NR,$7}' $1 > $1.macsBed
