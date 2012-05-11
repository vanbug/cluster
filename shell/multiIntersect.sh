# multiIntersectBed for the three files - WORKS FOR THREE FILES ONLY
if [ $# -ge 5 ]; then
	sed -e '1d' $1 > $1_MIB.tmp
	sed -e '1d' $2 > $2_MIB.tmp
	sed -e '1d' $3 > $3_MIB.tmp
	multiIntersectBed -i $1_MIB.tmp $2_MIB.tmp $3_MIB.tmp | awk '$4==3' > A_B_C.MIB
	wc -l multiIntersect_$1
	rm $1_MIB.tmp $2_MIB.tmp $3_MIB.tmp
else 
	multiIntersectBed -i $1 $2 $3 | awk '$4==3' > A_B_C.MIB
fi
