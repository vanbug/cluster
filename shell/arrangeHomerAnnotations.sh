# rearranges homer annotation file to be used with intersectBed for intersections

cut -f2,3,4,5 $1 > $1_homerArrange_A.tmp
cut -f1,6 $1 >  $1_homerArrange_B.tmp
paste -d" " $1_homerArrange_A.tmp $1_homerArrange_B.tmp > $1_homerIntersections.bed
rm $1_homerArrange_A.tmp $1_homerArrange_B.tmp

