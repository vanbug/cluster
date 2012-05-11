# overlap summary
#!bin/sh
# usage : ~/src/shell/overlaps.sh A.bed B.bed

if [ -z "$1" ]; then
echo "~/src/shell/overlaps.sh A.bed B.bed"
exit 1
fi

echo -en "# Command used : sh ~/src/shell/overlaps.sh $1 $2\n" >> overlaps.log
echo -en "intersectBed -u -a $1 -b $2\n" >> overlaps.log
sed -e '1d' $1 > a.tmp
sed -e '1d' $2 > b.tmp
# biological factor names
#A=$3
#B=$4
# calculations
f1A=`intersectBed -u -a a.tmp -b b.tmp | wc -l`
f50A=`intersectBed -u -a a.tmp -b b.tmp -f 0.5 | wc -l`
f90A=`intersectBed -u -a a.tmp -b b.tmp -f 0.9 | wc -l`
f100A=`intersectBed -u -a a.tmp -b b.tmp -f 1.0 | wc -l`

# writing overlaps as bed file for 50% overlaps
intersectBed -u -a a.tmp -b b.tmp -f 0.5 > A_B_overlap.bed
# printing
echo -en "f0.1\t"$f1A"\n" >> overlaps.log
echo -en "f0.5\t"$f50A"\n" >> overlaps.log
echo -en "f0.9\t"$f90A"\n" >> overlaps.log
echo -en "f1.0\t"$f100A"\n" >> overlaps.log
echo -en "-------------------------------------------------------------------------------------\n" >> overlaps.log
# calculations
f1B=`intersectBed -u -a b.tmp -b a.tmp | wc -l`
f50B=`intersectBed -u -a b.tmp -b a.tmp -f 0.5 | wc -l`
f90B=`intersectBed -u -a b.tmp -b a.tmp -f 0.9 | wc -l`
f100B=`intersectBed -u -a b.tmp -b a.tmp -f 1.0 | wc -l`

# writing overlaps as bed file
intersectBed -u -a b.tmp -b a.tmp -f 0.5 > B_A_overlap.bed
totalA=`cat a.tmp | wc -l`
totalB=`cat b.tmp | wc -l`
# printing
echo -en "intersectBed -u -a $2 -b $1\n" >> overlaps.log
echo -en "f0.1\t"$f1B"\n" >> overlaps.log
echo -en "f0.5\t"$f50B"\n" >> overlaps.log
echo -en "f0.9\t"$f90B"\n" >> overlaps.log
echo -en "f1.0\t"$f100B"\n">> overlaps.log
echo -en "-------------------------------------------------------------------------------------\n" >> overlaps.log
echo -en "Total peaks for "$1" =" $totalA"\n" >> overlaps.log
echo -en "Total peaks for "$2" =" $totalB"\n" >> overlaps.log

echo -en $totalA"\t"$totalB"\t"$f1A"\n"  >> overlap_venn.R
echo -en $totalA"\t"$totalB"\t"$f50A"\n" >> overlap_venn.R
echo -en $totalA"\t"$totalB"\t"$f90A"\n" >> overlap_venn.R
echo -en $totalA"\t"$totalB"\t"$f100A"\n">> overlap_venn.R
echo -en $totalB"\t"$totalA"\t"$f1B"\n"  >> overlap_venn.R
echo -en $totalB"\t"$totalA"\t"$f50B"\n" >> overlap_venn.R
echo -en $totalB"\t"$totalA"\t"$f90B"\n" >> overlap_venn.R
echo -en $totalB"\t"$totalA"\t"$f100B"\n">> overlap_venn.R

rm a.tmp b.tmp
########################################### FINISH ############################################################
