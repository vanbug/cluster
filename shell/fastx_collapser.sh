# Author : Sukhdeep Singh
# Date : 13/March/2012
# clips(removes the adapter sequences) the fastq sequences
# usage : sh ~/src/shell/fastx_collapser.sh fastq.gz ATTCAT 49 2> log
# $1 = fastq.gz file
# $2 = adapter sequence
# $3 = length to be trimmed down (OPTIONAL PARAMETER)

if [ -z "$1" ]; then
echo "You need at least one parameter"
fi


if [ $# -ge 3 ]; then
	gunzip -c $1 | fastx_trimmer -Q33 -l $3 -o $1.trimmed_$3
	fastx_clipper -Q33 -i $1.trimmed_$3 -a $2 -C -o $1.adapterClipped
	fastx_collapser -Q33 -v -i $1.adapterClipped -o $1.collapsed
	gzip $1.adapterClipped
	# switch this ON, when you get more space
	# gzip $1.trimmed_$3
	rm $1.trimmed_$3
else
	gunzip -c $1 | fastx_clipper -Q33 -a $2 -C -o $1.adapterClipped
	fastx_collapser -Q33 -v -i $1.adapterClipped -o $1.collapsed
	gzip $1.adapterClipped
	# switch this ON, when you get more space
	# gzip $1.trimmed_49
fi

