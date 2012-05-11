#!/bin/sh
# USAGE [single sample run]      : ~/src/shell/macs14.sh file.bed
# USAGE2[sample and control run] : ~/src/shell/macs14.sh sample.bed control.bed

# running macs14 with and without control -- EXTEND IT WITH MACS2
if [ $# -ge 2 ]; then
	macs14 -t $1 -c $2 --bw 250 -g mm -n $1.macs14 2> $1_macs.log	
else
	macs14 -t $1 --bw 250 -g mm -n $1.macs14 2> $1_macs.log
fi

