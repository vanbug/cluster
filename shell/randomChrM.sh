#!/bin/sh
# removes random chromosomes and chrM from the bed file along with the negative starts
egrep -v -e '-[[:digit:]]' $1 > $1.tmp
egrep -v -w '.*random|chrM' $1.tmp > $1.removed
rm $1.tmp
