#!/bin/sh
# removes random chromosomes and chrM from the bed file
egrep -v -w '.*random|chrM' $1 > $1.removed
