# sort annotation file produced from HOMER
## .5 means focus on 5 characters, k means column
sed -e 1d $1 | sort -k2.5n -k3n > $1.sorted

