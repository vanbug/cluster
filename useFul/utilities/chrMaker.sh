# Shell script for producing empty chr named text files
#!bin/sh
for ((i=1; i<20; i++));
do
	touch chr$i
done

touch chrX
touch chrY
touch chrM


##########################
