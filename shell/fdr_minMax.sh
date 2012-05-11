# calulates the min and maximum of FDR ratios from a Macs peak file
sed -e '1,24d' $1 | cut -f9 | awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1< min) {min=$1}; total+=$1; count+=1} END {print total/count, min, max}'
