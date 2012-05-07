# USAGE : ~/src/shell/du.sh .sai
# just input the common match and the code will search for the copies

# prints the total disk occupied by the inputted files
du -hs *$1 | awk '{total=total+$1} END {print total}'
