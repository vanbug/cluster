# this script adds a timestamp to the file using awk
cat $1 | awk '{ print strftime("#%Y-%m-%d %H:%M:%S\n"),$0; }' > $1.TXT
