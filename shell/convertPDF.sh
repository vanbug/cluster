# converts pdf to jpg using imagemagic
for i in $1*pdf; do convert $i $i.jpg; done
