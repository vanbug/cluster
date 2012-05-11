# reports the cpg intersections , store the output to have the list
intersectBed -u -a ~/src/software/homer/data/genomes/mm9/annotations/basic/cpgIsland.ann.txt_homerIntersections.bed -b $1 -f 0.5 | wc -l

