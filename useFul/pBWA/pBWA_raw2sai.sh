# pBWA script for aligning the raw fastq short read file to the reference genome of mouse (mm9)
# Author : Sukhdeep Singh
# Organization : Max Planck
######################################################################

#!/bin/bash
# qsub parameters beginning with #$
# Use bash as the interpreter
#$ -S /bin/bash
#before executing the command thange to the working directory
#$ -cwd
#submit the complete environment to the execution host
#$ -V
# merge stdout and stderr
#$ -j y
# set your output file
#$ -o output.txt

# now initialize the parallel environment with 20 cpus
#$ -pe orte 20
file=$1
#now we can run our command here.
/opt/openmpi/bin/mpirun -np 20 --mca btl_tcp_if_include eth0 pBWA aln -f $file.sai /biodata/biodb/ABG/genomes/mouse/mm9/mm9 $1

