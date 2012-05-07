# reading args
args=commandArgs(TRUE)
# fetching fastq and adapter files

fastqFiles=args[1]
system(paste('cp template.qsub ',fastqFiles[1],'.qsub',sep=''))
#system(com)
