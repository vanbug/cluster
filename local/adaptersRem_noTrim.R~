# Rscript for removing adapter from fastq files and using fastx_pipeline for each sample on qsub
# command : ~/src/scripts/adapterRem.R folder/*fq folder/adapter.txt 50
############################################################################################################

# initializing parameters
fastCommand<-c(); fastQsubEcho<-c() ; fastQsubFile<-c()

# reading args
args=commandArgs(TRUE)

# fetching fastq and adapter files
fastqFiles=args[1:(length(args)-1)]

# reading adapters per file and optional trim length
adapters=read.csv(tail(args,1),header=FALSE)

for (i in 1:length(fastqFiles)){

        # making fastx calling command
        fastCommand[i]=paste('sh ~/src/shell/fastx_collapser.sh',fastqFiles[i],adapters$V1[i],trimLength)

        # copying template qsub script for each file
        fastQsubFile[i]=paste('cp ~/src/qsub/template.qsub ',fastqFiles[i],'.tmp',sep='')
        system(fastQsubFile[i])

        # writing fastCommand to a temp file for later concatanation using cat to qsub script
        write(fastCommand[i],file=paste('tempa_',i,'.txt',sep=''))
        system(paste("cat",paste(fastqFiles[i],'.tmp',sep=''),paste('tempa_',i,'.txt',sep=''),">",paste(fastqFiles[i],'.qsub',sep='')))

        # removing tmp files
        system(paste("rm",paste('tempa_',i,'.txt',sep='')))
        system(paste("rm",paste(fastqFiles[i],'.tmp',sep='')))

        # executing qsub script
        system(paste("qsub ",fastqFiles[i],'.qsub',sep=''))
}
