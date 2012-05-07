print("Enter the coverage file")
covFile=scan("stdin",what="character",n=1,quiet=TRUE)
cov=read.table(covFile,stringsAsFactors=FALSE)
boo=cov[,3]==cov[,2]
cov[,3][boo]<-cov[,3][boo]+1

# we want data in non-scientific format, without 'e'
cov$V2=format(cov$V2,scientific=F)
cov$V3=format(cov$V3,scientific=F)
cov$V4=format(cov$V4,scientific=F)
write.table(cov,file=paste("covFile",".bedDiff",sep=""),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
print ("done")
