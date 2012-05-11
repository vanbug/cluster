# reading file from command Args
homerFile=commandArgs(TRUE)[1]

# invoke system command to sort the homer file
system(paste("sort -k2.5,2 -k3",homerFile,"> homerFile_R.tmp"))

# adding na to empty places
genes=read.csv('homerFile_R.tmp',sep="\t",na.strings="",header=TRUE)

# removing na's - gene transcripts with empty gene names but names ucsc transcripts
genes=genes[which(is.na(genes$Gene.Name)==FALSE),]

# ordering and sorting list
#geneOrder=mixedorder(genes$Gene.Name)
#genesSorted=genes[geneOrder,]

# removing duplicated ones - actually taking the isoform with longest length [list is sorted with for ]
uniqGenes=genes[!duplicated(genes$Gene.Name),]

system('rm homerFile_R.tmp')

write.table(uniqGenes,paste(homerFile,".R_uniqSorted",sep=''),sep='\t',row.names=FALSE,quote=FALSE,col.names=FALSE)

