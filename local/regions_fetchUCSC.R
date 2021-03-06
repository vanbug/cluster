# Rscript to produce a regions list to be called by fetch_ucsc to make snapshot fo UCSC with custom tracks

# read top genes
genesFile=commandArgs(TRUE)[1]
genes=read.csv(genesFile,sep='\t',header=FALSE)

# filter unique non-repetitive genes (as two peaks can be highly significant and binding to the same gene)
uniqueGenes=genes[which(duplicated(genes$V11)==FALSE),]

# print the number of unique genes removed
print(paste(length(which(duplicated(genes$V11)==TRUE)),"genes were identical, meaning there were multiple highly significant peaks within the same genes at different postions, gene might be big"))

# export this list an R formatted csv table in accordance to fetch_ucsc.py script [https://bitbucket.org/dalloliogm/ucsc-fetch/src/ef07f97a35a8/fetch_ucsc.py]
# added and removed 200 from the start and the end of the gene co-ordinates so as get a little zoomed snapshot
genesExport=cbind(as.character(uniqueGenes$V11),rep("mouse",length(uniqueGenes$V1)),rep("mm9",length(uniqueGenes$V1)),as.character(uniqueGenes$V6),(uniqueGenes$V7-4000),(uniqueGenes$V8+4000),rep("NULL",length(uniqueGenes$V1)),rep(0,length(uniqueGenes$V1)),rep(0,length(uniqueGenes$V1)))

# write the table in csv format
write.table(genesExport,"regionsFile_tmp",sep=',',quote=FALSE,row.names=FALSE,col.names=FALSE)

# we need to add header to the this exported file
system ('echo "#label,organism,chromosome,start,end,description,upstream,downstream" > regions_head')
system ('cat regions_head regionsFile_tmp > regionsFile')
system ('rm regionsFile_tmp')
