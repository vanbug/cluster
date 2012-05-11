#  Rscript for the bed2cvg
# Input chip and control file
chipFile=commandArgs()[5]
controlFile=commandArgs()[6]

# intializing libraries
library(BSgenome.Mmusculus.UCSC.mm9);library(ShortRead);library(chipseq);library(ff)
#data<-read.table(file,skip=1)
chip=read.table.ffdf(file=chipFile)
print ("Read Chip")
control=read.table.ffdf(file=controlFile)
print ("Read Control")

# setting parameters
fragment.size=200;chr.size<-seqlengths(Mmusculus)

# SMALL HACK FOR WORKING WITH ffdf
chip.subset=chip
control.subset=control

#GRanges, used a hack to work with ffdf table
chip.GRanges<-GRanges(seqnames=as.character(chip.subset$V1[1:length(chip.subset$V1)]),IRanges(start=chip.subset$V2[1:length(chip.subset$V1)],
end=chip.subset$V3[1:length(chip.subset$V1)]),strand=chip.subset$V4[1:length(chip.subset$V1)])
print ("Chip Granged")
control.GRanges<-GRanges(seqnames=as.character(control.subset$V1[1:length(control.subset$V1)]),IRanges(start=control.subset$V2[1:length(control.subset$V1)],end=control.subset$V3[1:length(control.subset$V1)]),strand=control.subset$V4[1:length(control.subset$V1)])
print ("Control Granged")
seqlengths(chip.GRanges)<-chr.size[c(1,10:19,2:9,20:21)]
seqlengths(control.GRanges)<-chr.size[c(1,10:19,2:9,20:21)]
chip.GRanges <- resize(chip.GRanges,width=fragment.size)
control.GRanges <- resize(control.GRanges,width=fragment.size)
print ("resized")

save(chip.GRanges,file=paste(chipFile,'GRanges',sep=''))
#save(control.GRanges,file=paste(controlFile,'GRanges',sep=''))

#Removing duplicates for the same start positions and same strand orientation
# cleaning has no big effect, just the duplicates explained above are removed
chip.clean<-chip.GRanges[!duplicated(paste(start(chip.GRanges),strand(chip.GRanges),sep=":"))]
control.clean<-control.GRanges[!duplicated(paste(start(control.GRanges),strand(control.GRanges),sep=":"))]
print ("Data Cleaned")

# NORMALIZING LIBRARY SIZE FROM CHIP AND CONTROL using sampling
sizeChip<-length(chip.clean)
sizeControl<-length(control.clean)
smallest.size<-min(sizeControl,sizeChip)

print ("Sampling Datasets")

# normalizing chip and control library size
chip.norm.aln<-sample(chip.clean,smallest.size,replace=FALSE)
control.norm.aln<-sample(control.clean,smallest.size,replace=FALSE)

# fragment length extension using  #SISSR- Site Identification from Short Sequence Reads
#chip.fragment.size=round(estimate.mean.fraglen(chip.norm.aln,method="SISSR"))
#control.fragment.size=round(estimate.mean.fraglen(control.norm.aln,method="SISSR"))

print ("Obtaining Coverage")

# Coverage - Extending read length like this doesn't gives peaks for some reason
# calculating extend length - length of the read to be extended from 36 bp to 200 or so, to which 
# it was broken down in start
#chip.extend.length=chip.fragment.size-mean(unique(width(chip.norm.aln)))
#control.extend.length=control.fragment.size-mean(unique(width(control.norm.aln)))

# Coverage 
# NOTE : WE ALREADY EXTENDED THE FRAGMENTS SO THERE IS NO NEED TO DO IT AGAIN
# calculating extend length - length of the read to be extended from 36 bp to 200 or so, to which it was broken down in start
#chip.extend.length=chip.fragment.size-width(chip.norm.aln)
#control.extend.length=control.fragment.size-width(control.norm.aln)

# calculating coverage, extend - extends the read reference from 36 to 200
#chip.cov<-coverage(chip.norm.aln,width=seqlengths(Mmusculus),extend=as.integer(chip.extend.length))
#control.cov<-coverage(control.norm.aln,width=seqlengths(Mmusculus),extend=as.integer(control.extend.length))
print ("Starting Coverage")
chip.cov<-coverage(chip.norm.aln,width=seqlengths(Mmusculus))
control.cov<-coverage(control.norm.aln,width=seqlengths(Mmusculus))

#Coverage
#pdf('data.pdf')
#plot(start(data.cov[[1]]),runValue(data.cov[[1]]),type='l')
#dev.off()

print ("Writing Coverage")

# writing coverage externally
bedGraph=function(cov,chrs,type){
        pf<-runValue(cov)>0
        write.table(cbind(chrs,start(cov)[pf],
        format(end(cov),scientific=F)[pf],        
	runValue(cov)[pf]),file=paste(type,".bedGraph",sep=''),append=TRUE,sep="\t",
        quote=FALSE,row.names=FALSE,col.names=FALSE)
        }
############## CONTROL IS OFF - SWITCH IT ON TO PRODUCE CONTROL BED FILE ##############################
for (chrs in c(paste('chr',seq(1,19),sep=''),'chrX','chrY')){
bedGraph(chip.cov[[chrs]],chrs,type=chipFile)
#bedGraph(control.cov[[chrs]],chrs,type=controlFile)
print (paste("Processing",chrs))
}


# writing a log file
write(paste(paste("Chip Library Size:",sizeChip),paste("Control Library Size:",sizeControl),sep='\n'),file=paste(chipFile,'.log',sep=''))
