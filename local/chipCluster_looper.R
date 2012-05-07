# looper clone of chipCluster.R for processing each chromosome
# Author : Sukhdeep Singh
# Organization : Biotec, Dresden
#######################################################################
chip=list();control<-list()
# import bs genome package
library(BSgenome.Mmusculus.UCSC.mm9)
library(chipseq)
chrs=c(paste("chr",seq(1,19),sep=''),paste("chr",c('X','Y','M'),sep=''))
what=scanBamWhat()
number=seq(1,22,by=1)
for(i in 1:length(number)){
	which<-GRanges((chrs[number[i]]),ranges=IRanges(number[i],c(seqlengths(Mmusculus)[number[i]])))
params<-ScanBamParam(which=which,what=what)
chip[[i]]=readAligned('/projects/globalscratch/sukhi/rotterdam/pipe/uniqueRot/s_4_1.fastq.chipseq.bowtie.SR.mm9.sam.srt.bam.unique',type="BAM",param=params)
control[[i]]=readAligned('/projects/globalscratch/sukhi/rotterdam/pipe/uniqueRot/s_5_1.fastq.chipseq.bowtie.SR.mm9.sam.srt.bam.unique',type="BAM",param=params)
print(paste("Chromosome",i,"imported"))
}
# variable declarations
chipP<-list(); controlP<-list(); controlF<-list(); chipF<-list()

# uniquely aligned data, dropping NA's
for (i in 1:length(chip)){
controlP[[i]]=control[[i]][which(is.na(control[[i]]@strand)==F)]
chipP[[i]]=chip[[i]][which(is.na(chip[[i]]@strand)==F)]
#chip2P[[i]]=chip2p[[i]][which(is.na(chip[[i]]@strand)==F)]
print(paste("chr",i,sep=''))
}


# inbuilt filtering for alignQuality occurrenceFilter
filt1<-chromosomeFilter("chr[0-9XYM]")
filt2<-alignQualityFilter(5) # random quality filter selected on quality scores by bwa
filt3<-occurrenceFilter(withSread=FALSE)
filt<-compose(filt1,filt2,filt3)

# applying filters on chip and control
for (i in 1:length(chip)){
controlF[[i]]<-controlP[[i]][filt(controlP[[i]])]
chipF[[i]]<-chipP[[i]][filt(chipP[[i]])]
#chip2F[[i]]<-chip2P[[i]][filt(chipP[[i]])]
print(paste("chr",i,sep=''))
}

# variable declarations
chip.norm.aln<-list(); control.norm.aln<-list(); chip2.norm.aln<-list()
chip.fragment.size<-c(); chip2.fragment.size<-c(); control.fragment.size<-c()
chip.cov<-c(); chip2.cov<-c(); control.cov<-c()

# NORMALIZING LIBRARY SIZE FROM CHIP AND CONTROL
for (i in 1:length(chipF)){
sizeControl<-length(controlF[[i]]@id)
sizeChip<-length(chipF[[i]]@id)
#sizeChip2<-length(chip2Chrs[[i]]@id)
smallest.size<-min(sizeControl,sizeChip)
#smallest.size2<-min(sizeControl,sizeChip2)
chip.norm.aln[[i]]<-sample(chipF[[i]],smallest.size,replace=FALSE)
#chip2.norm.aln[[i]]<-sample(chip2Chrs[[i]],smallest.size2,replace=FALSE)
control.norm.aln[[i]]<-sample(controlF[[i]],smallest.size,replace=FALSE)
print (paste("Done with Chr",i,sep=''))
#print (chip.norm.aln[[i]])
#print (control.norm.aln[[i]])
chip.fragment.size=c(chip.fragment.size,round(estimate.mean.fraglen(chip.norm.aln[[i]],method="SISSR")))
#chip2.fragment.size=c(chip2.fragment.size,round(estimate.mean.fraglen(chip2.norm.aln[[i]],method="SISSR")))
control.fragment.size=c(control.fragment.size,round(estimate.mean.fraglen(control.norm.aln[[i]],method="SISSR")))
}

# Investigating Coverage
# fetching the sequence lengths of mm9 - the chromosome size
library(BSgenome.Mmusculus.UCSC.mm9)
seqlengths(Mmusculus)

# calculating extend length
chip.extend.length=chip.fragment.size-width(chip.norm.aln[[1]][1])
#chip2.extend.length=chip2.fragment.size-width(chip2.norm.aln[[1]][1])
control.extend.length=control.fragment.size-width(control.norm.aln[[1]][1])


# calculating coverage
for (i in 1:length(chipF)){
chip.cov<-c(chip.cov,coverage(chip.norm.aln[[i]],width=seqlengths(Mmusculus),extend=as.integer(chip.extend.length[i])))
#chip2.cov<-c(chip2.cov,coverage(chip.norm.aln[[i]],width=seqlengths(Mmusculus),extend=as.integer(chip2.extend.length[i])))
control.cov<-c(control.cov,coverage(control.norm.aln[[i]],width=seqlengths(Mmusculus),extend=as.integer(control.extend.length[i])))
print (paste((length(chipF)-i),"chromosomes left",sep=''))
}

# for cov.no.zero , score is present at chip.cov.zero@values@unlistData@listData$score
chrs=c(paste('chr',seq(1,19),sep=''),'chrX','chrY','chrM')

# writing bedGraphs using write.table for all chromosomes, x= chip.cov, type = chip
bedGraph=function(cov,chrs,type){
for (i in 1:length(cov)){
	pf<-runValue(cov[[i]][[1]])>0
	write.table(cbind(chrs[i],start(cov[[i]][[1]])[pf],
	format(end(cov[[i]][[1]]),scientific=F)[pf],
	runValue(cov[[i]][[1]])[pf]),file=paste(type,".bedGraph",sep=''),append=TRUE,sep="\t",
	quote=FALSE,row.names=FALSE,col.names=FALSE)
}
}
##############################################################################
# optional plotting of coverage
# JUST WITH INITIAL ANALYSIS - NO NEED OF RECURRANCE
# plotting function for coverage visualization of each chromosome
plotChIP.Coverage<-function(x,xlab="Position",ylab="Coverage",main="ChIP Coverage",sub){
	plot(c(start(x),length(x)),c(runValue(x),1),type="l",col="blue",xlab=xlab,ylab=ylab,main=main,sub=sub)
}

plotControl.Coverage<-function(x,xlab="Position",ylab="Coverage",main="Control Coverage",sub){
	plot(c(start(x),length(x)),c(runValue(x),1),type="l",col="red",xlab=xlab,ylab=ylab,main=main,sub=sub)
}

# plots coverage for control and chip single or multiple chromosomes, x=chip.cov
plotChr=function(x,type){
	for (i in 1:length(x)) {
	if (type=="control") {jpeg(paste("chr",i,".jpg",sep=''))
				plotControl.Coverage(x[[i]][[1]],sub=paste(type,"chr",i))} 
	if (type=="chip")    {jpeg(paste("chr",i,".jpg",sep=''))
				plotChIP.Coverage(x[[i]][[1]],sub=paste(type,"chr",i))}
	if (type=="chip2")   {jpeg(paste("chr",i,".jpg",sep='')) 
				plotChIP.Coverage(x[[i]][[1]],sub=paste(type,"chr",i))}
	dev.off()
	}
}

## plotting chrs
plotChr(chip.cov,"chip")
plotChr(control.cov,"control")
###############################################################
# finding islands - regions of continuous reads >5 which represent the real binding sites


