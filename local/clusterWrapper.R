##########################################
# looping data per chromosome
chip=list();control<-list()

# import bs genome package
#library(BSgenome.Mmusculus.UCSC.mm9)
#library(chipseq)
library(ShortRead)

chrs=c(paste("chr",seq(1,19),sep=''),paste("chr",c('X','Y','M'),sep=''))
what=scanBamWhat()
number=seq(1,22,by=1)

for(i in 1:length(number)){
	which<-GRanges((chrs[number[i]]),ranges=IRanges(number[i],c(seqlengths(Mmusculus)[number[i]])))
	params<-ScanBamParam(which=which,what=what)
chip=readAligned('/projects/globalscratch/sukhi/rotterdam/pipe/uniqueRot/s_4_1.fastq.chipseq.bowtie.SR.mm9.sam.srt.bam.unique',type="BAM",param=params)
control=readAligned('/projects/globalscratch/sukhi/rotterdam/pipe/uniqueRot/s_5_1.fastq.chipseq.bowtie.SR.mm9.sam.srt.bam.unique',type="BAM",param=params)
print(paste("Working with",i))

# removing NA's
control=control[which(is.na(control@strand)==F)]
chip=chip[which(is.na(chip@strand)==F)]

# filtering data
filt<-occurrenceFilter(withSread=FALSE)
control<-control[filt(control)]
chip<-chip[filt(chip)]

# variable declarations
chip.norm.aln<-list(); control.norm.aln<-list();
chip.fragment.size<-c(); control.fragment.size<-c();
chip.cov<-c(); control.cov<-c();

# NORMALIZING LIBRARY SIZE FROM CHIP AND CONTROL - in GFP & MEN, control has more reads
sizeControl<-length(control@id)
sizeChip<-length(chip@id)
smallest.size<-min(sizeControl,sizeChip)

# normalizing chip and control library size
chip.norm.aln<-sample(chip,smallest.size,replace=FALSE)
control.norm.aln<-sample(control,smallest.size,replace=FALSE)

# fragment length extension
chip.fragment.size=round(estimate.mean.fraglen(chip.norm.aln,method="SISSR"))
control.fragment.size=round(estimate.mean.fraglen(control.norm.aln,method="SISSR"))

# Investigating Coverage
# calculating extend length
chip.extend.length=chip.fragment.size-width(chip.norm.aln)
control.extend.length=control.fragment.size-width(control.norm.aln)

# calculating coverage
chip.cov<-coverage(chip.norm.aln,width=seqlengths(Mmusculus),extend=as.integer(chip.extend.length))
control.cov<-coverage(control.norm.aln,width=seqlengths(Mmusculus),extend=as.integer(control.extend.length))

	# writing coverage externally
	bedGraph=function(cov,chrs,type){
		pf<-runValue(cov[[1]])>0
		write.table(cbind(chrs[i],start(cov[[1]])[pf],
		format(end(cov[[1]]),scientific=F)[pf],
		runValue(cov[[1]])[pf]),file=paste(type,".bedGraph",sep=''),append=TRUE,sep="\t",
		quote=FALSE,row.names=FALSE,col.names=FALSE)
	}
}
