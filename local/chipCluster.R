###################################################################################################################################
# R script for analysis of Chip Seq data after Shell pre-processingaligned data
# Author : Sukhdeep Singh
# Organization : Max Planck Dresden
###################################################################################################################################

# importing libraries
library('ShortRead');
library('chipseq'); # required for estimate.mean.fraglen

# declaring variables
chrs<-c();chipChr<-list();controlChr<-list();

# reading the user inputted file
control<-readline('Please enter the control location! ')
control<-gsub('\'','',control)
chip<-readline('Please enter the chip location! ')
chip<-gsub('\'','',chip)
chip2<-readline('Please enter the chip II location! ')
chip2<-gsub('\'','',chip2)



# reading control and chip data
control<-readAligned(control,type="BAM");
chip<-readAligned(chip,type="BAM");
chip2<-readAligned(chip2,type="BAM");

# printing info
print (control);
print (chip);
print (chip2);

# uniquely aligned data, dropping NA's
controlP=control[which(is.na(control@strand)==F)]
chipP=chip[which(is.na(chip@strand)==F)]
chip2P=chip2[which(is.na(chip2@strand)==F)]

# inbuilt filtering for alignQuality occurrenceFilter
filt1<-chromosomeFilter("chr[0-9XYM]")
filt2<-alignQualityFilter(25) # random quality filter selected on quality scores by bwa
filt3<-occurrenceFilter(withSread=FALSE)
filt<-compose(filt1,filt2,filt3)
controlF<-controlP[filt(controlP)]
chipF<-chipP[filt(chipP)]
chip2F<-chip2P[filt(chip2P)]

# sort chromosome levels 
sortChr=function(x){
inputlevel=levels(x)
sortedchrlevel=paste("chr",sort(as.numeric(gsub("chr","",as.character(inputlevel)))),sep="")
diff=length(inputlevel)-length(sortedchrlevel)
if (length(inputlevel)>length(sortedchrlevel)&&(diff!=3)) {print ("Investigate")} 
else {sortedchrlevel[(length(sortedchrlevel)+1)]="chrX";sortedchrlevel[(length(sortedchrlevel)+1)]="chrY";sortedchrlevel[(length(sortedchrlevel)+1)]="chrM"}
sortedchrlevel=paste(sortedchrlevel,"$",sep='')
#sortedchrlevel[(length(sortedchrlevel)+1):(length(sortedchrlevel)+3)]=c("chrX","chrY","chrM")
return (sortedchrlevel)
}

# sort chromosome levels
controlLevel=sortChr(controlF@chromosome);
chipLevel=sortChr(chipF@chromosome);
chip2Level=sortChr(chip2F@chromosome)

# simplified chromosome extraction - single chromosome, multiple sample extraction
extractChr=function(control1level=NULL,control2level=NULL,chip1level=NULL,chip2level=NULL,chip3level=NULL,control1=NULL,control2=NULL,chip1=NULL,chip2=NULL,chip3=NULL){
if (is.null(control1)==FALSE) {controlChr1=control1[(chromosomeFilter(control1level))(control1)]} else {controlChr1=NULL}
if (is.null(control2)==FALSE) {controlChr2=control2[(chromosomeFilter(control2level))(control2)]} else {controlChr2=NULL}
if (is.null(chip1)==FALSE)    {chipChr1=chip1[(chromosomeFilter(chip1level))(chip1)]} else {chipChr1=NULL}
if (is.null(chip2)==FALSE) {chip2Chr=chip2[(chromosomeFilter(chip2level))(chip2)]} else {chip2Chr=NULL}
if (is.null(chip3)==FALSE) {chip3Chr=chip3[(chromosomeFilter(chip3level))(chip3)]} else {chip3Chr=NULL}
return (list(control1=controlChr1,control2=controlChr2,chip1=chipChr1,chip2=chip2Chr,chip3=chip3Chr))
}

# lappy usage - multiple chromosome, multiple sample extraction
extractNucleosome=function(control1level=NULL,chip1level=NULL,chip2level=NULL,control1=NULL,chip1=NULL,chip2=NULL){
chrs<-mapply(extractChr,control1level=control1level,chip1level=chip1level,chip2level=chip2level,MoreArgs=list(control1=control1,chip1=chip1,chip2=chip2))
return (chrs)
}

# extract chip's n control's chromosome specific data
extractedChrs<-extractNucleosome(controlLevel,chipLevel,chip2Level,control1=controlF,chip1=chipF,chip2=chip2F)

# for only chip and control
#extractedChrs<-extractNucleosome(controlLevel,chipLevel,control1=controlF,chip1=chipF)

#fetching chromosome data using seq
controlChrs=extractedChrs[seq(1,length(extractedChrs),by=5)]
chipChrs=extractedChrs[seq(3,length(extractedChrs),by=5)]
chip2Chrs=extractedChrs[seq(4,length(extractedChrs),by=5)]

# unlisting chip and control chromosome data to remove head list - NO NEED WITH CURRENT CODE
controlChrs=unlist(controlChrs)
chipChrs=unlist(chipChrs)
chip2Chrs=unlist(chip2Chrs)

# variable declarations
chip.norm.aln<-list(); control.norm.aln<-list(); chip2.norm.aln<-list()
chip.fragment.size<-c(); chip2.fragment.size<-c(); control.fragment.size<-c()
chip.cov<-c(); chip2.cov<-list(); control.cov<-c()

# NORMALIZING LIBRARY SIZE FROM CHIP AND CONTROL - in GFP & MEN, control has more reads
for (i in 1:length(chipChrs)){
sizeControl<-length(controlChrs[[i]]@id)
sizeChip<-length(chipChrs[[i]]@id)
#sizeChip2<-length(chip2Chrs[[i]]@id)
smallest.size<-min(sizeControl,sizeChip)
#smallest.size2<-min(sizeControl,sizeChip2)
chip.norm.aln[[i]]<-sample(chipChrs[[i]],smallest.size,replace=FALSE)
#chip2.norm.aln[[i]]<-sample(chip2Chrs[[i]],smallest.size2,replace=FALSE)
control.norm.aln[[i]]<-sample(controlChrs[[i]],smallest.size,replace=FALSE)
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
chip2.extend.length=chip2.fragment.size-width(chip2.norm.aln[[1]][1])
control.extend.length=control.fragment.size-width(control.norm.aln[[1]][1])

# calculating coverage
for (i in 1:length(chipChrs)){
chip.cov<-c(chip.cov,coverage(chip.norm.aln[[i]],width=seqlengths(Mmusculus),extend=as.integer(chip.extend.length[i])))
#chip2.cov<-c(chip2.cov,coverage(chip.norm.aln[[i]],width=seqlengths(Mmusculus),extend=as.integer(chip2.extend.length[i])))
control.cov<-c(control.cov,coverage(control.norm.aln[[i]],width=seqlengths(Mmusculus),extend=as.integer(control.extend.length[i])))
print (paste((length(chipChrs)-i),"chromosomes left",sep=''))
}

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

# coverage to track data
tracks=function(x){return (as(x,"RangedData"))}
nozero=function(x){return (subset(x,score>0))}
# obtaining tracks
chip.cov.Track=lapply(chip.cov,tracks)
chip2.cov.Track=lapply(chip2.cov,tracks)
control.cov.Track=lapply(control.cov,tracks)
chip.cov.no.zero=lapply(chip.cov.Track,nozero)
chip2.cov.no.zero=lapply(chip2.cov.Track,nozero)
control.cov.no.zero=lapply(control.cov.Track,nozero)

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
 
# my way of writing bed files
# exporting the coverage no zero files as bedGraph files for each chromosome
#exportBed=function(cov.no.zero,bedFileName){
#export(cov.no.zero,bedFileName,"bedGraph")
#}
# have to run a loop instead of mapply so as to get distinct names for different chrmosome bedFiles
#for (i in 1:length(control.cov.no.zero)){
#exportBed(control.cov.no.zero[[i]],paste(controlLevels[i],".control",sep=''))
#exportBed(chip.cov.no.zero[[i]],paste(chipLevels[i],".chip",sep=''))
#}

# single line production of bed files is, distinct filename feature is sorted.
#mapply(exportBed,chip.cov.no.zero,MoreArgs=list(bedFileName="chip"))

# Saving chip and control coverage data as RData for session restore
save(chip.cov,file="chip.cov.RData")
save(control.cov,file="control.cov.RData")

# islands - regions of interest are contiguous segments of non-zero coverage - selecting with atleast with 5 (regions less then 5 overlapping reads are not likely to indicate real binding site)
chip.islands=mapply(slice,chip.cov,MoreArgs=list(lower=5))
control.islands=mapply(slice,control.cov,MoreArgs=list(lower=5))

# unlisting the top [[1]] list, use to [[1]] or the chromosome name as [['chr1']]
for (i in 1:length(chip.cov)){
chip.islands[[i]]<-chip.islands[[i]][[1]]
control.islands[[i]]<-control.islands[[i]][[1]]
}

# Views - container for storing a set of views on an arbitrary Vector object called subject. we need to compare chip islands for the same location os control islands on the same reference position.
#control.islands1=mapply(Views,control.cov,MoreArgs=list(start=start(chip.islands),end=end(chip.islands)))
for (i in 1:length(chip.islands)){
control.islands[[i]]<-Views(control.cov[[i]][[1]],start=start(chip.islands[[i]]),end=end(chip.islands[[i]]))
}

# finding the peak height (maximum number of overlapping reads)
chip.peaks<-lapply(chip.islands,viewMaxs)
control.peaks<-lapply(control.islands,viewMaxs)

# making frequency distribution table
chip.table<-lapply(chip.peaks,table)
control.table<-lapply(control.peaks,table)

# plotting chip and control distribution tables
for (i in 1:length(chip.table)){
jpeg(paste("chr",i,".jpg",sep=''))
plot(names(chip.table[[i]]),log10(chip.table[[i]]),main="Frequency Distributions",xlab="Peak Height",ylab='Log10 Frequency',type='l',col='red',sub=paste("chr",i))
lines(names(control.table[[i]]),log10(control.table[[i]]),col="blue")
legend('topright',legend=c('ChIP','Control'),col=c('red','blue'),lwd=1)
dev.off()
}
# plotting level of enrichment - peak heights
for (i in 1:length(chip.peaks)){
jpeg(paste("chr",i,".jpg",sep=''))
plot(chip.peaks[[i]],control.peaks[[i]],xlim=c(0,400),ylim=c(0,400),pch=20,col="black",main="ChIP vs Control",xlab="ChIP",ylab="Control",sub=paste("chr",i))
dev.off()
}

# x= chip.peaks, y=control.peaks
enrich=function(x,y){return(x>y)}
mapply(enrich,chip.peaks,MoreArgs=list(y=control.peaks))
##########################################
#Dummy code for reading in large bam files
chip=list()
# import bs genome package
library(BSgenome.Mmusculus.UCSC.mm9)
chrs=c(paste("chr",seq(1,19),sep=''),paste("chr",c('X','Y','M'),sep=''))
what=scanBamWhat()
number=seq(1,27,by=3)
number=number[-length(number)]
for (i in 1:(length(number)-1)){
which<-GRanges((chrs[number[i]:number[i+1]]),ranges=IRanges((number[i]:number[i+1]),c(seqlengths(Mmusculus)[number[i]:number[i+1]])))
if (i==2|i==3) {which=which[-1]}
params<-ScanBamParam(which=which,what=what)
print (which)
chip[[i]]=readAligned('/projects/globalscratch/sukhi/rotterdam/uniqueRot/s_1_1.fastq.chipseq.bowtie.SR.mm9.sam.srt.bam.unique',type="BAM",param=params)
}

##########################################
# looping data per chromosome
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

# length of each chromosome
chrLen=function(x){return(length(x@id))}
chipchrLength=unlist(lapply(chip,chrLen))
controlchrLength=unlist(lapply(control,chrLen))




##########################################

# End of Code

