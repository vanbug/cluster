# Rscript to normalize chip and control and find coverage from the bam files 
# Author : Sukhdeep Singh
# BUGS : Warnings in the read aligned params, shorten the what parameters to required ones, set TRUE for bamSimpleCigar(params) & bamReverseComplement(params)
# USAGE: R --vanilla --slave --args "chip.bam" "control.bam" "chr7" "outputDir" < bam2cvg_singleChr.R
##########################################

# Fetching chip and control file along with the chromosome name to be processed
chipFile=commandArgs()[5]
controlFile=commandArgs()[6]
chr=commandArgs()[7]

# setting output directory
outputDir=commandArgs()[8]
setwd(outputDir)

# importing libraries
#install.packages('/home/genom/sukhdeeps/BSgenome.Mmusculus.UCSC.mm9_1.3.16.tar.gz','/home/genom/sukhdeeps/R/x86_64-redhat-linux-gnu-library/2.13/')
library(BSgenome.Mmusculus.UCSC.mm9);library(chipseq);library(ShortRead)

# readAligned parameters
what=scanBamWhat()
which<-GRanges(chr,ranges=IRanges(1,c(seqlengths(Mmusculus)[chr])))

# important to fill in the NA with correct chr length
seqlengths(which)=seqlengths(Mmusculus)[chr]

# reading in files with assigned parameters
params<-ScanBamParam(which=which,what=what)
chip=readAligned(chipFile,type="BAM",param=params)
control=readAligned(controlFile,type="BAM",param=params)

# removing NA's
control=control[which(is.na(control@strand)==F)]
chip=chip[which(is.na(chip@strand)==F)]

# filtering data - no alignQuality FILTER
#filt<-occurrenceFilter(withSread=FALSE) # occurrence Filter reduces the amount of reads significantly as it retains the min=1 and max=1 (min and max no. of times a read occurs after filter)
#control<-control[filt(control)]
#chip<-chip[filt(chip)]

#print ("Data Filtered")

# variable declarations
chip.norm.aln<-list(); control.norm.aln<-list();
chip.fragment.size<-c(); control.fragment.size<-c();
chip.cov<-c(); control.cov<-c();

# NORMALIZING LIBRARY SIZE FROM CHIP AND CONTROL - in GFP & MEN, control has more reads
sizeControl<-length(control@id)
sizeChip<-length(chip@id)
smallest.size<-min(sizeControl,sizeChip)

print ("Starting Normalization")

# normalizing chip and control library size
chip.norm.aln<-sample(chip,smallest.size,replace=FALSE)
control.norm.aln<-sample(control,smallest.size,replace=FALSE)

# fragment length extension
chip.fragment.size=round(estimate.mean.fraglen(chip.norm.aln,method="SISSR")) # method=coverage (#120) & method=correlation (#270)
control.fragment.size=round(estimate.mean.fraglen(control.norm.aln,method="SISSR"))

print (chip.fragment.size)
print ("Working on Coverage")

# Investigating Coverage
# calculating extend length
chip.extend.length=chip.fragment.size-width(chip.norm.aln)
control.extend.length=control.fragment.size-width(control.norm.aln)

# calculating coverage
chip.cov<-coverage(chip.norm.aln,width=seqlengths(Mmusculus),extend=as.integer(chip.extend.length)) # extend extends the read coverage display from length of 36 to 200 in reference to the human genome
control.cov<-coverage(control.norm.aln,width=seqlengths(Mmusculus),extend=as.integer(control.extend.length))


# writing coverage externally
bedGraph=function(cov,chrs,type){
	pf<-runValue(cov[[1]])>0
	write.table(cbind(chrs,format(start(cov[[1]]),scientific=F)[pf],
	format(end(cov[[1]]),scientific=F)[pf],
	runValue(cov[[1]])[pf]),file=paste(type,".bedGraph",sep=''),append=TRUE,sep="\t",
	quote=FALSE,row.names=FALSE,col.names=FALSE)
	}

bedGraph(chip.cov,chr,type=chipFile)
bedGraph(control.cov,chr,type=controlFile)
print ("Done")


## END
