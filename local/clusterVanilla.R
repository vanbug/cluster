##########################################
# looping data per chromosome
chip=list();control<-list()

# import bs genome package
library(BSgenome.Mmusculus.UCSC.mm9)
library(chipseq)
library(ShortRead)

# chip and control file user input
print("Enter the chip file")
chipFile=scan("stdin",what="character",n=1,quiet=TRUE)
print("Enter the control file")
controlFile=scan("stdin",what="character",n=1,quiet=TRUE)

chrs=c(paste("chr",seq(1,19),sep=''),paste("chr",c('X','Y','M'),sep=''))
what=scanBamWhat()
number=seq(1,22,by=1)

for(i in 1:length(number)){

which<-GRanges((chrs[number[i]]),ranges=IRanges(number[i],c(seqlengths(Mmusculus)[number[i]])))

# replacing NA with chr length
seqlengths(which)=seqlengths(Mmusculus)[number[i]]

params<-ScanBamParam(which=which,what=what)

# reading in files
chip=readAligned(chipFile,type="BAM",param=params)
control=readAligned(controlFile,type="BAM",param=params)
print(paste("Working with",i))

# removing NA's
control=control[which(is.na(control@strand)==F)]
chip=chip[which(is.na(chip@strand)==F)]

# filtering data
filt<-occurrenceFilter(withSread=FALSE)
control<-control[filt(control)]
chip<-chip[filt(chip)]

print ("Data Filtered")

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
chip.fragment.size=round(estimate.mean.fraglen(chip.norm.aln,method="SISSR"))
control.fragment.size=round(estimate.mean.fraglen(control.norm.aln,method="SISSR"))

print ("Working on Coverage")

# Investigating Coverage
# calculating extend length
chip.extend.length=chip.fragment.size-width(chip.norm.aln)
control.extend.length=control.fragment.size-width(control.norm.aln)

# calculating coverage
chip.cov<-coverage(chip.norm.aln,width=seqlengths(Mmusculus),extend=as.integer(chip.extend.length))
control.cov<-coverage(control.norm.aln,width=seqlengths(Mmusculus),extend=as.integer(control.extend.length))

print (paste("Writing cov for Chr",i))

# writing coverage externally
bedGraph=function(cov,chrs,type){
        pf<-runValue(cov[[1]])>0
        write.table(cbind(chrs,start(cov[[1]])[pf],
        format(end(cov[[1]]),scientific=F)[pf],
        runValue(cov[[1]])[pf]),file=paste(type,".bedGraph",sep=''),append=TRUE,sep="\t",
        quote=FALSE,row.names=FALSE,col.names=FALSE)
        }

bedGraph(chip.cov,chrs[i],type="chip")
bedGraph(control.cov,chrs[i],type="control")

print (paste("Done with",i))
}
