#  Rscript for the bed2cvg
# Input chip and control file
chipFile=commandArgs()[5]
controlFile=commandArgs()[6]

# intializing libraries
library(BSgenome.Mmusculus.UCSC.mm9);library(ShortRead);library(chipseq);library(ff)
chip=read.table.ffdf(file=chipFile)
print ("Read Chip")
control=read.table.ffdf(file=controlFile)
print ("Read Control")

# setting parameters
fragment.size=200;chr.size<-seqlengths(Mmusculus)

#First one has to collect only the 22 chromosomes from your table 
    #1st option
    #chip.subset<-subset(chip,V1 %in% 
#c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chrX","chrY"))
#print ("Subsetting Chip")
    #control.subset<-subset(control,V1 %in% 
#c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chrX","chrY"))
#print ("Subsetting Control")
    #2nd option
 #    data.subset<-subset(data, V1!='chr1_random' & V1!='chr13_random' & V1!='chr4_random' & V1!='chr5_random' & V1!='chr7_random' & V1!='chr8_random' & V1!='chr9_random' & V1!='chrM' & V1!='chrUn_random' & V1!='chrX_random' & V1!='chrY_random')
    #checking if it worked
# SMALL HACK FOR WORKING WITH ffdf
chip.subset=chip
control.subset=control

#check<-subset(chip.subset,!duplicated(V1))
#    check

#GRanges
chip.GRanges<-GRanges(seqnames=as.character(chip.subset$V1[1:length(chip.subset$V1)]),IRanges(start=chip.subset$V2[1:length(chip.subset$V1)],
end=chip.subset$V3[1:length(chip.subset$V1)]),strand=chip.subset$V4[1:length(chip.subset$V1)])
print ("Chip Granged")
control.GRanges<-GRanges(seqnames=as.character(control.subset$V1[1:length(control.subset$V1)]),IRanges(start=control.subset$V2[1:length(control.subset$V1)],
end=control.subset$V3[1:length(control.subset$V1)]),strand=control.subset$V4[1:length(control.subset$V1)])
print ("Control Granged")
seqlengths(chip.GRanges)<-chr.size[c(1,10:19,2:9,20:21)]
seqlengths(control.GRanges)<-chr.size[c(1,10:19,2:9,20:21)]
chip.GRanges <- resize(chip.GRanges,width=fragment.size)
control.GRanges <- resize(control.GRanges,width=fragment.size)
print ("resized")

#save(data.GRanges,file='data.GRanges')

#Removing duplicates for the same start positions and same strand orientation
# cleaning has no big effect, just the duplicates explained above are removed
chip.clean<-chip.GRanges[!duplicated(paste(start(chip.GRanges),strand(chip.GRanges),sep=":"))]
control.clean<-control.GRanges[!duplicated(paste(start(control.GRanges),strand(control.GRanges),sep=":"))]

print ("Data Cleaned")
# NORMALIZING LIBRARY SIZE FROM CHIP AND CONTROL using sampling
#sizeChip<-length(chip.clean)
#sizeControl<-length(control.clean)
#smallest.size<-min(sizeControl,sizeChip)

#print ("Sampling Datasets")

# normalizing chip and control library size
#chip.norm.aln<-sample(chip.clean,smallest.size,replace=FALSE)
#control.norm.aln<-sample(control.clean,smallest.size,replace=FALSE)

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
chip.cov<-coverage(chip.clean,width=seqlengths(Mmusculus))
control.cov<-coverage(control.clean,width=seqlengths(Mmusculus))

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
	runValue(cov)[pf]),file=paste(type,".yahoo",sep=''),append=TRUE,sep="\t",
        quote=FALSE,row.names=FALSE,col.names=FALSE)
        }
############## CONTROL IS OFF - SWITCH IT ON TO PRODUCE CONTROL BED FILE ##############################
for (chrs in c(paste('chr',seq(1,19),sep=''),'chrX','chrY')){
bedGraph(chip.cov[[chrs]],chrs,type=chipFile)
bedGraph(control.cov[[chrs]],chrs,type=controlFile)
print (paste("Processing",chrs))
}
