####R scripts#####
#checking the final lines of a chromosome in a dataframe
tail(data[data[,1]=='chr1',])

#subsetting columns
h3ash2l<-h3ash2l[,c(1,2,3,5,6)]

#Grepping collors for plots
colors()[grep("red",colors())]

##Plot example
plot(c(-2.5e5,2.5e5),c(1,y.exp.max), type= "n", ylab="Reads/Million",xaxt='n',xlab="Chromosomal position relative to the viewpoint (100Kb)",
main=paste("Interaction regions close to viewpoint (zoom in 500 KB)"),ylim=c(1,y.exp.max))
points(exp.500k$p_start,exp.500k$expRPMs,col='red',pch=20)
lines(exp.500k$p_start,exp.500k$expRPMs,col='blue')
abline(v=0,lty=3,col="grey",lwd=2)
axis(1, at = c(-2.5e5,-2e5,-1.5e5,-1e5,-0.5e5,0,0.5e5,1e5,1.5e5,2e5,2.5e5),lab=c(-2.5,-2.0,-1.5,-1.0,-0.5,0,0.5,1,1.5,2,2.5))
legend("topright", legend = c(expLabeled), fill=c("blue"), cex=0.65 )


##Counting reads per million
expRPMs  <- round(your n reads / (your librarysize/1e6))

####Processing Bed files to obtain the coverage and upload in UCSC genome browser####
#Main libraries
library(BSgenome.Mmusculus.UCSC.mm9)
library(ShortRead)

#Reading tables from Bed files in R
chip<-read.table("your file",skip=1)

#Preparing GRanges files for one chromosome (example: chromosome 1)
chr="chr1"
fragment.size=200
chr.size<-seqlengths(Mmusculus)[chr]
chip.GRanges<-GRanges(seqnames=as.vector(chr1.data$V1),IRanges(start=chr1.data$V2,end=chr1.data$V3),strand=chr1.data$V6)
seqlengths(chip.GRanges)<-chr.size[1:21]
chip.GRanges <- resize(chip.GRanges, width = fragment.size)
save(chip.GRange,file="filename")

##Coverage
chip.cov<-coverage(chip.GRanges)
plot(start(chip.cov[[1]]),runValue(chip.cov[[1]]),type='l')
pf <-runValue(chip.cov[[chr]]) >0
write.table(cbind(chr, start(chip.cov[[chr]])[pf], format((end(chip.cov[[chr]])+1), scientific=F)[pf],runValue(chip.cov[[chr]])[pf]), file="Your file name",append=TRUE, sep="\t", quote=FALSE,row.names=FALSE, col.names=FALSE)
===================================================================================================
##Pipeline for Bed files##

library(BSgenome.Mmusculus.UCSC.mm9)
library(ShortRead)
data<-read.table("file",skip=1)
fragment.size=200
chr.size<-seqlengths(Mmusculus)

#First one has to collect only the 21 chromossomes from your table 
    #1st option
    data.subset<-subset(data,V1 %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chrX","chrY"))
    #2nd option
    data.subset<-subset(data, V1!='chr1_random' & V1!='chr13_random' & V1!='chr4_random' & V1!='chr5_random' & V1!='chr7_random' & V1!='chr8_random' & V1!='chr9_random' & V1!='chrM' & V1!='chrUn_random' & V1!='chrX_random' & V1!='chrY_random')
    #checking if it worked - gives the unique elements - good one
    check<-subset(data.subset,!duplicated(V1))
    check
#GRanges
data.GRanges<-GRanges(seqnames=as.character(data$V1),IRanges(start=data$V2,end=data$V3),strand=data$V6)
seqlengths(data.GRanges)<-chr.size[c(1,10:19,2:9,20:21)]
data.GRanges <- resize(data.GRanges, width = fragment.size)
save(data.GRanges,file='data.GRanges')
#Removing duplicates
#chr.GRanges<-chr.GRanges[!duplicated(paste(start(chr.GRanges),strand(chr.GRanges),sep=":")),]
dataclean.GRanges<-data.GRanges[!duplicated(paste(start(data.GRanges),strand(data.GRanges),sep=":"))]
#Coverage
data.cov<-coverage(data.GRanges)
pdf('data.pdf')
plot(start(data.cov[[1]]),runValue(data.cov[[1]]),type='l')
dev.off()

##writing the table to upload in UCSC genome browser
for (chr in c(paste('chr',seq(1,19),sep=''),'chrX','chrY')){
	pf <-runValue(data.cov[[chr]]) >0
	write.table(cbind(chr, start(data.cov[[chr]])[pf], format((end(data.cov[[chr]])+1), scientific=F)[pf],runValue(data.cov[[chr]])[pf]), file="dataCov.bedGraph",append=TRUE, sep="\t", quote=FALSE,row.names=FALSE, col.names=FALSE)
}
================
#altering the start and end positions when they are equal by adding 1 at the end position
data<-read.table("dataCov.bedGraph")
bo<- data[,3] == data[,2]
data[,3][bo] <- data[,3][bo] +1
write.table(data,file="dataCov.bedGraph",sep="\t",quote=FALSE,row.names=FALSE, col.names=FALSE,format(data,scientific=F))

#inserting a header in the bedGraph
vi dataCov.bedGraph
track type=bedGraph name="dataCov.bedGraph" description="data_coverage"
:wq ##save and quit

 ##Making the BigWig file from the bedGraph##
fetchChromSizes mm9 > mouseChr	###add 1bp at the end of each chromossome end position using vi editor
bedGraphToBigWig dataCov.bedGraph mouseChr+1 data.bw

#collecting all the runvalues in a vector file 
pf <- NULL
for (i in 1:length(data.cov)){
tmp <- runValue(data.cov[[i]]) >0
pf <- c(pf, tmp)
}

#Adding tracks in UCSC genome browser from a http server using BigWig files
track type=bigWig name="ash2l" description="ash2l_coverage" bigDataUrl=http://panda:panda@shire.bccs.uib.no/ucsc_bam/BigWig/H3K4_subunits_francis/ash2l.bw
track type=bigWig name="cxxc1" description="cxxc1_coverage" bigDataUrl=http://panda:panda@shire.bccs.uib.no/ucsc_bam/BigWig/H3K4_subunits_francis/cxxc1.bw
track type=bigWig name="wdr82" description="wdr82_coverage" bigDataUrl=http://panda:panda@shire.bccs.uib.no/ucsc_bam/BigWig/H3K4_subunits_francis/wdr82.bw
track type=bigWig name="ptip" description="ptip_coverage" bigDataUrl=http://panda:panda@shire.bccs.uib.no/ucsc_bam/BigWig/H3K4_subunits_francis/ptip.bw
track type=bigWig name="menin" description="menin_coverage" bigDataUrl=http://panda:panda@shire.bccs.uib.no/ucsc_bam/BigWig/H3K4_subunits_francis/menin.bw
track type=bigWig name="h3k4" description="h3k4_coverage" bigDataUrl=http://panda:panda@shire.bccs.uib.no/ucsc_bam/BigWig/H3K4_subunits_francis/h3k4.bw
track type=bigWig name="h3k27" description="h3k27_coverage" bigDataUrl=http://panda:panda@shire.bccs.uib.no/ucsc_bam/BigWig/H3K4_subunits_francis/h3k27.bw
===================================================================================================================================================
###Calculating CpG islands content
cpg<-read.table('cpg_mouse.txt') ##CpG positions were retrieved from UCSC genome browser Tables
cpg<-cpg[,2:4]
cpg<-subset(cpg,V2 %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chrX","chrY"))
cpg.GRanges<-GRanges(seqnames=as.character(cpg$V2),IRanges(start=cpg$V3,end=cpg$V4))
seqlengths(cpg.GRanges)<-chr.size[1:21]

ash2l.cpg<-as.matrix(findOverlaps(ash2l.5peak.GRanges,cpg.GRanges))
ash2l.cpg.percent<-(length(ash2l.cpg[,1])/length(ash2l.5peak.GRanges))*100

cxxc.cpg<-as.matrix(findOverlaps(cxxc.5peak.GRanges,cpg.GRanges))
cxxc.cpg.percent<-(length(cxxc.cpg[,1])/length(cxxc.5peak.GRanges))*100

wdr82.cpg<-as.matrix(findOverlaps(wdr82.5peak.GRanges,cpg.GRanges))
wdr82.cpg.percent<-(length(wdr82.cpg[,1])/length(wdr82.5peak.GRanges))*100

ptip.cpg<-as.matrix(findOverlaps(ptip.5peak.GRanges,cpg.GRanges))
ptip.cpg.percent<-(length(ptip.cpg[,1])/length(ptip.5peak.GRanges))*100

k4.cpg<-as.matrix(findOverlaps(k4.5peak.GRanges,cpg.GRanges))
k4.cpg.percent<-(length(k4.cpg[,1])/length(k4.5peak.GRanges))*100

k27.cpg<-as.matrix(findOverlaps(k27.5peak.GRanges,cpg.GRanges))
k27.cpg.percent<-(length(k27.cpg[,1])/length(k27.5peak.GRanges))*100

> wdr82.cpg.percent
[1] 32.73092
> cxxc.cpg.percent
[1] 61.14554
> ash2l.cpg.percent
[1] 31.71073
> ptip.cpg.percent
[1] 3.373058
> k4.cpg.percent
[1] 46.74151
> k27.cpg.percent
[1] 36.99185

#Finding bivalency regions
bimat<-as.matrix(findOverlaps(h3k4.5peak.GRanges,h3k27.5peak.GRanges))
bih3k4.GRanges<-h3k4.5peak.GRanges[bimat[,1],]
bih3k27.GRanges<-h3k27.5peak.GRanges[bimat[,2],]

    ##Converting GRanges in data.frame
        bih3k27data<-as.data.frame(bih3k27.GRanges,row.names = NULL, optional = FALSE)
        bih3k4data<-as.data.frame(bih3k4.GRanges,row.names = NULL, optional = FALSE)
    ##Collecting the lines based on the the NoReads
        bik4reads<-bik4data[bih3k4data$NoReads>=bih3k27data$NoReads,]
        bik27reads<-bik27data[bih3k27data$NoReads>bih3k4data$NoReads,]
        bireads<-rbind(bik4reads,bik27reads)
        bireads.GRanges<-GRanges(seqnames=as.character(bireads$seqnames),IRanges(start=bireads$PeakPosition-200,end=bireads$PeakPosition+200),NoReads=bireads$NoReads,pvalue=bireads$pvalue,PeakPosition=bireads$PeakPosition)
        seqlengths(bireads.GRanges)<-chr.size[1:21]
    ##Bivalent positions of the subunits
        ash2l.bireads<-as.matrix(findOverlaps(ash2l.5peak.GRanges,bireads.GRanges))
        ash2l.bireads.percent<-(length(ash2l.bireads[,1])/length(ash2l.5peak.GRanges))*100

        cxxc.bireads<-as.matrix(findOverlaps(cxxc.5peak.GRanges,bireads.GRanges))
        cxxc.bireads.percent<-(length(cxxc.bireads[,1])/length(cxxc.5peak.GRanges))*100

        wdr82.bireads<-as.matrix(findOverlaps(wdr82.5peak.GRanges,bireads.GRanges))
        wdr82.bireads.percent<-(length(wdr82.bireads[,1])/length(wdr82.5peak.GRanges))*100

        ptip.bireads<-as.matrix(findOverlaps(ptip.5peak.GRanges,bireads.GRanges))
        ptip.bireads.percent<-(length(ptip.bireads[,1])/length(ptip.5peak.GRanges))*100

        > ash2l.bireads.percent
        [1] 3.427530
        > cxxc.bireads.percent
        [1] 9.278678
        > wdr82.bireads.percent
        [1] 1.204819
        > ptip.bireads.percent
        [1] 0.4681847

    
#subsetting NoReads in GRanges
###don't use loops for big files in R###
tmp<-NULL
for (i in 1:length(bik4)){
    if (values(bik4)[i,1]>=values(bik27)[i,1]) {
        tmp<-bik4[i,]
        bi<-append(bi,tmp)
} else {
        tmp<-bik27[i,]
        bi<-append(bi,tmp)

}

###Processing the Bed files significant peaks after using CCAT 3.0 peak caller 
datapeak<-read.table('filepath',col.names=c('chromosome','PeakPosition','start','end','NoReads','pvalue'))
datapeak.GRanges<-GRanges(seqnames=as.character(datapeak$chromosome),IRanges(start=datapeak$start,end=datapeak$end),NoReads=datapeak$NoReads,pvalue=datapeak$'pvalue',PeakPosition=datapeak$PeakPosition)
seqlengths(datapeak.GRanges)<-chr.size[1:21]

data.5peak<-datapeak[datapeak[,6]<=0.05,]
data.5peak.GRanges<-GRanges(seqnames=as.character(data.5peak$chromosome),IRanges(start=data.5peak$PeakPosition-200,end=data.5peak$PeakPosition+200),NoReads=data.5peak$NoReads,pvalue=data.5peak$'pvalue',PeakPosition=data.5peak$PeakPosition)
seqlengths(data.5peak.GRanges)<-chr.size[1:21]

##Getting ensembl atributes and TSS position
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
ensembl.genes <- getBM(attributes = c('chromosome_name','start_position','end_position','ensembl_gene_id','external_gene_id','strand'), mart = ensembl)
head(ensembl.genes)

##change_chr_name
ensembl.genes$chromosome_name=paste('chr',ensembl.genes$chromosome_name,sep='')

## add_tss
ensembl.genes$tss <- ifelse(ensembl.genes$strand==1,ensembl.genes$start_position,ensembl.genes$end_position)
head(ensembl.genes)

##GRanges of tss
ensembl.genes.tss<-subset(ensembl.genes,chromosome_name %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chrX","chrY"))
tss.GRanges<-GRanges(seqnames=as.character(ensembl.genes.tss$chromosome_name),IRanges(start=ensembl.genes.tss$tss-5000,end=ensembl.genes.tss$tss+1000),strand=ensembl.genes.tss$strand,ensembl_gene_id=ensembl.genes.tss$ensembl_gene_id,gene_name=ensembl.genes.tss$external_gene_id,TSS=ensembl.genes.tss$tss)
seqlengths(tss.GRanges)<-chr.size[1:21]
length(tss.GRanges)
[1] 36404

##Getting active promoter regions overlapping tss against h3k4
actpromat<-as.matrix(findOverlaps(tss.GRanges,h3k4.5peak.GRanges))
dim(actpromat)
tss.actpro.GRanges<-tss.GRanges[actpromat[,1],]
length(tss.actpro.GRanges)
[1] 18113

####Getting active promoter regions overlapping tss against h3k4 and then against h3k27
bimat<-as.matrix(findOverlaps(h3k4.5peak.GRanges,h3k27.5peak.GRanges))
h3k4bi.GRanges<-h3k4.5peak.GRanges[-bimat[,1]]

actpro2mat<-as.matrix(findOverlaps(tss.GRanges,h3k4bi.GRanges,minoverlap=201))
dim(actpro2mat)
tss.actpro2.GRanges<-tss.GRanges[actpro2mat[,1],]
length(tss.actpro2.GRanges)
[1] 13261

##Getting active promoters (h3k4) that seats in CpG islands
actprocpgmat<-as.matrix(findOverlaps(tss.actpro.GRanges,cpg.GRanges))
dim(actprocpgmat)
tss.actprocpg.GRanges<-tss.actpro.GRanges[actprocpgmat[,1],]
length(tss.actprocpg.GRanges)
[1] 16651

##Getting active promoters (h3k4 and not h3k27) that seats in CpG islands
actpro2cpgmat<-as.matrix(findOverlaps(tss.actpro2.GRanges,cpg.GRanges))
dim(actpro2cpgmat)
tss.actpro2cpg.GRanges<-tss.actpro2.GRanges[actpro2cpgmat[,1],]
length(tss.actpro2cpg.GRanges)
[1] 11984


##Calculating how many peaks are not in promoter regions
ash2l.promat<-as.matrix(findOverlaps(ash2l.5peak.GRanges,tss.GRanges))
dim(ash2l.promat)
ash2l.nonpro<-length(ash2l.5peak.GRanges)-length(ash2l.promat[,1])
(ash2l.nonpro/length(ash2l.5peak.GRanges))*100


##correlation between CXXC1 and Wdr82
load('cxxxG')
load('wdr82')
cxwdmat<-as.matrix(findOverlaps(cxxc.5peak.GRanges,wdr82.5peak.GRanges))
cxwd.GRanges<-cxxc.5peak.GRanges[cxwdmat[,1],]
wdcx.GRanges<-wdr82.5peak.GRanges[cxwdmat[,2],]
cxwd<-((values(cxwd.GRanges)[,1])*1000000)/length(cxxc.GRanges)
wdcx<-((values(wdcx.GRanges)[,1])*1000000)/length(wdr82.GRanges)
pdf('cxwdCorrNorm.pdf')
plot(cxwd,wdcx)
dev.off()

##Retrieving gene names for the significant peaks
data.annotation.frame<-data.frame()
for (chr in paste('chr',c(seq(1,19),'X','Y'),sep='')){
	
	peaks.chr<-subset(data,chromosome==chr)
	peaks.chr.irange<-IRanges(
	start=peaks.chr$PeakPosition,
	end=peaks.chr$PeakPosition)
	ensembl.genes.chr<-subset(ensembl.genes,chromosome_name==chr)
	ensembl.genes.chr.irange<-IRanges(
				start=ensembl.genes.chr$tss,
				end=ensembl.genes.chr$tss)

	gene.index<-nearest(peaks.chr.irange,ensembl.genes.chr.irange)

	target.gene.name <-ensembl.genes.chr$external_gene_id[gene.index]
	target.ensembl.id <-ensembl.genes.chr$ensembl_gene_id[gene.index]
	target.ensembl.tss <-ensembl.genes.chr$tss[gene.index]
	distance.to.gene <-target.ensembl.tss-peaks.chr$PeakPosition
	chr.frame <-cbind(peaks.chr,target.ensembl.id,
			target.gene.name,distance.to.gene)
	data.annotation.frame <-rbind(data.annotation.frame,
			chr.frame)
	gc()
}
head(data.annotation.frame)


###Plotting the coverage in TSS
ensembl.genes.tss

load('ash2lG')
ash2l.cov<-coverage(ash2l.GRanges)

###Plotting one chromosome
ash2l.views <- Views(ash2l.cov[['chr1']], start=ensembl.genes.tss[ensembl.genes.tss[,1]=='chr1',]$tss-1000, end=ensembl.genes.tss[ensembl.genes.tss[,1]=='chr1',]$tss+1000)
a <- as.matrix (viewApply(ash2l.views, as.vector))
b<-t(a)
pdf('ashchr1view.pdf')
plot(1, type="n", axes=T, xlim=c(1,2001),ylim=c(500,20000),xlab="TSS", ylab="SumReads")
lines(colSums(b), type="l", col="blue")
dev.off()

###Plotting all chromosomes in all promoters
ash2l.views<-NULL
for (chr in paste('chr',c(seq(1,19),'X','Y'),sep='')){

    temp <- Views(ash2l.cov[[chr]], start=ensembl.genes.tss[ensembl.genes.tss[,1]==chr,]$tss-2000, end=ensembl.genes.tss[ensembl.genes.tss[,1]==chr,]$tss+2000)
    ash2l.views<-c(ash2l.views,temp)
    
    }

ash2l.views.matrix<-NULL
for (i in 1:length(ash2l.views))  {
    temp<-t(as.matrix(viewApply(ash2l.views[[i]],as.vector)))
    ash2l.views.matrix<-rbind(ash2l.views.matrix,temp)
}


###Plotting using number of reads per million in order to compare among samples
ash2lRPMs  <- round(colSums(ash2l.views.matrix) / (length(ash2l.GRanges)/1e6))
pdf('ash2lTSSviewRPM.pdf')
plot(1, type="n", xaxt='n', xlim=c(1,2001),ylim=c(min(ash2lRPMs)-500,max(ash2lRPMs)),xlab="TSS position", ylab="Sum of RPMs",main="Ash2l binding to promoter regions")
lines(ash2lRPMs, type="l", col="blue")
axis(1,pos=min(ash2lRPMs)-500,at=c(1,500,1001,1500,2001),lab=c(-1000,-500,0,500,1000))
dev.off()

###Plotting Views using just colSums
pdf('ash2lTSSview.pdf')
plot(1, type="n", xaxt='n', xlim=c(1,2001),ylim=c(500,400000),xlab="TSS position", ylab="Sum of the number of reads",main="Ash2l binding to promoter regions")
lines(colSums (ash2l.views.matrix), type="l", col="blue")
axis(1,pos=0,at=c(1,500,1001,1500,2001),lab=c(-1000,-500,0,500,1000))
dev.off()



    
