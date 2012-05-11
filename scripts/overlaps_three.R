# overlappings in R
suppressPackageStartupMessages(library('ChIPpeakAnno'))

file1=commandArgs(TRUE)[1]
file2=commandArgs(TRUE)[2]
file3=commandArgs(TRUE)[3]

# getting tf names
tf1=commandArgs(TRUE)[4]; tf2=commandArgs(TRUE)[5]; tf3=commandArgs(TRUE)[6]

# reading peak files
peaks1=read.csv(file1,sep='\t',header=FALSE)
peaks2=read.csv(file2,sep='\t',header=FALSE)
peaks3=read.csv(file3,sep='\t',header=FALSE)

# filtering important info
peaks1=peaks1[,1:3]; peaks2=peaks2[,1:3]; peaks3=peaks3[,1:3];

# converting them to ranges
peakRanges1=BED2RangedData(peaks1)
peakRanges2=BED2RangedData(peaks2)
peakRanges3=BED2RangedData(peaks3)

# finding overlaps
#overlaps=findOverlappingPeaks(peakRange1,peakRange2,maxgap=100,multiple=c(TRUE,TRUE),NameOfPeaks1='bod1L',NameOfPeaks2='cxxc1')
# save graphs as pdf
#pdf(paste(tf1,"-",tf2,"-",tf3,"_","overlaps.pdf",sep=""))
# save as high quality tiff
tiff(file=paste(tf1,"-",tf2,"-",tf3,"_","overlaps.tiff",sep=""),width=17.15,height=17.15,units="cm",res=1200, pointsize=10, compression = "lzw")
makeVennDiagram(RangedDataList(peakRanges1,peakRanges2,peakRanges3),NameOfPeaks=c(tf1,tf2,tf3),totalTest=max(length(peaks1$V1),length(peaks2$V1),length(peaks3$V1)),maxgap=0,cex=1,counts.col="purple")
dev.off()

# pie chart of overlap features
#pie(table(overlaps$OverlappingPeaks$overlapFeature))
