# overlappings in R
suppressPackageStartupMessages(library('ChIPpeakAnno'))

file1=commandArgs(TRUE)[1]
file2=commandArgs(TRUE)[2]
file3=commandArgs(TRUE)[3]
file4=commandArgs(TRUE)[4]

# getting tf names
tf1=commandArgs(TRUE)[5]; tf2=commandArgs(TRUE)[6]; tf3=commandArgs(TRUE)[7]; tf4=commandArgs(TRUE)[8]

# reading peak files
peaks1=read.csv(file1,sep='\t',header=FALSE)
peaks2=read.csv(file2,sep='\t',header=FALSE)
peaks3=read.csv(file3,sep='\t',header=FALSE)
peaks4=read.csv(file4,sep='\t',header=FALSE)

# filtering important info
peaks1=peaks1[,1:3]; peaks2=peaks2[,1:3]; peaks3=peaks3[,1:3]; peaks4=peaks4[,1:3]

# converting them to ranges
peakRanges1=BED2RangedData(peaks1)
peakRanges2=BED2RangedData(peaks2)
peakRanges3=BED2RangedData(peaks3)
peakRanges4=BED2RangedData(peaks4)

# finding overlaps
#overlaps=findOverlappingPeaks(peakRange1,peakRange2,maxgap=100,multiple=c(TRUE,TRUE),NameOfPeaks1='bod1L',NameOfPeaks2='cxxc1')
    
        
# save as high quality tiff
tiff(file=paste(tf1,"-",tf2,"-",tf3,"_",tf4,"_","overlaps.tiff",sep=""),width=17.15,height=17.15,units="cm",res=1200, pointsize=10, compression = "lzw")
makeVennDiagram(RangedDataList(peakRanges1,peakRanges2,peakRanges3,peakRanges4),NameOfPeaks=c(tf1,tf2,tf3,tf4),totalTest=max(length(peaks1$V1),length(peaks2$V1),length(peaks3$V1),length(peaks4$V1)),cex=1,maxgap=0,counts.col="purple")
dev.off()


# pie chart of overlap features
#pie(table(overlaps$OverlappingPeaks$overlapFeature))
