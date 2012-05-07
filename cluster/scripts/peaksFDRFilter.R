# R script to filter the peak xls file with peaks given FDR (default 10%)
# Author : Sukhdeep Singh
# usage : Rscript peaksFDRFilter.R peaksByMacs.xls 5
#####################################################################################################

# reading arguments
file=commandArgs(TRUE)[1]
inputFDR=commandArgs(TRUE)[2]

# setting up working directory according to the file names

# reading file with header skipped
data=read.csv(file,sep='\t',skip=23)

# renaming columns
colnames(data)[7]="-10*log10(pvalue)"
colnames(data)[9]="FDR(%)"

if (is.na(inputFDR)==TRUE){inputFDR=10} else {inputFDR=inputFDR}

# filtering applied
passed=which(data[,9]<=as.numeric(inputFDR))
filtered=data[passed,]
homer=cbind(filtered[,c(1,2,3)],paste("filteredPeak",seq(1:length(filtered[,1])),sep=""),filtered[,7])
colnames(homer)[4]="PeakNumber"
colnames(homer)[5]="-10*log10(pvalue)"

# outputting message
print(paste(length(passed),"peaks passed FDR filter of",inputFDR,"out of",length(data[,9])))
print(paste("Total fold enrichment of whole sample profile peaks w.r.t to used control:",sum(filtered[8])))

# writing xls file
write.table(filtered,paste(dirname(file),'/fdrFilterPeaks_',inputFDR,'.bed',sep=''),sep='\t',quote=FALSE,row.names=FALSE)
write.table(homer,paste(dirname(file),'/homerPeaks_',inputFDR,'.bed',sep=''),sep='\t',quote=FALSE,row.names=FALSE)

