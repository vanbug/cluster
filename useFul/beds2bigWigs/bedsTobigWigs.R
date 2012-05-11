# R script sorting and converting the bed to bigWig conversion
# requirements folder should have only the bed files to be converted, so make a new folder and put your bed files in

# fetching a list of bed files
allBeds=list.files()
beds=grep(".bedGraph",allBeds,value=TRUE)

system('mkdir bigwigBattlefield')
for (i in 1:length(beds)){
system(paste("grep track",beds[i],"> bigwigBattlefield/grepTrack"))

# sorting bed file and removing track line
if (file.info('bigwigBattlefield/grepTrack')$size>0) {system(paste("sed -e 1d",beds[i],">bigwigBattlefield/customBedA"));system(paste("sort -k1,1 -k2,2n bigwigBattlefield/customBedA >bigwigBattlefield/sortedBed"));system('rm bigwigBattlefield/customBedA')} else 
{system(paste("sort -k1,1 -k2,2n",beds[i],">bigwigBattlefield/sortedBed"))}

# bedClip for the removing big end coordinates
system(paste("~/src/useFul/ucsc/utilities/bedClip bigwigBattlefield/sortedBed ~/src/useFul/ucsc/genomeIndex/mm9.chrom bigwigBattlefield/clipped.bed"))
# bedtobigWig conversion
system(paste("~/src/useFul/ucsc/utilities/bedGraphToBigWig bigwigBattlefield/clipped.bed ~/src/useFul/ucsc/genomeIndex/mm9.chrom bigwigBattlefield/",beds[i],".bw",sep=""))

# removing cached files
system('rm bigwigBattlefield/sortedBed bigwigBattlefield/grepTrack bigwigBattlefield/clipped.bed bigwigBattlefield/bedsTobigWigs.Rout.bw')
print (paste("Done with",beds[i]))
print ("BigWigs are located in bigwigBattlefield folder in the current dir")
}
############################################################################################
