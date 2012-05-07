# R script sorting the bed to bigWig conversion
# bedFile input
bedFile=commandArgs()[5]

# checking if dir is made and grepping for track line
ff=file.path(getwd(),"bigwigBattlefield")
if (file_test("-d",ff)==FALSE) {system('mkdir bigwigBattlefield');system(paste("grep track",bedFile,"> bigwigBattlefield/grepTrack"))} else 
{system(paste("grep track",bedFile,"> bigwigBattlefield/grepTrack"))}

# sorting bed file and removing track line
if (file.info('bigwigBattlefield/grepTrack')$size>0) {system(paste("sed -e 1d",bedFile,">bigwigBattlefield/customBedA"));system(paste("sort -k1,1 -k2,2n bigwigBattlefield/customBedA >bigwigBattlefield/sortedBed"))} else 
{system(paste("sort -k1,1 -k2,2n",bedFile,">bigwigBattlefield/sortedBed"))}

# bedtobigWig conversion
system(paste("./bedGraphToBigWig bigwigBattlefield/sortedBed mm9.chrom bigwigBattlefield/",bedFile,".bw",sep=""))

# removing cached files
system('rm bigwigBattlefield/customBedA bigwigBattlefield/sortedBed bigwigBattlefieldgrepTrack')
print ("Done")	
