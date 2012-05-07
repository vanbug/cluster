# Rscript to plot pie charts using ggplot2 library for the inputted 2-cloumn region annotation files
# Author : Sukhdeep Singh
# Organization : Biotec
########################################################################################################

library('ggplot2');

# reading files in
#fileName=commandArgs[5]
files=grep('*regions',list.files(),value=T)
dev.new(width=18, height=10)
for (i in 1:length(files)){
sub=read.csv(files[i],sep='\t',header=FALSE)
x_sub=strsplit(files[i],'\\.')[[1]][1]
ggplot(sub,aes(x=factor(V1),y=factor(V2),fill=factor(V1)))+geom_bar(width=0.8,alpha=0.8)+coord_polar(theta='x')+geom_text(aes(y=(factor(V2)+0.5),label=factor(V2),colour='No. of Gene Targets'))+scale_x_discrete(paste(x_sub,"binding profile pie chart"))+scale_y_discrete("Number of targets")
ggsave(paste(files[i],".png",sep=""))
}
