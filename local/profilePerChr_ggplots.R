# read file
ctr9=read.csv('Galaxy30-[L317_nijCon(peaks__bed)].bed.anno',sep='\t')
library('gtools','ggplot2','RColorBrewer')

# chr list collapsed
chr=lapply(mixedsort(unique(ctr9$Chr)),function(x){return(ctr9[which(ctr9$Chr==x),])})

# x is the read annotation file
collapser=function(x){
x$Annotation=gsub("*intron.*","intron",x$Annotation)
x$Annotation=gsub("*promoter.*","promoter",x$Annotation)
x$Annotation=gsub("*exon.*","exon",x$Annotation)
x$Annotation=gsub("*5'.*","5'",x$Annotation)
x$Annotation=gsub("*3'.*","3'",x$Annotation)
x$Annotation=gsub("*TTS.*","TTS",x$Annotation)
return(data.frame(table(x$Annotation)))
}

chr_collapsed=lapply(chr,collapser)

# melting the list for a specific id
chr_melted=melt(chr_collapsed,id='Var1')

# adding chr names to the melted list
plot_list=cbind(chr_melted,unlist(mapply(function(x,y){return(rep(x,y))},mixedsort(unique(ctr9$Chr)),data.frame(table(chr_melted$L1))$Freq)))

# removing unwanted columns
plot_list=plot_list[,-c(2,4)]

# renaming colnames
colnames(plot_list)[3]='chrNames'
colnames(plot_list)[1]='Annotations'
colnames(plot_list)[2]='Bindings'

#reordering columns of the plot_list data.frame
plot_list=data.frame(plot_list$chrNames,plot_list$Annotations,plot_list$Bindings)

# renaming columns
colnames(plot_list)[3]="Bindings"
colnames(plot_list)[2]="Annotations"
colnames(plot_list)[1]="chrNames"

# bindings per chromsome list
chrBinding=data.frame(mixedsort(unique(ctr9$Chr)),sapply(chr,function(x)length(x[[1]])))
colnames(chrBinding)[1]='Chr'
colnames(chrBinding)[2]='Bindings'

# chrBinding plot
ggplot(chrBinding,aes(x=mixedsort(Chr),y=Bindings))+geom_bar()

# plot using qplot
qplot(factor(chrNames),y=factor(Bindings),data=plot_list,fill=Annotations)+geom_bar()+scale_fill_manual(values=c(brewer.pal(8,"Dark2")))+coord_polar(theta='y')

# facets graph
qplot(factor(chrNames),y=factor(Bindings),data=plot_list,fill=Annotations)+facet_grid(facets=.~Annotations)

# histogram + facets
qplot(x=factor(chrNames),y=factor(Bindings),z=factor(Annotations),data=plot_list,fill=Annotations)+facet_grid(facets=.~Annotations)+geom_bar(binwidth=12)

# amazing pattern graph
qplot(x=factor(chrNames),z=factor(Annotations),data=plot_list)+facet_grid(facets=.~Annotations)+coord_polar(theta='y')

# nice bar chart 
ggplot(plot_list,aes(factor(Bindings)))+geom_histogram(binwidth=0.2)+facet_grid(chrNames~Annotations)

# beta phase bar chart (good)
ggplot(plot_list,aes(x=chrNames,y=Bindings,fill=Annotations))+geom_bar()

# color palatte applied
ggplot(plot_list,aes(x=chrNames,y=Bindings,fill=Annotations))+geom_bar(alpha=0.8)+scale_fill_manual(values=(brewer.pal(8,"Dark2")))

