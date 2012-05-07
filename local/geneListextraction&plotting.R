# ctr9 gene list plotting
ctr9=read.csv('Galaxy30-[L317_nijCon(peaks__bed)].bed.anno',sep='\t')
geneName_freqTable=as.data.frame(table(ctr9$Gene.Name))
head(geneName_freqTable)
geneName_ordered=geneName_freqTable[with(geneName_freqTable,order(-Freq)),]
head(geneName_ordered)
history()
# plot command
ggplot(geneName_ordered,aes(x=factor(Var1[1:50]),y=factor(Freq[1:50]),fill=factor(Var1[1:50])))+geom_bar()+geom_text(aes(label=factor(Var1[1:50]),colour="Gene Name"))+opts(axis.title.x = theme_text(face="bold", colour="Black", size=10),axis.text.x = theme_text(angle=90,size=10))


ggplot(geneName_ordered,aes(x=factor(Var1[1:50]),y=factor(Freq[1:50]),fill=factor(Var1[1:50])))+geom_bar()+opts(title="ctr9 top 50 binding (bindings per gene) profile",axis.title.x = theme_text(face="bold", colour="Black", size=10),axis.text.x = theme_text(angle=90,size=10))+scale_x_discrete("Gene Names")+scale_y_discrete('Number of Bindings per Gene')


# promoter subset from whole annotation
promoter_subset=ctr9[grep("*promoter",ctr9$Annotation),]
promoter_geneList_table=data.frame(table(promoter_subset$Gene.Name))
promoter_geneList_table=promoter_geneList_table[with(promoter_geneList_table,order(-Freq)),]

# ctr9_promoterBindings.png
ggplot(promoter_geneList_table,aes(x=factor(Var1[1:65]),y=factor(Freq[1:65])))+geom_bar(fill="#515151")+coord_polar(start=1)+geom_text(aes(y=factor(Freq[1:65])+0.5,label=factor(Freq[1:65])))+scale_x_discrete('Maximum Peaks per Promoter Gene List',expand=c(10,0))+scale_y_discrete('Number of Peaks per Gene',expand=c(10,10))+opts(title="Ctr9 Promoter Binding Profile (top 65 Genes)",legend.position="top")

