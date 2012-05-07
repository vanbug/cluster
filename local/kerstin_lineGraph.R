# kerstin's graphs
a=read.csv('~/Desktop/Weight_curve_Kerstin.csv',sep='\t')
b=data.frame(t(a))
meltedA=melt(a,id=("Mouse"))
meltedB=data.frame(meltedA[,c(1,2)],data.frame(gsub('g$','',meltedA$value)))
colnames(meltedB)[3]="weight"
colnames(meltedB)[2]="days"

# a is the csv.file
a=read.csv('~/Desktop/Weight_curve_Kerstin2.csv',sep='\t',dec=',')
#ggplot(meltedB,aes(days,weight,group=Mouse))+geom_line(aes(color=Mouse),size=2)
meltedA=melt(a,id=c("Mouse","Type"))
d=as.Date(gsub('^X','',meltedA$variable),format="%d.%m.%Y")
meltedB=data.frame(meltedA[,c(1,2)],d,data.frame(gsub('g$','',meltedA$value)))
colnames(meltedB)[3]="days"
colnames(meltedB)[4]="weight"
ggplot(na.omit(meltedB),aes(factor(days),weight,group=Mouse))+geom_line(aes(color=Mouse),size=2)+facet_grid(Type~.)

# Kerstin1
ggplot(na.omit(meltedB),aes(factor(days),weight,group=Mouse))+geom_line(aes(color=Mouse),size=2)+facet_grid(Type~Mouse)+scale_y_discrete(breaks=seq(1,40,by=1),name="Weight (in gms)")+scale_color_discrete(name="Mouse")+scale_x_discrete("Date")+opts(title="Control and Knockout Mouse Weight Comparisons")
ggsave('~/Desktop/Kerstin1.png')

# Kerstin2
ggplot(na.omit(meltedB),aes(factor(days),weight,group=Mouse))+geom_line(aes(color=Mouse),size=2)+facet_grid(Type~.)+scale_y_discrete(breaks=seq(15,40,by=1),name="Weight (in gms)")+scale_color_discrete(name="Mouse")+scale_x_discrete("Date")+opts(title="Control and Knockout Mouse Weight Comparisons")
ggsave('~/Desktop/Kerstin2.png')
# for control and knockout graphs
# meltedK=data.frame(meltedK[,c(1,2,3)],data.frame(gsub('g$','',meltedK$value)))
# colnames(meltedK)[4]="weight"
# colnames(meltedK)[3]="days"
# qplot example
qplot(date, value,data=graph1,
      geom="line",
      colour=variable,xlab="",
      ylab="",
      size= I(1))+
   scale_y_continuous(limits = c(-0.3,0.3))+
   scale_colour_manual("Variable",c(Line1="red",Line2="blue"))+ 
   opts(legend.size="none",
        aspect.ratio = 2/(1+sqrt(5)))


# helmut graphs
ggplot(data=melted,aes(x=Regions,fill=factor(variable),y=value))+geom_histogram(binwidth=0.5,position=position_identity())+opts(title="N vs C terminal binding profile comparisons")+scale_y_continuous('No. of significant peaks')
