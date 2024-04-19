
snps=read.table('/Users/awilder/Documents/RhinoWGS/GATK/DicBic_AdptTrimPaired_bt2_CerSim1.0_raw_snps.tab.gz',header=T)

variable='QD'
plot(density(snps[,variable],na.rm=T),xlab=variable,main='',lwd=2)
abline(v=10,lwd=2)
abline(v=21,lwd=2)
text(x=8,y=0.05,labels='10')
text(x=18.25,y=0.05,labels='21')
text(x=8,y=0.04,labels=length(which(snps[,variable]<10)))
text(x=18.25,y=0.04,labels=length(which(snps[,variable]<21)))

variable='FS'
plot(density(log10(snps[,variable]),na.rm=T),xlab=variable,main='',lwd=2,xaxt='n')
axis(1, at=(-1:6)*0.5, labels=c(-0.32,1,10^((1:6)*0.5)))
abline(v=log10(32),lwd=2)
text(x=log10(20),y=0.07,labels='32')
text(x=log10(20),y=0.065,labels=length(which(snps[,variable]>32)))
abline(v=log10(60),lwd=2)
text(x=log10(80),y=0.07,labels='60')
text(x=log10(100),y=0.065,labels=length(which(snps[,variable]>100)))

variable='SOR'
plot(density(snps[,variable],na.rm=T,bw=0.2),xlab=variable,main='',lwd=2)
abline(v=3,lwd=2)
text(x=2.75,y=1,labels='3')
text(x=2,y=0.9,labels=length(which(snps[,variable]>3)))
abline(v=4,lwd=2)
text(x=4.75,y=1,labels='4')
text(x=5.2,y=0.9,labels=length(which(snps[,variable]>4)))

variable='MQ'
plot(density(snps[,variable],na.rm=T),xlab=variable,main='',lwd=2,xlim=c(20,45))
abline(v=32,lwd=2)
text(x=35,y=0.25,labels='32')
text(x=35,y=0.2,labels=length(which(snps[,variable]<32)))
abline(v=35,lwd=2)
text(x=30,y=3,labels='28')
text(x=28,y=2.75,labels=length(which(snps[,variable]<28)))

variable='MQRankSum'
plot(density(snps[,variable],na.rm=T,bw=0.2),xlab=variable,main='',lwd=2)#,xlim=c(-5,5))
# abline(v=-2,lwd=2)
# text(x=-1,y=25,labels='-0.5')
# text(x=-1.25,y=23,labels=length(which(snps[,variable]<(-0.5))))
abline(v=-7,lwd=2)
text(x=-6,y=0.2,labels='-5')
text(x=-7,y=0.15,labels=length(which(snps[,variable]<(-5))))
abline(v=-2,lwd=2)
text(x=1,y=0.2,labels='-2')
text(x=2,y=0.15,labels=length(which(snps[,variable]<(-2))))
# abline(v=2,lwd=2)
# text(x=2.5,y=25,labels='2')
# text(x=2.75,y=23,labels=length(which(snps[,variable]>(2))))

variable='ReadPosRankSum'
plot(density(snps[,variable],na.rm=T,bw=0.1),xlab=variable,main='',lwd=2,xlim=c(-10,10))
abline(v=-2,lwd=2)
text(x=-2.5,y=0.5,labels='-2')
text(x=-2.75,y=0.4,labels=length(which(snps[,variable]<(-2))))
abline(v=-1,lwd=2)
text(x=1.4,y=0.5,labels='-1')
text(x=1.5,y=0.4,labels=length(which(snps[,variable]<(-1))))

library(VennDiagram)
QDfilt=21
FSfilt=32
SORfilt=3
MQfilt=32
MQRSfilt=-2
RPRSfilt=-2

QD=which(snps$QD<QDfilt)
FS=which(snps$FS>FSfilt)
SOR=which(snps$SOR>SORfilt)
MQ=which(snps$MQ<MQfilt)#using 40 is a lot of SNPs!
MQRankSum=which(snps$MQRankSum<MQRSfilt)
ReadPosRankSum=which(snps$ReadPosRankSum <RPRSfilt)

barplot(sapply(list(QD, FS,SOR,MQ,MQRankSum, ReadPosRankSum),length),names.arg=apply(cbind(c('QD<','FS>','SOR>','MQ<','MQRS>','RPRS>'),as.character(c(QDfilt,FSfilt,SORfilt,MQfilt,MQRSfilt,RPRSfilt))),1,function(x) paste0(x[1],x[2])),ylab='No. of SNPs failing filter',main='Conservative')

venn.diagram(x = list(QD, MQ,MQRankSum),category.names = apply(cbind(c('QD<','MQ<','MQRS>'),as.character(c(QDfilt,MQfilt,MQRSfilt))),1,function(x) paste0(x[1],x[2])),filename = '#14_venn_diagramm.png', output = TRUE , imagetype="png" ,lwd = 2, lty = 'blank', fill = c('yellow', 'purple', 'green'),main='Conservative')#, height = 480 , width = 480 , resolution = 300, compression = "lzw",  cex = 1, fontface = "bold", fontfamily = "sans", cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "outer", cat.pos = c(-27, 27, 135), cat.dist = c(0.055, 0.055, 0.085), cat.fontfamily = "sans", rotation = 1 )

indels=read.table('/Users/awilder/Documents/RhinoWGS/GATK/DicBic_AdptTrimPaired_bt2_CerSim1.0_raw_indels.tab.gz',header=T)

variable='QD'
plot(density(indels[,variable],na.rm=T),xlab=variable,main='',lwd=2)
abline(v=10,lwd=2)
abline(v=21,lwd=2)
text(x=8,y=0.05,labels='10')
text(x=18.25,y=0.05,labels='21')
text(x=8,y=0.04,labels=length(which(indels[,variable]<10)))
text(x=18.25,y=0.04,labels=length(which(indels[,variable]<21)))

variable='FS'
plot(density(log10(indels[,variable]),na.rm=T),xlab=variable,main='',lwd=2,xaxt='n')
axis(1, at=(-1:6)*0.5, labels=c(-0.32,1,10^((1:6)*0.5)))
abline(v=log10(32),lwd=2)
text(x=log10(20),y=0.07,labels='32')
text(x=log10(20),y=0.065,labels=length(which(indels[,variable]>32)))
abline(v=log10(60),lwd=2)
text(x=log10(80),y=0.07,labels='60')
text(x=log10(100),y=0.065,labels=length(which(indels[,variable]>100)))

variable='SOR'
plot(density(indels[,variable],na.rm=T,bw=0.2),xlab=variable,main='',lwd=2)
abline(v=3,lwd=2)
text(x=2.75,y=1,labels='3')
text(x=2,y=0.9,labels=length(which(indels[,variable]>3)))
abline(v=4,lwd=2)
text(x=4.75,y=1,labels='4')
text(x=5.2,y=0.9,labels=length(which(indels[,variable]>4)))

variable='MQ'
plot(density(indels[,variable],na.rm=T),xlab=variable,main='',lwd=2,xlim=c(20,45))
abline(v=32,lwd=2)
text(x=35,y=0.12,labels='32')
text(x=35,y=0.08,labels=length(which(snps[,variable]<32)))
abline(v=35,lwd=2)
text(x=30,y=.12,labels='28')
text(x=28,y=0.08,labels=length(which(snps[,variable]<28)))

variable='MQRankSum'
plot(density(indels[,variable],na.rm=T,bw=0.2),xlab=variable,main='',lwd=2)#,xlim=c(-5,5))
# abline(v=-2,lwd=2)
# text(x=-1,y=25,labels='-0.5')
# text(x=-1.25,y=23,labels=length(which(indels[,variable]<(-0.5))))
abline(v=-7,lwd=2)
text(x=-6,y=0.2,labels='-5')
text(x=-7,y=0.15,labels=length(which(indels[,variable]<(-5))))
abline(v=-2,lwd=2)
text(x=1,y=0.2,labels='-2')
text(x=2,y=0.15,labels=length(which(indels[,variable]<(-2))))
# abline(v=2,lwd=2)
# text(x=2.5,y=25,labels='2')
# text(x=2.75,y=23,labels=length(which(indels[,variable]>(2))))

variable='ReadPosRankSum'
plot(density(indels[,variable],na.rm=T,bw=0.1),xlab=variable,main='',lwd=2,xlim=c(-10,10))
abline(v=-2,lwd=2)
text(x=-2.5,y=0.5,labels='-2')
text(x=-2.75,y=0.4,labels=length(which(indels[,variable]<(-2))))
abline(v=-1,lwd=2)
text(x=1.4,y=0.5,labels='-1')
text(x=1.5,y=0.4,labels=length(which(indels[,variable]<(-1))))

QDfilt=22
FSfilt=32
SORfilt=3
MQRSfilt=-5
RPRSfilt=-2

QD=which(indels$QD<QDfilt)
FS=which(indels $FS>FSfilt)
SOR=which(indels $SOR>SORfilt)
MQRankSum=which(indels $MQRankSum<MQRSfilt)
ReadPosRankSum=which(indels $ReadPosRankSum <RPRSfilt)

barplot(sapply(list(QD, FS,SOR,MQRankSum, ReadPosRankSum),length),names.arg=apply(cbind(c('QD<','FS>','SOR>','MQRS>','RPRS>'),as.character(c(QDfilt,FSfilt,SORfilt,MQRSfilt,RPRSfilt))),1,function(x) paste0(x[1],x[2])),ylab='No. of SNPs failing filter',main='Conservative')

