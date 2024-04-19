
##############SNPs##############
snps=read.table('/Users/awilder/Documents/RhinoWGS/GATK/AllWR_CerSimSim1.0_raw_snps.tab.gz',header=T)

variable='QD'
plot(density(snps[,variable],na.rm=T),xlab=variable,main='',lwd=2)
abline(v=2,lwd=2)
abline(v=5,lwd=2)
text(x=1.25,y=0.0425,labels='2')
text(x=4.25,y=0.0425,labels='5')
text(x=0,y=0.04,labels=length(which(snps[,variable]<2)))
text(x=8.25,y=0.04,labels=length(which(snps[,variable]<5)))

variable='FS'
plot(density(log10(snps[,variable]),na.rm=T),xlab=variable,main='',lwd=2,xaxt='n')
axis(1, at=(-1:6)*0.5, labels=c(-0.32,1,10^((1:6)*0.5)))
abline(v=log10(32),lwd=2)
text(x=log10(20),y=0.7,labels='32')
text(x=log10(20),y=0.65,labels=length(which(snps[,variable]>32)))
abline(v=log10(60),lwd=2)
text(x=log10(80),y=0.7,labels='60')
text(x=log10(100),y=0.65,labels=length(which(snps[,variable]>100)))

variable='SOR'
plot(density(snps[,variable],na.rm=T,bw=0.2),xlab=variable,main='',lwd=2)
abline(v=3,lwd=2)
text(x=2.75,y=1,labels='3')
text(x=2,y=0.9,labels=length(which(snps[,variable]>3)))
abline(v=4,lwd=2)
text(x=4.75,y=1,labels='4')
text(x=5.2,y=0.9,labels=length(which(snps[,variable]>4)))

variable='MQ'
plot(density(snps[,variable],na.rm=T),xlab=variable,main='',lwd=2,xlim=c(0,55))
abline(v=40,lwd=2)
text(x=45,y=3,labels='40')
text(x=48,y=2.75,labels=length(which(snps[,variable]<40)))
abline(v=35,lwd=2)
text(x=30,y=3,labels='28')
text(x=28,y=2.75,labels=length(which(snps[,variable]<28)))

variable='MQRankSum'
plot(density(snps[,variable],na.rm=T,bw=0.2),xlab=variable,main='',lwd=2)#,xlim=c(-5,5))
# abline(v=-2,lwd=2)
# text(x=-1,y=25,labels='-0.5')
# text(x=-1.25,y=23,labels=length(which(snps[,variable]<(-0.5))))
abline(v=-3,lwd=2)
text(x=-5,y=1.2,labels='-5')
text(x=-6,y=1.1,labels=length(which(snps[,variable]<(-3))))
abline(v=-2,lwd=2)
text(x=1,y=1.2,labels='-2')
text(x=2,y=1.1,labels=length(which(snps[,variable]<(-2))))
# abline(v=2,lwd=2)
# text(x=2.5,y=25,labels='2')
# text(x=2.75,y=23,labels=length(which(snps[,variable]>(2))))

variable='ReadPosRankSum'
plot(density(snps[,variable],na.rm=T,bw=0.1),xlab=variable,main='',lwd=2,xlim=c(-5,5))
abline(v=-2,lwd=2)
text(x=-2.5,y=0.7,labels='-2')
text(x=-2.75,y=0.6,labels=length(which(snps[,variable]<(-2))))
abline(v=-1,lwd=2)
text(x=-1.4,y=0.7,labels='-1')
text(x=-1.5,y=0.6,labels=length(which(snps[,variable]<(-1))))


library(VennDiagram)
QDfilt=5
FSfilt=32
SORfilt=3
MQfilt=28
MQRSfilt=-2
RPRSfilt=-2

QD=which(snps$QD<QDfilt)
FS=which(snps$FS>FSfilt)
SOR=which(snps$SOR>SORfilt)
MQ=which(snps$MQ<MQfilt)#using 40 is a lot of SNPs!
MQRankSum=which(snps$MQRankSum<MQRSfilt)
ReadPosRankSum=which(snps$ReadPosRankSum <RPRSfilt)

barplot(sapply(list(QD, FS,SOR,MQ,MQRankSum, ReadPosRankSum),length),names.arg=apply(cbind(c('QD<','FS>','SOR>','MQ<','MQRS>','RPRS>'),as.character(c(QDfilt,FSfilt,SORfilt,MQfilt,MQRSfilt,RPRSfilt))),1,function(x) paste0(x[1],x[2])),ylab='No. of SNPs failing filter',main='Conservative')

venn.diagram(x = list(QD, MQ,MQRankSum),category.names = apply(cbind(c('QD<','MQ<','MQRS>'),as.character(c(QDfilt,MQfilt,MQRSfilt))),1,function(x) paste0(x[1],x[2])),filename = '/Users/awilder/Documents/RhinoWGS/GATK/snp_venn_diagramm.png', output = TRUE , imagetype="png" ,lwd = 2, lty = 'blank', fill = c('yellow', 'purple', 'green'),main='Conservative')


############Indels####################
# Starter recommendations:
# QD < 2.0indels
# FS > 200.0
# ReadPosRankSum < -20.0

indels=read.table('AllWR_CerSimSim1.0_raw_indels.tab.gz',header=T)

variable='QD'
plot(density(indels[,variable],na.rm=T),xlab=variable,main='',lwd=2)
abline(v=2,lwd=2)
abline(v=5,lwd=2)
text(x=0,y=0.035,labels='2')
text(x=6.25,y=0.035,labels='5')
text(x=-1,y=0.03,labels=length(which(indels[,variable]<2)))
text(x=8.25,y=0.03,labels=length(which(indels[,variable]<5)))

variable='FS'
plot(density(log10(indels[,variable]),na.rm=T),xlab=variable,main='',lwd=2,xaxt='n')
axis(1, at=(-1:6)*0.5, labels=c(-0.32,1,10^((1:6)*0.5)))
abline(v=log10(32),lwd=2)
text(x=log10(20),y=0.6,labels='32')
text(x=log10(20),y=0.55,labels=length(which(indels[,variable]>32)))
abline(v=log10(60),lwd=2)
text(x=log10(80),y=0.6,labels='60')
text(x=log10(100),y=0.55,labels=length(which(indels[,variable]>100)))

variable='SOR'
plot(density(indels[,variable],na.rm=T,bw=0.2),xlab=variable,main='',lwd=2)
abline(v=3,lwd=2)
text(x=2.75,y=0.85,labels='3')
text(x=2,y=0.8,labels=length(which(indels[,variable]>3)))
abline(v=4,lwd=2)
text(x=4.75,y=0.85,labels='4')
text(x=5.2,y=0.8,labels=length(which(indels[,variable]>4)))

variable='MQ'
plot(density(indels[,variable],na.rm=T),xlab=variable,main='',lwd=2,xlim=c(0,55))
abline(v=40,lwd=2)
text(x=45,y=0.45,labels='40')
text(x=48,y=0.4,labels=length(which(indels[,variable]<40)))
abline(v=28,lwd=2)
text(x=30,y=0.45,labels='28')
text(x=32,y=0.4,labels=length(which(indels[,variable]<28)))

variable='MQRankSum'
plot(density(indels[,variable],na.rm=T,bw=0.2),xlab=variable,main='',lwd=2)#,xlim=c(-5,5))
# abline(v=-2,lwd=2)
# text(x=-1,y=25,labels='-0.5')
# text(x=-1.25,y=23,labels=length(which(indels[,variable]<(-0.5))))
abline(v=-3,lwd=2)
text(x=-5,y=1,labels='-5')
text(x=-6,y=0.9,labels=length(which(indels[,variable]<(-3))))
abline(v=-2,lwd=2)
text(x=1,y=1,labels='-2')
text(x=2,y=0.9,labels=length(which(indels[,variable]<(-2))))
# abline(v=2,lwd=2)
# text(x=2.5,y=25,labels='2')
# text(x=2.75,y=23,labels=length(which(indels[,variable]>(2))))

variable='ReadPosRankSum'
plot(density(indels[,variable],na.rm=T,bw=0.1),xlab=variable,main='',lwd=2,xlim=c(-5,5))
abline(v=-2,lwd=2)
text(x=-2.5,y=0.55,labels='-2')
text(x=-2.75,y=0.5,labels=length(which(indels[,variable]<(-2))))
abline(v=-1,lwd=2)
text(x=-1.4,y=0.55,labels='-1')
text(x=-1.5,y=0.5,labels=length(which(indels[,variable]<(-1))))

indels=merge(indels,wrbc,by.x=c("CHROM","POS"),by.y=c('Chr','Pos'),all.x=T)

keepinds=which(indels$QD>5 & indels $FS<32 & indels $SOR<3 & (indels$ReadPosRankSum>-2 | is.na(indels$ReadPosRankSum)))
exclinds=which(1:nrow(indels) %in% keepinds==FALSE)

plot(density(indels[keepinds,'NWR_AltF'],na.rm=T,bw=0.05),col='red',ylim=c(0,6),xlab='AltFreq',main='Indels',lwd=2)
lines(density(indels[exclinds,'NWR_AltF'],na.rm=T,bw=0.05),col='red',lty=2,lwd=2)
lines(density(indels[keepinds,'SWR_AltF'],na.rm=T,bw=0.05),col='dodgerblue',lwd=2)
lines(density(indels[exclinds,'SWR_AltF'],na.rm=T,bw=0.05),col='dodgerblue',lty=2,lwd=2)
legend('topright',legend=c('SWR','NWR','SWR filtered','NWR filtered'),lty=c(1,1,2,2),col=rep(c('dodgerblue','red'),2),lwd=2)

