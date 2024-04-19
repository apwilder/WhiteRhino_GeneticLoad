library(RColorBrewer)
library(reshape)
library(ggplot2)
library(stringr)
library(venneuler)

#READ IN DATA TABLE
gtsd=read.table('~/USS/conservation-genetics/nwr/Tables/AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_GERP_Ann_2BR_DerGT.txt',header=T,colClasses=c('character','integer','character','character','numeric','character',rep('integer',22),rep('character',3),'numeric'),sep='\t')

gtsd=gtsd[complete.cases(gtsd[,7:28]),]
gtsd=gtsd[substr(gtsd$Chr, 1, 1)=='J',]
gtsd=gtsd[which(!is.na(gtsd$DicBicMic) & !is.na(gtsd$DicBicMin)),]

cols=c(brewer.pal(9,'Oranges')[c(6)],'dodgerblue')

gtsd$NWRderF=rowSums(gtsd[,7:15])/(rowSums(!is.na(gtsd[,7:15]))*2)
gtsd$SWRderF=rowSums(gtsd[,16:28])/(rowSums(!is.na(gtsd[,16:28]))*2)
gtsd$derF=rowSums(gtsd[,7:28])/(rowSums(!is.na(gtsd[,7:28]))*2)

gtsd$RS[which(gtsd$RS<0)]=0
gtsd$phyloP[which(gtsd$phyloP<0)]=0

#mammalia-odb busco gene list
busco=read.table('/home/centos/USS/conservation-genetics/nwr/Tables/buscogenes.txt')
busco$busco='Y'
busco$V1=as.character(busco$V1)
gtsd=merge(gtsd,busco,by.x='GeneName',by.y='V1',all.x=T)
gtsd$busco[is.na(gtsd$busco)]='N'
gtsd=gtsd[,c(2:31,1,32:38)]

#boxplot function
scatterbox=function(df,title,ylab=NULL){
  boxplot(df[,1]~df[,2],data=df,main=title,ylab=ylab,xlab=NULL)
  stripchart(df[1:9,1]~as.factor(df[1:9,2]),data=df,col=brewer.pal(9,'Oranges')[c(6)],add=T, vertical = TRUE, method = "jitter",pch=16)
  stripchart(c(rep(NA,length(1:9)),df[10:18,1])~as.factor(df[1:18,2]),data= df,col='dodgerblue',add=T, vertical = TRUE, method = "jitter",pch=1)
  stripchart(c(rep(NA,length(1:18)),df[19:22,1])~as.factor(df[1:22,2]),data= df,col='dodgerblue',add=T, vertical = TRUE, method = "jitter",pch=16)
  siglev=c(1,0.05,0.01,0.001,0.0001)
  names(siglev)=c("","*","**","***","****")
  lmh=aov(df[,1]~df[,2])
  thsd=TukeyHSD(lmh)
  p=thsd$df[4]
  mtext(names(siglev)[max(which(siglev>p))],at=1.5,cex=2,line=-2)
  print(p)
}

calledalleles=colSums(!is.na(gtsd[,7:28]))*2

###Heterozygosity#######
het=colSums(gtsd[,7:28]==1,na.rm=T)/(2426503069*(nrow(gtsd)/9549126))#9549126=raw SNPs
het=as.data.frame(het)
het$pop=factor(c(rep('NWR',9),rep('SWR',13)),levels=c('NWR','SWR'))
scatterbox(het,'Heterozygous sites',ylab='Heterozygosity')

#####private alleles#####
NWR_derT=colSums(gtsd[which(gtsd$SWRderF==0 & gtsd$NWRderF>0),7:15],na.rm=T)
SWR_derT=colSums(gtsd[which(gtsd$NWRderF==0 & gtsd$SWRderF>0),16:28],na.rm=T)
derT=c(NWR_derT, SWR_derT)
derT=as.data.frame(derT)
derT$pop=factor(c(rep('NWR',9),rep('SWR',13)),levels=c('NWR','SWR'))
scatterbox(derT,'Number of private deleterious alleles (RS>4)')

shared=sum(gtsd$SWRderF>0 & gtsd$NWRderF>0)
NWR_priv=sum(gtsd$SWRderF==0 & gtsd$NWRderF>0)
SWR_priv=sum(gtsd$SWRderF>0 & gtsd$NWRderF==0)
v <- venneuler(c(NWR=NWR_priv, SWR=SWR_priv, "NWR&SWR"=shared))
plot(v,col=cols)

#####deleterious alleles#####
#private deleterious
NWR_derT=colSums(gtsd[which(gtsd$SWRderF==0 & gtsd$NWRderF>0 & gtsd$RS>4),7:15],na.rm=T)
SWR_derT=colSums(gtsd[which(gtsd$NWRderF==0 & gtsd$SWRderF>0 & gtsd$RS>4),16:28],na.rm=T)
derT=c(NWR_derT, SWR_derT)
derT=as.data.frame(derT)
derT$pop=factor(c(rep('NWR',9),rep('SWR',13)),levels=c('NWR','SWR'))
scatterbox(derT,'Number of private deleterious alleles (RS>4)')

#Deleterious hom/het boxplots
derT1=as.data.frame(colSums(gtsd[which(gtsd$RS>4),7:28]==2,na.rm=T))
derT1=cbind(derT1,colSums(gtsd[which(gtsd$RS>4),7:28]==1,na.rm=T))
derT1$pop=factor(c(rep('NWR',9),rep('SWR',13)),levels=c('NWR','SWR'))
names(derT1)[1:2]=c('hom','het')
derT1$type='Conserved (GERP>4)'

derT2=as.data.frame(colSums(gtsd[which(gtsd$RS>4),7:28]==2,na.rm=T)/colSums(gtsd[!is.na(gtsd$RS),7:28]==2,na.rm=T))
derT2=cbind(derT2,colSums(gtsd[which(gtsd$RS>4),7:28]==1,na.rm=T)/colSums(gtsd[!is.na(gtsd$RS),7:28]==1,na.rm=T))
derT2$pop=factor(c(rep('NWR',9),rep('SWR',13)),levels=c('NWR','SWR'))
names(derT2)[1:2]=c('hom','het')
derT2$type='Proportion conserved (GERP>4)'

derT3=as.data.frame(colSums(gtsd[which(gtsd$Impact=='MODERATE' & gtsd$RS>4 & gtsd$busco=='Y'),7:28]==2,na.rm=T))
derT3=cbind(derT3,colSums(gtsd[which(gtsd$Impact=='MODERATE' & gtsd$RS>4 & gtsd$busco=='Y'),7:28]==1,na.rm=T))
derT3$pop=factor(c(rep('NWR',9),rep('SWR',13)),levels=c('NWR','SWR'))
names(derT3)[1:2]=c('hom','het')
derT3$type='Missense conserved'
TukeyHSD(aov(derT3$hom~derT3$pop))
TukeyHSD(aov(derT3$het~derT3$pop))

derT4=as.data.frame(colSums(gtsd[which(gtsd$Impact=='MODERATE' & gtsd$RS>4 & gtsd$busco=='Y'),7:28]==2,na.rm=T)/colSums(gtsd[which(!is.na(gtsd$Impact) & !is.na(gtsd$RS) & gtsd$busco=='Y'),7:28]==2,na.rm=T))
derT4=cbind(derT4,colSums(gtsd[which(gtsd$Impact=='MODERATE' & gtsd$RS>4 & gtsd$busco=='Y'),7:28]==1,na.rm=T)/colSums(gtsd[which(!is.na(gtsd$Impact) & !is.na(gtsd$RS) & gtsd$busco=='Y'),7:28]==2,na.rm=T))
derT4$pop=factor(c(rep('NWR',9),rep('SWR',13)),levels=c('NWR','SWR'))
names(derT4)[1:2]=c('hom','het')
derT4$type='Proportion missense conserved'
TukeyHSD(aov(derT4$hom~derT4$pop))
TukeyHSD(aov(derT4$het~derT4$pop))

derT5=as.data.frame(colSums(gtsd[which(gtsd$Impact=='HIGH' & gtsd$busco=='Y'),7:28]==2,na.rm=T))
derT5=cbind(derT5,colSums(gtsd[which(gtsd$Impact=='HIGH' & gtsd$busco=='Y'),7:28]==1,na.rm=T))
derT5$pop=factor(c(rep('NWR',9),rep('SWR',13)),levels=c('NWR','SWR'))
names(derT5)[1:2]=c('hom','het')
derT5$type='LOF'
TukeyHSD(aov(derT5$hom~derT5$pop))
TukeyHSD(aov(derT5$het~derT5$pop))

derT6=as.data.frame(colSums(gtsd[which(gtsd$Impact=='HIGH' & gtsd$busco=='Y'),7:28]==2,na.rm=T)/colSums(gtsd[which(!is.na(gtsd$Impact) & gtsd$busco=='Y'),7:28]==2,na.rm=T))
derT6=cbind(derT6,colSums(gtsd[which(gtsd$Impact=='HIGH' & gtsd$busco=='Y'),7:28]==1,na.rm=T)/colSums(gtsd[which(!is.na(gtsd$Impact) & !is.na(gtsd$RS) & gtsd$busco=='Y'),7:28]==2,na.rm=T))
derT6$pop=factor(c(rep('NWR',9),rep('SWR',13)),levels=c('NWR','SWR'))
names(derT6)[1:2]=c('hom','het')
derT6$type='Proportion LOF'
TukeyHSD(aov(derT6$hom~derT6$pop))
TukeyHSD(aov(derT6$het~derT6$pop))


derT=melt(rbind(derT1,derT3,derT5,derT2,derT4,derT6),id.vars=c('pop','type'))
derT$type=factor(derT$type,levels=c("Conserved (GERP>4)","Missense conserved",'LOF',"Proportion conserved (GERP>4)","Proportion missense conserved",'Proportion LOF'))

ggplot(data=derT,aes(x=variable,y=value,col=pop,group = interaction(variable, pop))) +facet_wrap(.~type,scale='free',nrow=2)+geom_boxplot(color='black')+geom_point(position = position_jitterdodge())+scale_color_manual(values=cols) +theme(text = element_text(size=10),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(),axis.ticks = element_blank(),legend.title=element_blank(),legend.key=element_blank(),legend.key.height = unit(0.5, "cm"),legend.key.width = unit(0.4, "cm"))+scale_x_discrete('', labels=c('Homozygous','Heterozygous'))+ylab('')

#######genetic load enrichment in ROH######
library(ggplot2)
library(reshape)
library(stringr)
library(RColorBrewer)

nwror=brewer.pal(9,'Oranges')[c(6)]
swrbl='dodgerblue'

rohn=read.table('NWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_roh_NWRfilt_RG_VarCounts.txt',header=T)
rohs=read.table('SWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_roh_SWRfilt_RG_VarCounts.txt',header=T)
rohs$Sample=as.factor(rohs$Sample)
rohn$Sample=as.factor(rohn$Sample)

rohdels=function(df,intervals){
  rohdel=matrix(nrow=length(intervals)-1,ncol=nlevels(df$Sample))
  for (sample in levels(df$Sample)){
    dfs=df[which(df$Sample==sample),]
    for (i in 1:(length(intervals)-1)){
      lower=intervals[i]*1e6
      upper=intervals[i+1]*1e6
      #rohdel[i,which(levels(dfs$Sample)==sample)]=mean(dfs$MODERATERS4[which(dfs$Length>lower & dfs$Length<=upper & dfs$LOWRS>0)]/dfs$LOWRS[which(dfs$Length>lower & dfs$Length<=upper & dfs$LOWRS>0)],na.rm=T)
      rohdel[i,which(levels(dfs$Sample)==sample)]=mean(dfs$RS4[which(dfs$Length>lower & dfs$Length<=upper & dfs$RS>10)]/dfs$RS[which(dfs$Length>lower & dfs$Length<=upper & dfs$RS>10)],na.rm=T)
    }
  }
  return(rohdel)
}

intervals=c(0.1,0.5,1,3,8)
rohdel=rohdels(rohs,intervals=intervals)
rohdel=as.data.frame(rohdel)
names(rohdel)=levels(rohs$Sample)
rohdel=as.data.frame(t(rohdel))
rohdel$Sample=row.names(rohdel)
names(rohdel)[1:(ncol(rohdel)-1)]=as.character(intervals[1:(length(intervals)-1)])
rohdelm=melt(rohdel,id.vars='Sample')

rohdel=rohdels(rohn,intervals=intervals)
rohdel=as.data.frame(rohdel)
names(rohdel)=levels(rohn$Sample)
rohdel=as.data.frame(t(rohdel))
rohdel$Sample=row.names(rohdel)
names(rohdel)[1:(ncol(rohdel)-1)]=as.character(intervals[1:(length(intervals)-1)])
rohdel=melt(rohdel,id.vars='Sample')
rohdelm=rbind(rohdel,rohdelm)
rohdelm$Pop=str_split_fixed(rohdelm$Sample,'_',Inf)[,1]
rohdelm$Factor=as.factor(apply(rohdelm[,c(2,4)],1,function(x) paste(c(x[2],as.character(x[1]),'MB'),collapse=' ')))

cols=c(brewer.pal(9,'Oranges')[c(6)],'dodgerblue')


roh=rbind(rohn,rohs)
gw=matrix(nrow=nlevels(roh$Sample),ncol=4)
i=1
for (sample in levels(roh$Sample)){
  dfs=roh[which(roh$Sample==sample),]
  gw[i,1]=sample
  rd=sum(dfs$RS4[which(dfs$RS>0 & dfs$Length>=1e6)])
  rdenom=sum(dfs$RS[which(dfs$RS>0 & dfs$Length>=1e6)])
  # rd=sum(dfs$MODERATERS4[which(dfs$LOWRS>0 & dfs$Length>=1e6)])
  # rdenom=sum(dfs$LOWRS[which(dfs$LOWRS>0 & dfs$Length>=1e6)])
  gw[i,2]=rd/rdenom
  # gw[i,3]=(sum(gtsd[which(gtsd$Impact=="MODERATE" & gtsd$RS>=4),sample]==2,na.rm=T)-rd)/(sum(gtsd[which(gtsd$Impact=='LOW' & !is.na(gtsd$RS)),sample]==2,na.rm=T)-rdenom)
  # gw[i,4]=sum(gtsd[which(gtsd$Impact=="MODERATE" & gtsd$RS>=4),sample]==1,na.rm=T)/sum(gtsd[which(gtsd$Impact=='LOW' & !is.na(gtsd$RS)),sample]==1,na.rm=T)
  gw[i,3]=(sum(gtsd[which(gtsd$RS>=4),sample]==2,na.rm=T)-rd)/(sum(gtsd[which(!is.na(gtsd$RS)),sample]==2,na.rm=T)-rdenom)
  gw[i,4]=sum(gtsd[which(gtsd$RS>=4),sample]==1,na.rm=T)/sum(gtsd[!is.na(gtsd$RS),sample]==1,na.rm=T)
  i=i+1
}
gw=as.data.frame(gw)
names(gw)=c('Sample','ROH','gw_hom','gw_het')
gw[,2:4]=apply(gw[,2:4],2,function(x) as.numeric(as.character(x)))
t.test(as.numeric(gw$ROH), gw$gw_hom,paired=T)
t.test(as.numeric(gw$ROH), gw$gw_het,paired=T)

t.test(as.numeric(gw$ROH[grep('NWR',gw$Sample)]), gw$gw_hom[grep('NWR',gw$Sample)],paired=T)#no diff, p=0.2637
t.test(as.numeric(gw$ROH[grep('SWR',gw$Sample)]), gw$gw_hom[grep('SWR',gw$Sample)],paired=T)#no diff, p=0.4614
t.test(as.numeric(gw$ROH[grep('NWR',gw$Sample)]), gw$gw_het[grep('NWR',gw$Sample)],paired=T)#het is greater, p=2.654e-05
t.test(as.numeric(gw$ROH[grep('SWR',gw$Sample)]), gw$gw_het[grep('SWR',gw$Sample)],paired=T)#het is greater, p=2.939e-08

ggplot(data= rohdelm,aes(x=variable,y=value,col=Pop))+geom_boxplot(color='black')+geom_point(position = position_jitterdodge())+scale_color_manual(values=cols)+facet_grid(.~Pop)+theme(text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(),axis.ticks = element_blank(),legend.position='none')+ylab('Proportion of homozygous\nmutations at conserved sites')+xlab('ROH length (MB)')+geom_hline(data=data.frame(yintercept=c(mean(gw$gw_hom[grep('NWR',gw$Sample)]),mean(gw$gw_hom[grep('SWR',gw$Sample)])),Pop=factor(c('NWR','SWR')),lty=3),aes(yintercept = yintercept),linetype = "dashed")+geom_hline(data=data.frame(yintercept=c(mean(gw$gw_het[grep('NWR',gw$Sample)]),mean(gw$gw_het[grep('SWR',gw$Sample)])),Pop=factor(c('NWR','SWR'))),aes(yintercept = yintercept),linetype = "twodash")

