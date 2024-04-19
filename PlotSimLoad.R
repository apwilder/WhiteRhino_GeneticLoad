library(stringr)
library(reshape)
library(RColorBrewer)
library(colorRamps)
library(ggplot2)
library(data.table)
library(datawizard)


cols=c(brewer.pal(9,'Oranges')[c(6)],'dodgerblue')
names(cols)=c('NWR','SWR')

grthcols=c("#8B0069","#A4473D","#A57E00","#86B13E","#53DCAB","#B0F4FA")
names(grthcols)=0:5

nwrfixedload=read.table('NWRfixedLoad_Hess_Feb.txt',sep='\t',header=T)
nwrfload=nwrfixedload$value[which(nwrfixedload$var=='nwrfload')]
nwrfm1=nwrfixedload$value[which(nwrfixedload$var=='nwrfm1')]
nwrfm2=nwrfixedload$value[which(nwrfixedload$var=='nwrfm2')]
nwrfm3=nwrfixedload$value[which(nwrfixedload$var=='nwrfm3')]
fixedNWRfitness=nwrfixedload$value[which(nwrfixedload$var=='fixedNWRfitness')]

Loadmat=NULL
Growth=matrix(nrow=11,ncol=7)
Growth[,1]=1:11
for (grth in c(0:5)){#c(0:5)
Load=read.table(paste(c('/home/centos/USS/conservation-genetics/nwr/slim/Reintroductions/NWR_RS2MODHInuet_AllChr_polym_ssm4_Henns_1.',grth,'grth_supp_SimIndivLoad.txt'),collapse=''),header=T,sep='\t')
Load2=read.table(paste(c('/home/centos/USS/conservation-genetics/nwr/slim/Reintroductions/Nov_files/NWR_RS2MODHInuet_AllChr_polym_ssm4_Henns_1.',grth,'grth_supp_SimIndivLoad.txt.gz'),collapse=''),header=T,sep='\t')
names(Load2)[2]='load_1'
nl2=names(Load2)
names(Load2)=sapply(nl2,function(x) paste(str_split_fixed(x,'_',2)[1,1],as.integer(str_split_fixed(x,'_',2)[1,2])+40,sep='_'))
Load=cbind(Load,Load2[,2:ncol(Load2)])
#rows are genomes (the # increases every gen with pop growth), cols are var types for different reps
Load[,grep('m1_',names(Load))]=Load[,grep('m1_',names(Load))]-Load[,grep('m1hom_',names(Load))]*2
Load[,grep('m2_',names(Load))]=Load[,grep('m2_',names(Load))]-Load[,grep('m2hom_',names(Load))]*2
Load[,grep('m3_',names(Load))]=Load[,grep('m3_',names(Load))]-Load[,grep('m3hom_',names(Load))]*2

Growth[,grth+2]=table(Load$Gen_1)

for (r in 1:50){
Load[,paste0('RL_',r)]=(s['m1'])*Load[,paste0('m1hom_',r)]+(s['m2'])*Load[,paste0('m2hom_',r)]+(s['m3'])*Load[,paste0('m3hom_',r)]+(s['m1']*h['m1'])*(Load[,paste0('m1_',r)])+(s['m2']*h['m2'])*(Load[,paste0('m2_',r)])+(s['m3']*h['m3'])*(Load[,paste0('m3_',r)])
}
fixedNWR_RL=(s['m1'])*nwrfm1+(s['m2'])*nwrfm2+(s['m3'])*nwrfm3

#add fixed fitness to w from slim
Load[,grep('w_',names(Load))]=Load[,grep('w_',names(Load))]*fixedNWRfitness
Load[,grep('RL_',names(Load))]=Load[,grep('RL_',names(Load))]+fixedNWR_RL
names(Load)[2]='load_1'
Load=as.data.frame(apply(Load,2,as.numeric))

loadmat=aggregate(Load,by=list(Load$Gen_1),median)#median across individuals
loadmat=reshape2::melt(loadmat[,2:ncol(loadmat)],id.vars=c('Gen_1'))
loadmat$Gen_1=as.integer(loadmat$Gen_1)
loadmat$rep=str_split_fixed(loadmat$variable,"_",Inf)[,2]
loadmat$var=str_split_fixed(loadmat$variable,"_",Inf)[,1]
loadmat$value[which(loadmat$var=='m1')]=loadmat$value[which(loadmat$var=='m1')]+nwrfm1*2
loadmat$value[which(loadmat$var=='m2')]=loadmat$value[which(loadmat$var=='m2')]+nwrfm2*2
loadmat$value[which(loadmat$var=='m3')]=loadmat$value[which(loadmat$var=='m3')]+nwrfm3*2
loadmat$value[which(loadmat$var=='m1hom')]=loadmat$value[which(loadmat$var=='m1hom')]+nwrfm1
loadmat$value[which(loadmat$var=='m2hom')]=loadmat$value[which(loadmat$var=='m2hom')]+nwrfm2
loadmat$value[which(loadmat$var=='m3hom')]=loadmat$value[which(loadmat$var=='m3hom')]+nwrfm3
loadmat$value[which(loadmat$var=='homload')]=loadmat$value[which(loadmat$var=='homload')]+nwrfload
loadmat$value[which(loadmat$var=='load')]=loadmat$value[which(loadmat$var=='load')]+nwrfload*2

loadmat$growth=grth
Loadmat=rbind(Loadmat,loadmat)
}
head(Loadmat)
Growth=as.data.frame(Growth)
names(Growth)=c('Gen_1',0:5)
Growth=reshape2::melt(Growth,id.vars='Gen_1')
Growth$growth=Growth$variable
Growth$var='growth'
Growth$rep=1
Loadmat=rbind(Loadmat,Growth)
Loadmat$Gen_1=Loadmat$Gen_1-1

Loadmat$var=factor(Loadmat$var,levels=c("m1hom","m2hom","m3hom","homload","m1","m2","m3","load",'m4','m4comp',"w",'RL','growth'),labels=c('Homozygous severe','Homozygous moderate','Homozygous mild','Homozygous load','Heterozygous severe','Heterozygous moderate','Heterozygous mild','Heterozygous load','Neutral','Neutral (ss)','Fitness','Realized load','Population size'))

swrmedian=read.table('/home/centos/USS/conservation-genetics/nwr/slim/Reintroductions/SWRmedianLoad_Hess_Feb.txt',header=T,sep='\t')
swrmedian$value[grep('m2$',swrmedian$var)]=swrmedian$value[grep('m2$',swrmedian$var)]-swrmedian$value[grep('m2hom',swrmedian$var)]*2
swrmedian$value[grep('m3$',swrmedian$var)]=swrmedian$value[grep('m3$',swrmedian$var)]-swrmedian$value[grep('m3hom',swrmedian$var)]*2
swrmedian$value[grep('m1$',swrmedian$var)]=swrmedian$value[grep('m1$',swrmedian$var)]-swrmedian$value[grep('m1hom',swrmedian$var)]*2

swrmedian$var=factor(swrmedian$var,levels=c("m1","m2","m3","load","m1hom","m2hom","m3hom","homload","w","RL"),labels=c('Heterozygous severe','Heterozygous moderate','Heterozygous mild','Heterozygous load','Homozygous severe','Homozygous moderate','Homozygous mild','Homozygous load','Fitness','Realized load'))

#relative to SWR median
Loadmat$value[which(Loadmat$var=="Fitness")]=Loadmat$value[which(Loadmat$var=="Fitness")]/swrmedian$value[which(swrmedian$var=="Fitness")]
swrmedian$value[which(swrmedian$var=="Fitness")]=1

gg=ggplot(data=Loadmat[which(Loadmat$var!="Neutral (ss)"),],aes(x=Gen_1,y=value,group=interaction(as.factor(rep),factor(growth,levels=5:0)),color = as.factor(growth)))+geom_line(alpha=0.2,lwd=0.5)+geom_smooth(lwd=0.75,se=F,aes(x=Gen_1,y=value,group=factor(growth,levels=5:0)))+scale_color_manual(values=grthcols, name="Growth (lambda)",labels=seq(1.0,1.5,0.1))+scale_x_continuous(name="Generation", breaks=seq(2,10,2))+facet_wrap(.~var,scales="free")+geom_hline(data=swrmedian,aes(yintercept=value),lty=2,lwd=0.9)#,nrow=3,ncol=4)#

gg+theme(text = element_text(size=10),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black",size=0.5),panel.spacing.y=unit(0.05,'in'),panel.spacing.x=unit(0.05,'lines'),strip.text.y = element_text(size = 6, angle = 0),legend.title=element_text(size=9),legend.position=c(0.9, 0.15),legend.key=element_blank(),legend.key.height = unit(0.4, "cm"),axis.title.y=element_blank())


