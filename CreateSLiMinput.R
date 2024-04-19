library(fields)
library(RColorBrewer)

nwridx=c(7,9:15)#NWR5763 has no cell line.

#read in all variants in the genome
gtsd=read.table('AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_GERP_Ann_2BR_DerGT.txt',header=T,colClasses=c('character','integer','character','character','numeric','character',rep('integer',22),rep('character',3),'numeric'),sep='\t')

gtsd=gtsd[complete.cases(gtsd[,7:28]),]
gtsd=gtsd[substr(gtsd$Chr, 1, 1)=='J',]
gtsd=gtsd[which(!is.na(gtsd$DicBicMic) & !is.na(gtsd$DicBicMin)),]

cols=c(brewer.pal(9,'Oranges')[c(6)],'dodgerblue')

gtsd$NWRderF=rowSums(gtsd[,nwridx])/(rowSums(!is.na(gtsd[,nwridx]))*2)
gtsd$SWRderF=rowSums(gtsd[,16:28])/(rowSums(!is.na(gtsd[,16:28]))*2)
gtsd$derF=rowSums(gtsd[,7:28])/(rowSums(!is.na(gtsd[,7:28]))*2)

#excluding these lines because we don't want very negative RS values
gtsd$RS[which(gtsd$RS<0)]=0
gtsd$phyloP[which(gtsd$phyloP<0)]=0

busco=read.table('/home/centos/USS/conservation-genetics/nwr/Tables/buscogenes.txt')
busco$busco='Y'
busco$V1=as.character(busco$V1)
gtsd=merge(gtsd,busco,by.x='GeneName',by.y='V1',all.x=T)
gtsd$busco[is.na(gtsd$busco)]='N'
gtsd=gtsd[,c(2:31,1,32:38)]

#choose allele function to go on each haplotype (for heterozygous sites)
chooseAllele=function(x){
  if(x==1){
    g1=sample(0:1,1)
    g2=(0:1)[which(0:1!=g1)]
  }
  else {
    g1=sum(x==2)
    g2=sum(x==2)
  }
  return(c(g1,g2))
}


#number of individuals (N*2)
Ne=16

#get the genome bedfile
bed=read.table('GCA_000283155.1_CerSimSim1.0_genomic.bed',sep='\t')
bigscaf=as.character(bed$V1[bed$V3>=1e6])
chromo=gtsd[which(gtsd$Chr %in% bigscaf & !is.na(gtsd$RS)),]
chromo$rowID=which(gtsd$Chr %in% bigscaf & !is.na(gtsd$RS))
bigscaf=bigscaf[which(bigscaf %in% chromo$Chr)]
bigscaf=bigscaf[order(bigscaf)]
totlgth=sum(bed$V3[which(bed$V1 %in% chromo$Chr)])
chrints=0
for (i in bigscaf){
  brk=bed$V3[which(bed$V1==i)]
  brk=chrints[length(chrints)]+brk
  chrints=c(chrints,brk-1,brk)
}
recints=NULL
for (i in 2:length(chrints)){
  recints=c(recints,ifelse(chrints[i]==chrints[i-1]+1,0.5,1e-8))
}
#for SLiM script initializeRecombinationRate:
paste(chrints[2:length(chrints)],collapse=',')
paste(recints,collapse=',')

muttypes=c('m1','m2','m3','m4')

s=c(-0.002,-0.001,-0.0001,0) #from Henn 2016
#s=c(-0.072,-0.032,-0.024,0) #from Pieschl 2018, not Henn 2016
names(s)=muttypes
h=c(0.0330204,0.0619497,0.292893,0.5) #from Henn 2016 (h=0 if completely recessive, h=0.5 if additive)
names(h)=muttypes

chromo$s=NA
chromo$h=NA
chromo$m=NA
chromo$m[which(chromo$RS>=-2 & chromo$RS<2)]='m4'
chromo$m[which(chromo$RS>=2 & chromo$RS<4)]='m3'
chromo$m[which(chromo$RS>=4 & chromo$RS<5.8)]='m2'
chromo$m[which(chromo$RS>=5.8)]='m1'
chromo=chromo[!is.na(chromo$m),]
chromo$s=s[chromo$m]
chromo$h=h[chromo$m]

#####Fitness#######
#estimate total fitness (including fixed SNPs)
#multiplicative
fitness=as.data.frame((1+s['m1']*h['m1'])^colSums(chromo[which(chromo$m=='m1'),7:28]==1)*(1+s['m2']*h['m2'])^colSums(chromo[which(chromo$m=='m2'),7:28]==1)*(1+s['m3']*h['m3'])^colSums(chromo[which(chromo$m=='m3'),7:28]==1)*(1+s['m1'])^colSums(chromo[which(chromo$m=='m1'),7:28]==2)*(1+s['m2'])^colSums(chromo[which(chromo$m=='m2'),7:28]==2)*(1+s['m3'])^colSums(chromo[which(chromo$m=='m3'),7:28]==2))
names(fitness)='fitness'

#additive
fitness$RL=(s['m1']*h['m1'])*colSums(chromo[which(chromo$m=='m1'),7:28]==1)+(s['m2']*h['m2'])*colSums(chromo[which(chromo$m=='m2'),7:28]==1)+(s['m3']*h['m3'])*colSums(chromo[which(chromo$m=='m3'),7:28]==1)+(s['m1'])*colSums(chromo[which(chromo$m=='m1'),7:28]==2)+(s['m2'])*colSums(chromo[which(chromo$m=='m2'),7:28]==2)+(s['m3'])*colSums(chromo[which(chromo$m=='m3'),7:28]==2)

fitness$pop=factor(c(rep('NWR',9),rep('SWR',13)),levels=c('NWR','SWR'))
swrfitness=median(fitness$fitness[which(fitness$pop=='SWR')])
median(fitness$fitness[which(fitness$pop=='NWR')])
swrRL=median(fitness$RL[which(fitness$pop=='SWR')])
median(fitness$RL[which(fitness$pop=='NWR')])

#fitness reduction from fixed SNPs (to add to w output by slim, which were excluded from simulations)
fixedNWRfitness=(1+s['m1'])^colSums(chromo[which(chromo$m=='m1' & chromo$NWRderF==1),nwridx]==2)*(1+s['m2'])^colSums(chromo[which(chromo$m=='m2' & chromo$NWRderF==1),nwridx]==2)*(1+s['m3'])^colSums(chromo[which(chromo$m=='m3' & chromo$NWRderF==1),nwridx]==2)#all values should be the same
fixedNWRfitness=unique(fixedNWRfitness)
fixedNWR_RL=(s['m1'])*sum(chromo$m=='m1' & chromo$NWRderF==1)+(s['m2'])*sum(chromo$m=='m2' & chromo$NWRderF==1)+(s['m3'])*sum(chromo$m=='m3' & chromo$NWRderF==1)

#get other fixed variants (to add to SLiM output)
nwrfm1=sum(chromo$m=='m1' & chromo$NWRderF==1)
nwrfm2=sum(chromo$m=='m2' & chromo$NWRderF==1)
nwrfm3=sum(chromo$m=='m3' & chromo$NWRderF==1)
nwrfload=sum(chromo$m!='m4' & chromo$NWRderF==1)

#SWR medians
swrm1=median(colSums(chromo[which(chromo$m=='m1'),16:28]))
swrm2=median(colSums(chromo[which(chromo$m=='m2'),16:28]))
swrm3=median(colSums(chromo[which(chromo$m=='m3'),16:28]))
swrload=median(colSums(chromo[which(chromo$m!='m4'),16:28]))
swrm1hom=median(colSums(chromo[which(chromo$m=='m1'),16:28]==2))
swrm2hom=median(colSums(chromo[which(chromo$m=='m2'),16:28]==2))
swrm3hom=median(colSums(chromo[which(chromo$m=='m3'),16:28]==2))
swrloadhom=median(colSums(chromo[which(chromo$m!='m4'),16:28]==2))

swrmedian=data.frame(var=c("load","m2","m1hom","m1","w","RL","m2hom","m3","homload","m3hom"),value=c(swrload,swrm2,swrm1hom,swrm1,swrfitness,swrRL,swrm2hom,swrm3,swrloadhom,swrm3hom))
write.table(swrmedian,'/home/centos/USS/conservation-genetics/nwr/slim/Reintroductions/SWRmedianLoad_Hess_Feb.txt',row.names=F,quote=F,sep='\t')

nwrfixedload=data.frame(var=c("nwrfload","nwrfm1","nwrfm2","nwrfm3","fixedNWRfitness","fixedNWR_RL"),value=c(nwrfload,nwrfm1,nwrfm2,nwrfm3,fixedNWRfitness,fixedNWR_RL))
write.table(nwrfixedload,'/home/centos/USS/conservation-genetics/nwr/slim/Reintroductions/NWRfixedLoad_Hess_Feb.txt',row.names=F,quote=F,sep='\t')

chromo=chromo[which(chromo$NWRderF>0 & chromo$NWRderF<1),]
idx=sample(which(chromo$m=='m4'),10000)
chromo=chromo[c(which(chromo$m!='m4'),idx),]
chromo=chromo[order(chromo$Chr,chromo$Pos),]

chromo$SNPNum=format(0:(nrow(chromo)-1),scientific=F)
#assign PosID as one continuous segment
chromo$PosID=apply(chromo,1,function(x) chrints[which(bigscaf==as.character(x[1]))*2-1]+as.numeric(x[2]))
chromo$PosID=format(chromo$PosID,scientific=F)

SNPNumID=which(names(chromo)=='SNPNum')
rowIDn=which(names(chromo)=='rowID')
PosID=which(names(chromo)=='PosID')
sID=which(names(chromo)=='s')
hID=which(names(chromo)=='h')
mID=which(names(chromo)=='m')
NWRderFID=which(names(chromo)=='NWRderF')
derFID=which(names(chromo)=='derF')

#choose alleles on haplotypes
haplotypes=apply(chromo[,nwridx],2,function(x) sapply(x,chooseAllele))
haplotypes =cbind(haplotypes[seq(1,nrow(haplotypes),2),],haplotypes[seq(2,nrow(haplotypes),2),])
colorder=1:length(nwridx)
names(colorder)=colnames(haplotypes)[1:length(nwridx)]
haplotypes=haplotypes[,order(colorder[colnames(haplotypes)])]
genotype=apply(haplotypes,2,function(x) format(which(x==1)-1,scientific=F))
chromo$s=format(chromo$s, scientific = FALSE)

#write outfiles
filename="NWR_RS2MODHInuet_AllChr_polym_ssm4_Henns_n8.txt"
write(paste0("#OUT: 1 1 A\nVersion: 4\nPopulations:\np1 ",Ne/2," H\nMutations:"),file=filename)
apply(chromo,1,function(x) write(paste0(c(x[SNPNumID],x[rowIDn],x[mID],x[PosID],x[sID],x[hID],' p1 1',x[NWRderFID]),collapse=' '),file=filename,append=T,ncol=9))

write("Individuals:",file= filename,append=T)
for (i in 0:(Ne/2-1)){
#  write(paste(c('p1:i',i,' H p1:',i*2,' p1:',i*2+1),collapse=''),file=filename,append=T)
  write(paste(c('p1:i',i,' H p1:',i*2,' p1:',i*2+1,' 0'),collapse=''),file=filename,append=T)
}

write("Genomes:",file= filename,append=T)
for (i in 1:length(genotype)){
  write(paste(c(paste(c('p1:',i-1,' A'),collapse=''),genotype[[i]]),collapse=' '),file= filename,append=T,ncol=max(rowSums(haplotypes)))
}
