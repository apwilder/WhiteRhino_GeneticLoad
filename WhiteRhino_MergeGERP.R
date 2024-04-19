
swr=read.table('SWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3_plink_AF.txt',colClasses=c('factor','character','integer','character','character','numeric','integer')) 
nwr=read.table('NWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3_plink_AF.txt',colClasses=c('factor','character','integer','character','character','numeric','integer'))

swr=swr[,2:7]
nwr=nwr[,2:7]
names(swr)=c('Chr','Pos','SWR_Ref','SWR_Alt','SWR_AltF','SWR_n')
names(nwr)=c('Chr','Pos','NWR_Ref','NWR_Alt','NWR_AltF','NWR_n')

varinfo=read.table('AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3_variantInfo.txt',colClasses=c('character','integer','character','character'))
names(varinfo)=c('Chr','Pos','Ref','Alt')

wr=merge(swr,nwr,by=c('Chr','Pos'),all=T)
wr=merge(wr,varinfo,by=c('Chr','Pos'),all.x=T)

gerp=read.table('AllChr_snp_GERPscores.txt',colClasses=c('character','integer','numeric','numeric')) #created with cat chr_*_snp_GERPscores.txt > AllChr_snp_GERPscores.txt, from GetGERPS.R
names(gerp)=c('Chr','Pos','neutRate','RS')
gerp$Chr=paste0(gerp$Chr,'.1')
wrg=merge(wrbc[,1:13],gerp,by=c('Chr','Pos'),all.x=T)

write.table(wrg,file='AllWR_CerSimSim1.0_AF_snpGERPs.txt',col.names=T,row.names=F,quote=F,sep='\t')

##Next polarize ancestral/derived alleles!
br=read.table('DicBic_2CerSimSim_Variants.txt',colClasses=c('character','integer','character','character','character','numeric'))
names(br)=c('Chr','Pos','ID','BRRef','BRAlt','BRqual')
br=br[,c(1,2,4:6)] 
wrb=merge(wr,br,by=c('Chr','Pos'),all.x=T) 
wrb=wrb[which(wrb$Chr!="Chr.1"),] 
write.table(wrb,file='AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GQ20_mm3_AF_snpGERPs.txt',col.names=T,row.names=F,quote=F,sep='\t')

brcov=read.table('DicBicMic_2CerSimSim_mappedSites_WRsnps.txt')
names(brcov)=c('Chr','Pos','BRcov')
wrbc=merge(wrb,brcov,by=c('Chr','Pos'),all.x=T)


write.table(wrbc,file='AllWR_CerSimSim1.0_raw_AF_snpGERPs_BRalleles.txt',col.names=T,row.names=F,quote=F,sep='\t')

#polarize allele frequencies on Alt from ref allele
wrbc$SWR_AltF[which(wrbc$SWR_Ref!=wrbc$Ref)]=1-wrbc$SWR_AltF[which(wrbc$SWR_Ref!=wrbc$Ref)] #polarize AF to reference genome's alt allele
wrbc$NWR_AltF[which(wrbc$NWR_Ref!=wrbc$Ref)]=1-wrbc$NWR_AltF[which(wrbc$NWR_Ref!=wrbc$Ref)] #polarize AF to reference genome's alt allele
wrbc=wrbc[,c(1,2,5,6,9:18)]

wrbc$BRallele=rep(NA,nrow(wrbc))

wrbc$BRallele[which(wrbc$BRcov==1 & is.na(wrbc$BRRef))]=as.character(wrbc$Ref[which(wrbc$BRcov==1 & is.na(wrbc$BRRef))]) #sites with BR coverage and where BR doesn't differ from WR Ref (NA means it's not in SNP list), set BR allele as WR Ref
wrbc$BRallele[which(wrbc$BRcov==1 & !is.na(wrbc$BRRef))]=as.character(wrbc$BRAlt[which(wrbc$BRcov==1 & !is.na(wrbc$BRRef))]) #sites with BR coverage that differ from WR Ref (not NA means it's in SNP list), BR allele (allele in BR genome) is Alt

wrbc$BRcov[which(is.na(wrbc$BRcov))]=0 #change NA's to no coverage
wrbc$BRallele[which(wrbc$BRcov>1)]=NA #change BR alleles to NA for sites with BR coverage >1 (where multiple BR maps to single WR site)
wrbc$BRallele[which(wrbc$BRcov==0)]=NA #change BR alleles to NA for sites with no coverage

#length(which(is.na(wrbc$BRallele) & (wrbc$BRcov>1 | is.na(wrbc$Ref) | wrbc$BRcov==0)))==sum(is.na(wrbc$BRallele)) TRUE
#length(which((wrbc$BRcov>1 | is.na(wrbc$Ref) | wrbc$BRcov==0)))==sum(is.na(wrbc$BRallele)) TRUE

#number of deleterious derived alleles in SWR
length(which(wrbc$SWR_AltF>0 & wrbc$Alt!=wrbc$BRallele & wrbc$RS>4))
length(which(wrbc$SWR_AltF>0 & wrbc$Alt!=wrbc$BRallele & wrbc$RS>4))/length(which(wrbc$SWR_AltF>0 & wrbc$Alt!=wrbc$BRallele))

length(which(wrbc$NWR_AltF>0 & wrbc$Alt!=wrbc$BRallele & wrbc$RS>4))
length(which(wrbc$NWR_AltF>0 & wrbc$Alt!=wrbc$BRallele & wrbc$RS>4))/length(which(wrbc$NWR_AltF>0 & wrbc$Alt!=wrbc$BRallele))

#calculate derived frequencies
wrbc$SWRder=rep(NA,nrow(wrbc))
wrbc$NWRder=rep(NA,nrow(wrbc))
wrbc$SWRder[which(wrbc$BRallele==wrbc$Alt)]=1-wrbc$SWR_AltF[which(wrbc$BRallele==wrbc$Alt)]
wrbc$SWRder[which(wrbc$BRallele==wrbc$Ref)]=wrbc$SWR_AltF[which(wrbc$BRallele==wrbc$Ref)]
wrbc$NWRder[which(wrbc$BRallele==wrbc$Alt)]=1-wrbc$NWR_AltF[which(wrbc$BRallele==wrbc$Alt)]
wrbc$NWRder[which(wrbc$BRallele==wrbc$Ref)]=wrbc$NWR_AltF[which(wrbc$BRallele==wrbc$Ref)]


write.table(wrbc,file='AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_derFreq_GERP.txt.gz',col.names=T,row.names=F,quote=F,sep='\t')

gt=read.table('AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GT.GT.FORMAT.gz',sep='\t',header=T)
#wrbc=read.table('AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_derFreq_GERP.txt.gz',header=T,colClasses=c('character','integer','numeric','integer','numeric','integer','character','character','numeric','numeric','character','character','numeric','integer','character','numeric','numeric'))

gt=merge(gt,wrbc[,c('Chr','Pos','Ref','Alt','RS','BRallele','SWRder','NWRder')],by.x=c('CHROM','POS'),by.y=c('Chr','Pos'))
names(gt)[1:2]=c('Chr','Pos')
map = setNames(c(0,1,2,NA),c("0/0", "0/1", "1/1","./."))#0=2 refs, 2=0 ref alleles
gts=apply(gt[,3:24],2,function(x) map[as.character(x)])
gts=cbind(gt[,c('Chr','Pos','Ref','Alt','RS','BRallele','SWRder','NWRder')],gts)
gts$Ref=as.character(gts$Ref)
gts$BRallele =as.character(gts$BRallele)
gts$Alt =as.character(gts$Alt)

#define derived allele
gts$Der=rep(NA,nrow(gts))
gts$Der[which(gts$Alt==gts$BRallele)]=gts$Ref[which(gts$Alt==gts$BRallele)]
gts$Der[which(gts$Ref==gts$BRallele)]=gts$Alt[which(gts$Ref==gts$BRallele)]

wridx=9:30
nwridx=9:17
swridx=18:30

#polarize on derived allele
gtsd=gts[which(!is.na(gts$Der)),]
gtsd[which(gtsd$Der==gtsd$Ref), wridx]=2-gtsd[which(gtsd$Der==gtsd$Ref), wridx]
derAlleles=colSums(gtsd[, wridx],na.rm=T)
deletDerAlleles=colSums(gtsd[which(gtsd$RS>=4), wridx],na.rm=T)
deletDerAlleles/derAlleles

write.table(gtsd2,'AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2_GERP_derGTs.txt',row.names=F,sep='\t',quote=F)
