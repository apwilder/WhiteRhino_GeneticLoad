library("rtracklayer")
library(stringr)

bed=read.table('AllWR_CerSimSim1.0_allFilt_QD5_FS32_SOR3_MQ28_MQRS2_RPRS2.bed',skip=1)
bed$V1=str_split_fixed(bed$V1,'.1',Inf)[,1]
grObject <- GRanges(seqnames=bed$V1, ranges=IRanges(start=bed$V2, end=bed$V3))
chainObject <- import.chain("cerSim1ToHg19.over.chain")
results <- as.data.frame(liftOver(grObject, chainObject))
bed$idx=1:nrow(bed)
bedmg=merge(bed,results,by.x='idx',by.y='group',all.x=T)#all sites in original bed
numAnnotatedSites=sum(!is.na(bedmg$start))
pctAnnotatedSites=numAnnotatedSites/nrow(bedmg)

newbed=bedmg[!is.na(bedmg$start),]#hg19 coords w/out NA (subset of vars that lifted over)

grObject <- GRanges(seqnames=newbed$seqnames, ranges=IRanges(start=newbed$start, end=newbed$end))
chainObject <- import.chain("hg19ToHg18.over.chain")
results2 <- as.data.frame(liftOver(grObject, chainObject))# subset of newbed that lifted over to hg18
newbed$idx2=1:nrow(newbed)#create index2 to refer to results2
newbedmg=merge(newbed,results2,by.x='idx2',by.y='group',all.x=T)
origbedmg=merge(bed,newbedmg[,c('V1','V2','V3','seqnames.y','start.y','end.y')],by=c('V1','V2','V3'),all.x=T)
names(origbedmg)[5:7]=c('seqnames','start','end')
numAnnotatedSites2=sum(!is.na(origbedmg $start))
pctAnnotatedSites2=numAnnotatedSites/nrow(origbedmg)

setwd('hg18/')
source_lines <- function(file, lines){
read.table(textConnection(readLines(file)[lines]),sep='\t')
}
gerplengths=read.table('GERPrateFileLines.txt')
for (i in c(2:22,"Y","X")){
	gerplength= gerplengths[which(gerplengths$V1==i),2]
	bedss=origbedmg[which(origbedmg$seqnames==paste0('chr',i)),]
	gerp=source_lines(paste(c('chr',i,'_GERP.rates.gz'),collapse=''),bedss$start[which(bedss$start<=gerplength)])
	gerp=rbind(gerp,as.data.frame(matrix(data=0,nrow=length(which(bedss$start>gerplength)),ncol=2)))
	gerps=cbind(origbedmg[which(origbedmg$seqnames==paste0('chr',i)),c('V1','V2','start')],gerp)
	names(gerps)[c(1,2,4,5)]=c('Chr','Pos','neutRate','RS')
	write.table(file=paste(c('chr',i,'snp_GERPscores.txt'),collapse='_'),gerps[,c(1,2,4,5)],row.names=F,sep='\t',quote=F)
}
