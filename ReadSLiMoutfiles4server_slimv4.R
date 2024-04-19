library(stringr)
args <- commandArgs(trailingOnly = TRUE)
wd=args[1] #e.g. '/home/centos/USS/conservation-genetics/nwr/slim/Reintroductions'
mainfile=args[2] #e.g. '/home/centos/USS/conservation-genetics/nwr/slim/NWR_RS2MODHInuet_AllChr_polym_ssm4_Henns_n8.txt' /RS2MODHInuet_AllChr_polym_ssm4_Pieschls.txt'
rootfile=args[3] # e.g.'RS2MODHInuet_AllChr_polym_ssm4_Pieschls_10grth_gen'
seedindex=as.numeric(args[4]) # 8 if slimout.txt files are split on "_", index of seed number for replicate
outbase=args[5]#'RS2MODHInuet_AllChr_polym_ssm4_Pieschls_10grth'
supp=args[6]#1 if '_supp_'
founders=as.numeric(args[7])#8 number of founders

setwd(wd)

muttypes=c('m1','m2','m3','m4')
s=c(-0.002,-0.001,-0.0001,0) #from Henn 2016
names(s)=muttypes
h=c(0.0330204,0.0619497,0.292893,0.5) #from Henn 2016 (h=0 if completely recessive, h=0.5 if additive)
names(h)=muttypes

supp=ifelse(supp==1,"_supp_","_nosupp_")
seedlist=list.files(path=".",pattern=rootfile)#,pattern='_slimout.txt')
seedlist=seedlist[grep(supp,seedlist)]
seedlist=unique(str_split_fixed(seedlist,"_",Inf)[,seedindex])

#header=6 lines (main =5)
#mutations (298914)
#Indiv header (1 line)
#founders+Ne lines (main=8)
#Genomes header
#founders*2+Ne*2 (main=16)

Ne=as.integer(scan(mainfile, skip = 3,nlines= 1,what=c('character','integer','character'))[2])
mutations=scan(mainfile, skip = 5,nlines= length(readLines(mainfile)) - (Ne*3+3+4),sep='\n',what='')
nmuts=length(mutations)
m1=matrix(nrow=length(mutations),ncol=9)
for (i in 1:length(mutations)){
  mi=unlist(strsplit(mutations[i]," +"))
  if (length(mi)==9){
    m1[i,]=mi
  } else {
    m1[i,]=mi[2:10]}
}
m1=as.data.frame(m1)
names(m1)=c('Index','pID','muttype','Pos','s','h','Pop','gen','Frq1')

genomes=scan(mainfile, skip = length(readLines(mainfile)) - (Ne*2),sep='\n',what='')
y <- strsplit(genomes, "[[:space:]]+")
names(y) <- sapply(y, `[[`, 1)

mfreq=matrix(ncol=Ne,nrow=nrow(m1),data=0)
row.names(mfreq)=as.character(m1$Index)
indivs=scan(mainfile, skip = length(readLines(mainfile)) - (Ne*3+1),nlines=Ne,sep='\n',what='')
for (i in 1:Ne){
  haps=strsplit(indivs[i],' ')[[1]][3:4]
  muts=c(y[[haps[1]]][3:length(y[[haps[1]]])],y[[haps[2]]][3:length(y[[haps[2]]])])
  mutstable=table(muts)
  mfreq[names(mutstable),i]=mutstable#64250:64500
}

allsnps=merge(m1,as.data.frame(mfreq),by.x='Index',by.y=0)
allsnps$Frq1=as.numeric(as.character(allsnps$Frq1))
names(allsnps)[(ncol(allsnps)-Ne+1):ncol(allsnps)]=paste0('Gen1_',0:(Ne-1))

Gen=c(2:11)
for (j in 1:length(seedlist)){
  print(paste0('Replicate ',j))
  #polymorphic loci
  for (f in Gen){
    print(paste0('Generation ',f))
    filename=paste(c(rootfile,f,supp,seedlist[j],'_slimout.txt'),collapse='')
    Ne=as.integer(scan(filename, skip = 4,nlines= 1,what=c('character','integer','character'))[2])
    print(Ne)
    mutations=scan(filename, skip = 6,nlines= nmuts,sep='\n',what='')#length(readLines(filename)) - (Ne*3+founders*3+2+6)was (Ne*3+2+5), added 9*3 because now including FZ lines
    m=matrix(nrow=length(mutations),ncol=9)
    for (i in 1:length(mutations)){
      m[i,]=unlist(strsplit(mutations[i]," "))
    }
    m=as.data.frame(m)
    
    genomes=scan(filename, skip = length(readLines(filename)) - (Ne*2),sep='\n',what='')
    y <- strsplit(genomes, "[[:space:]]+")
    names(y) <- sapply(y, `[[`, 1)
    
    mfreq=matrix(ncol=Ne,nrow=nrow(m),data=0)
    row.names(mfreq)=m$V1
    indivs=scan(filename, skip = length(readLines(filename)) - (Ne*3+founders*2+1),nlines=Ne,sep='\n',what='')
    for (i in 1:Ne){
      haps=strsplit(indivs[i],' ')[[1]][3:4]
      muts=c(y[[haps[1]]][3:length(y[[haps[1]]])],y[[haps[2]]][3:length(y[[haps[2]]])])
      mutstable=table(muts)
      mfreq[names(mutstable),i]=mutstable
    }
    
    mfreq=as.data.frame(mfreq)
    mfreq$pID=as.integer(as.character(m$V2))
    
    allsnps=merge(allsnps,mfreq,all.x=T,by='pID')
    allsnps[is.na(allsnps$V1),(ncol(allsnps)-Ne+1):ncol(allsnps)]=0
    names(allsnps)[(ncol(allsnps)-Ne+1):ncol(allsnps)]=paste0('Gen',f,'_',0:(Ne-1))
    allsnps$Frq=rowSums(allsnps[,(ncol(allsnps)-Ne+1):ncol(allsnps)]/(Ne*2))
    names(allsnps)[which(names(allsnps)=='Frq')]=paste(c('Frq',f,'_',j),collapse='')
  }
  
  load=as.data.frame(colSums(allsnps[which(allsnps$muttype %in% c('m1','m2','m3')),grep('Gen',names(allsnps))]))
  names(load)='load'
  if(j==1){
    load$Gen=str_split_fixed(row.names(load),'_',Inf)[,1]
    load$Gen=sub('Gen','',load$Gen)}
  # boxplot(load$load~factor(load$Gen))
  load$m1=colSums(allsnps[which(allsnps$muttype=='m1'),grep('Gen',names(allsnps))])
  load$m2=colSums(allsnps[which(allsnps$muttype=='m2'),grep('Gen',names(allsnps))])
  load$m3=colSums(allsnps[which(allsnps$muttype=='m3'),grep('Gen',names(allsnps))])
  load$m4=colSums(allsnps[which(allsnps$muttype=='m4'),grep('Gen',names(allsnps))])
  
  m1AF=table(allsnps$Frq1[which(allsnps$muttype=='m1')])
  neut=NULL
  for (frq in names(m1AF)){
    neut=c(neut,sample(which(allsnps$muttype=='m4' & allsnps$Frq1==frq),m1AF[frq])) 
  }
  load$m4comp=colSums(allsnps[neut,grep('Gen',names(allsnps))])
  load$m1hom=colSums(allsnps[which(allsnps$muttype=='m1'),grep('Gen',names(allsnps))]==2)
  load$m2hom=colSums(allsnps[which(allsnps$muttype=='m2'),grep('Gen',names(allsnps))]==2)
  load$m3hom=colSums(allsnps[which(allsnps$muttype=='m3'),grep('Gen',names(allsnps))]==2)
  load$homload=colSums(allsnps[which(allsnps$muttype %in% c('m1','m2','m3')),grep('Gen',names(allsnps))]==2)
  load$w=(1+s['m1']*h['m1'])^(load$m1-load$m1hom*2)*(1+s['m2']*h['m2'])^(load$m2-load$m2hom*2)*(1+s['m3']*h['m3'])^(load$m3-load$m3hom*2)*(1+s['m1'])^(load$m1hom)*(1+s['m2'])^(load$m2hom)*(1+s['m3'])^(load$m3hom)
    
  vars=11
  names(load)[(ncol(load)-vars+1):ncol(load)]=paste(names(load)[(ncol(load)-vars+1):ncol(load)],j,sep='_')
  
  allsnps=allsnps[,grep('Gen',names(allsnps),invert=T)]
  if(j==1){
    Load=load
  } else{
    Load=cbind(Load,rbind(setNames(Load[1:founders,c(1,3:(vars+1))],names(load)),load))}
}

Load=Load[,c(2,1,3:ncol(Load))]


write.table(Load,file=paste0(outbase,supp,'SimIndivLoad.txt'),sep='\t',quote=F)
write.table(allsnps,file=paste0(outbase,supp,'SimAF.txt'),row.names=F,sep='\t',quote=F)

