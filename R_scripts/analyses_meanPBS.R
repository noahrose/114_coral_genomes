options(stringsAsFactors=F)
library(zoo)
library(ape)
library(pegas)
library(data.table)
library(permute)

setwd('~/Dropbox/GradWork/114_coral_genomes/processed_data')

amil<-read.delim('amil.picard')
sum(amil$PF_HQ_ALIGNED_READS)
mean(amil$PF_HQ_ERROR_RATE)
altref<-read.delim('altref.picard')
sum(altref$PF_HQ_ALIGNED_READS)
mean(altref$PF_HQ_ERROR_RATE)

admix<-read.delim('admix.qopt',header=F,sep=' ')[,-5]
admix<-as.matrix(t(admix))
samps<-read.delim('samps.txt',header=F)[,1]
colnames(admix)<-samps
species<-apply(admix,2,which.max)
names(species)<-samps
meta<-read.csv('meta.csv')
meta$species<-species[match(meta$filename,names(species))]

meta<-meta[order(meta$cryp),]
admix<-admix[,match(meta$filename,colnames(admix))]
admix<-admix[order(meta$cryp[apply(admix,1,which.max)]),]

meta$cryp<-factor(meta$cryp)

idx<-read.delim('all.idxstats.gz',header=T)
csym<-idx[grep('c_sym',idx$chr),]
ccov<-125*rowSums(csym[,3:116])/csym[,2]
csym<-csym[ccov<2*median(ccov)&ccov>0.5*median(ccov),]
cnorm<-125*colSums(csym[,3:116])/sum(csym[,2])

dsym<-idx[grep('d_sym', idx$chr),]
dcov<-125*rowSums(dsym[,3:116])/dsym[,2]
dsym<-dsym[dcov<2*median(dcov)&dcov>0.5*median(dcov),]
dnorm<-125*colSums(dsym[,3:116])/sum(dsym[,2])

propd<-dnorm/(dnorm+cnorm)
meta$propd<-propd[match(meta$filename ,gsub('X','',names(propd)))]

meta$binsym<-'C'
meta$binsym[meta$propd>0.5]<-'D'
meta$pool<-factor(meta$pool)
write.csv(meta,file='meta.csv',row.names=F,quote=F)

glm.out1<-glm(binbl~pool+cryp+propd,data=meta,family='binomial')
glm.out2<-glm(binbl~pool+propd,data=meta,family='binomial')


####################
#GENOMIC ANALYSES          
####################

HA.HC<-read.delim('HA.HC.fst.components.txt',row.names=NULL,header=F,col.names=c('chr','start','stop','num','den','nSites'),na.strings='.')
HA.HD<-read.delim('HA.HD.fst.components.txt',row.names=NULL,header=F,col.names=c('chr','start','stop','num','den','nSites'),na.strings='.')
HA.HE<-read.delim('HA.HE.fst.components.txt',row.names=NULL,header=F,col.names=c('chr','start','stop','num','den','nSites'),na.strings='.')
HC.HD<-read.delim('HC.HD.fst.components.txt',row.names=NULL,header=F,col.names=c('chr','start','stop','num','den','nSites'),na.strings='.')
HC.HE<-read.delim('HC.HE.fst.components.txt',row.names=NULL,header=F,col.names=c('chr','start','stop','num','den','nSites'),na.strings='.')
HD.HE<-read.delim('HD.HE.fst.components.txt',row.names=NULL,header=F,col.names=c('chr','start','stop','num','den','nSites'),na.strings='.')
HA.HC$pos<-paste(HA.HC$chr,HA.HC$start,sep=':')
HA.HD$pos<-paste(HA.HD$chr,HA.HD$start,sep=':')
HA.HE$pos<-paste(HA.HE$chr,HA.HE$start,sep=':')
HC.HD$pos<-paste(HC.HD$chr,HC.HD$start,sep=':')
HC.HE$pos<-paste(HC.HE$chr,HC.HE$start,sep=':')
HD.HE$pos<-paste(HD.HE$chr,HD.HE$start,sep=':')
postab<-table(c(HA.HC$pos,HA.HD$pos,HA.HE$pos,HC.HD$pos,HC.HE$pos,HD.HE$pos))
overlap<-names(postab[postab==6])
HA.HC<-HA.HC[HA.HC$pos%in%overlap,]
HA.HD<-HA.HD[HA.HD$pos%in%overlap,]
HA.HE<-HA.HE[HA.HE$pos%in%overlap,]
HC.HD<-HC.HD[HC.HD$pos%in%overlap,]
HC.HE<-HC.HE[HC.HE$pos%in%overlap,]
HD.HE<-HD.HE[HD.HE$pos%in%overlap,]
winpos<-HA.HC$pos

getPBS<-function(P1P2,P1P3,P1P4,P2P3,P2P4,P3P4){
	TP1P2=-log(1-P1P2)
	TP1P3=-log(1-P1P3)
	TP1P4=-log(1-P1P4)
	TP2P3=-log(1-P2P3)
	TP2P4=-log(1-P2P4)
	TP3P4=-log(1-P3P4)
	P1=((TP1P2+TP1P3-TP2P3)+(TP1P2+TP1P4-TP2P4)+(TP1P3+TP1P4-TP3P4))/6
	P2=((TP1P2+TP2P3-TP1P3)+(TP1P2+TP2P4-TP1P4)+(TP2P3+TP2P4-TP3P4))/6
	P3=((TP1P3+TP2P3-TP1P2)+(TP1P3+TP3P4-TP1P4)+(TP2P3+TP3P4-TP2P4))/6
	P4=((TP1P4+TP2P4-TP1P2)+(TP1P4+TP3P4-TP1P3)+(TP2P4+TP3P4-TP2P3))/6
	return(cbind(P1,P2,P3,P4))
}

getTrees<-function(vec){
	v<-(-log(1-vec))
	dmat<-matrix(nrow=4,ncol=4)
	colnames(dmat)<-c('HA','HC','HD','HE')
	rownames(dmat)<-c('HA','HC','HD','HE')
	dmat[lower.tri(dmat)]<-v
	currdist<-as.dist(dmat)
	nj.out<-bionj(currdist)
	return(nj.out)
}

Fst10kb<-data.frame(row.names=winpos,HA.HC=HA.HC$num/HA.HC$den)
Fst10kb$HA.HD<-HA.HD$num/HA.HD$den
Fst10kb$HA.HE<-HA.HE$num/HA.HE$den
Fst10kb$HC.HD<-HC.HD$num/HC.HD$den
Fst10kb$HC.HE<-HC.HE$num/HC.HE$den
Fst10kb$HD.HE<-HD.HE$num/HD.HE$den
Fst10kb<-na.omit(Fst10kb)
Fst10kb<-as.matrix(Fst10kb)
trees<-apply(Fst10kb,1,getTrees)

#####################
####################

#Figure S1
#Table S2

sns<-fread('snpfst_all_info.txt.gz')
colnames(sns)<-c('SNP','maf','HA.HC','HA.HD','HA.HE','HC.HD','HC.HE','HD.HE','Eff','gene','nuc','prot','warn')
sns$missense<-grepl('missense',sns$Eff)
sns$chr<-gsub(':.*','',sns$SNP)
sns$pos<-gsub('.*:','',sns$SNP)
sns$maf<-round(sns$maf,2)

table(sns[sns$HA.HE>0.99&sns$HC.HE>0.99&sns$HD.HE>0.99,'gene'])
table(sns[sns$HA.HE>0.99&sns$HC.HE>0.99&sns$HD.HE>0.99&sns$missense,'gene'])
dim(sns[sns$HA.HE>0.99&sns$HC.HE>0.99&sns$HD.HE>0.99&sns$missense,'gene'])
table(sns[sns$HA.HE>0.99&sns$HC.HE>0.99&sns$HD.HE>0.99,'chr'])

snsPBS<-getPBS(sns$HA.HC,sns$HA.HD,sns$HA.HE,sns$HC.HD,sns$HC.HE,sns$HD.HE)
sns$HA<-snsPBS[,1]
sns$HC<-snsPBS[,2]
sns$HD<-snsPBS[,3]
sns$HE<-snsPBS[,4]

s2<-sns[sns$missense&sns$HE>quantile(sns$HE,0.9),c(1:2,10:13,17:20)]
write.table(s2,sep='\t',quote=F,row.names=F,file='../Tables/TableS2.txt')

# shufs<-replicate(100,{
# shuffle(nrow(sns),control=how(blocks=sns$maf))
# })

# HAperms90<-apply(shufs,2,function(x) table(sns$HA[x][which(sns$missense)]>quantile(sns$HA,0.9))[2])
# HAobs90<-table(sns$HA[which(sns$missense)]>quantile(sns$HA,0.9))[2]
# HAperms95<-apply(shufs,2,function(x) table(sns$HA[x][which(sns$missense)]>quantile(sns$HA,0.95))[2])
# HAobs95<-table(sns$HA[which(sns$missense)]>quantile(sns$HA,0.95))[2]
# HAperms99<-apply(shufs,2,function(x) table(sns$HA[x][which(sns$missense)]>quantile(sns$HA,0.99))[2])
# HAobs99<-table(sns$HA[which(sns$missense)]>quantile(sns$HA,0.99))[2]

# HCperms90<-apply(shufs,2,function(x) table(sns$HC[x][which(sns$missense)]>quantile(sns$HC,0.9))[2])
# HCobs90<-table(sns$HC[which(sns$missense)]>quantile(sns$HC,0.9))[2]
# HCperms95<-apply(shufs,2,function(x) table(sns$HC[x][which(sns$missense)]>quantile(sns$HC,0.95))[2])
# HCobs95<-table(sns$HC[which(sns$missense)]>quantile(sns$HC,0.95))[2]
# HCperms99<-apply(shufs,2,function(x) table(sns$HC[x][which(sns$missense)]>quantile(sns$HC,0.99))[2])
# HCobs99<-table(sns$HC[which(sns$missense)]>quantile(sns$HC,0.99))[2]

# HDperms90<-apply(shufs,2,function(x) table(sns$HD[x][which(sns$missense)]>quantile(sns$HD,0.9))[2])
# HDobs90<-table(sns$HD[which(sns$missense)]>quantile(sns$HD,0.9))[2]
# HDperms95<-apply(shufs,2,function(x) table(sns$HD[x][which(sns$missense)]>quantile(sns$HD,0.95))[2])
# HDobs95<-table(sns$HD[which(sns$missense)]>quantile(sns$HD,0.95))[2]
# HDperms99<-apply(shufs,2,function(x) table(sns$HD[x][which(sns$missense)]>quantile(sns$HD,0.99))[2])
# HDobs99<-table(sns$HD[which(sns$missense)]>quantile(sns$HD,0.99))[2]

# HEperms90<-apply(shufs,2,function(x) table(sns$HE[x][which(sns$missense)]>quantile(sns$HE,0.9))[2])
# HEobs90<-table(sns$HE[which(sns$missense)]>quantile(sns$HE,0.9))[2]
# HEperms95<-apply(shufs,2,function(x) table(sns$HE[x][which(sns$missense)]>quantile(sns$HE,0.95))[2])
# HEobs95<-table(sns$HE[which(sns$missense)]>quantile(sns$HE,0.95))[2]
# HEperms99<-apply(shufs,2,function(x) table(sns$HE[x][which(sns$missense)]>quantile(sns$HE,0.99))[2])
# HEobs99<-table(sns$HE[which(sns$missense)]>quantile(sns$HE,0.99))[2]

# sns$loc1<-(sns$chr=='chr7'&sns$pos> 20350000&sns$pos<20570000)
# sns$loc2<-(sns$chr=='Sc0000015'&sns$pos> 2380000&sns$pos<2580000)
# sns2<-sns[sns$loc1==F & sns$loc2==F,]

# shufs2<-replicate(100,{
# shuffle(nrow(sns2),control=how(blocks=sns2$maf))
# })

# HEperms90_2<-apply(shufs2,2,function(x) table(sns2$HE[x][which(sns2$missense)]>quantile(sns2$HE,0.9))[2])
# HEobs90_2<-table(sns2$HE[which(sns2$missense)]>quantile(sns2$HE,0.9))[2]
# HEperms95_2<-apply(shufs2,2,function(x) table(sns2$HE[x][which(sns2$missense)]>quantile(sns2$HE,0.95))[2])
# HEobs95_2<-table(sns2$HE[which(sns2$missense)]>quantile(sns2$HE,0.95))[2]
# HEperms99_2<-apply(shufs2,2,function(x) table(sns2$HE[x][which(sns2$missense)]>quantile(sns2$HE,0.99))[2])
# HEobs99_2<-table(sns2$HE[which(sns2$missense)]>quantile(sns2$HE,0.99))[2]


# save(shufs,shufs2,HAperms90,HAperms95,HAperms99,HAobs90,HAobs95,HAobs99,
# HCperms90,HCperms95,HCperms99,HCobs90,HCobs95,HCobs99,
# HDperms90,HDperms95,HDperms99,HDobs90,HDobs95,HDobs99,
# HEperms90,HEperms95,HEperms99,HEobs90,HEobs95,HEobs99,
# HEperms90_2,HEperms95_2,HEperms99_2,HEobs90_2,HEobs95_2,HEobs99_2,
# file='sns_perms.Rdata')

getP<-function(perms,obs){
	currp<-(min(length(which(perms>obs)),length(which(perms<obs)))*2)/length(perms)
	if(currp==0) return(paste0('P<',1/length(perms)))
	return(paste0('P=',currp))
}

load('sns_perms.Rdata')

png('../Figures/FigureS2.png',width=8,height=6,res=300,units='in')
layout(rbind(c(1,2,2,2,2),c(3,4,5,6,7),c(8,9,10,11,12),c(13,14,15,16,17)))
par(mar=c(3,3,1,1),mgp=c(1.5,0.2,0),bty='n',tck=-0.01,oma=c(0,1,0,0),cex.main=1)

plot(density(sns$maf[sns$missense]),col='red',main='',xlab='MAF')
lines(density(sns$maf[!sns$missense]))
legend('top',fill=c('black','red'),legend=c('synonymous','missense'),bty='n',cex=0.7)
title('A',adj=0,font.main=1,cex.main=1.5)

snsPlot<-sns[sns$missense&sns$HE>0.1,]
snsPlot$col<-c(grey(0.2),grey(0.6))[1+as.numeric(factor(snsPlot$chr,levels=unique(snsPlot$chr)))%%2]
snsPlot$col[snsPlot$HE>quantile(sns$HE,0.9)]<-'tomato4'
snsPlot$col[snsPlot$HE>quantile(sns$HE,0.95)]<-'red4'
snsPlot$col[snsPlot$HE>quantile(sns$HE,0.99)]<-'red1'
plot(snsPlot$HE,cex=0.2,xlab='Missense SNPs',ylab='HE PBS',col=snsPlot$col,xlim=c(0,130000))
title('B',adj=0,font.main=1,cex.main=1.5)
legend('topright',fill=c('tomato4','red4','red1'),legend=c('top 10%','top 5%','top 1%'),bty='n',cex=0.7)

hist(HAperms90,col='grey',xlab='',ylab='Counts',xlim=c(50000,55000),main='HA',breaks=6)
abline(v=HAobs90,col='red')
title(main=getP(HAperms90,HAobs90),font.main=1,cex.main=0.8,adj=0.3,line=-1)
title('C',adj=0,font.main=1,cex.main=1.5)
mtext('top 10%',side=2,las=3,line=2.5)
hist(HCperms90,col='grey',xlab='',ylab='',xlim=c(50000,55000),main='HC',breaks=6)
abline(v=HCobs90,col='red')
title(main=getP(HCperms90,HCobs90),font.main=1,cex.main=0.8,adj=0.2,line=-1)
hist(HDperms90,col='grey',xlab='',ylab='',xlim=c(50000,55000),main='HD',breaks=6)
abline(v=HDobs90,col='red')
title(main=getP(HDperms90,HDobs90),font.main=1,cex.main=0.8,adj=0.08,line=-1)
hist(HEperms90,col='grey',xlab='',ylab='',xlim=c(50000,55000),main='HE',breaks=6)
abline(v=HEobs90,col='red')
title(main=getP(HEperms90,HEobs90),font.main=1,cex.main=0.8,adj=0.2,line=-1)
hist(HEperms90_2,col='grey',xlab='',ylab='Counts',xlim=c(50000,55000),main='HE (no L1,L2)',breaks=6)
abline(v=HEobs90_2,col='red')
title(main=getP(HEperms90_2,HEobs90_2),font.main=1,cex.main=0.8,adj=0.2,line=-1)

hist(HAperms95,col='grey',xlab='',ylab='Counts',xlim=c(24500,27500),main='',breaks=6)
abline(v=HAobs95,col='red')
title(main=getP(HAperms95,HAobs95),font.main=1,cex.main=0.8,adj=0.4,line=-1)
mtext('top 5%',side=2,las=3,line=2.5)
hist(HCperms95,col='grey',xlab='',ylab='',xlim=c(24500,27500),main='',breaks=6)
abline(v=HCobs95,col='red')
title(main=getP(HCperms95,HCobs95),font.main=1,cex.main=0.8,adj=0.2,line=-1)
hist(HDperms95,col='grey',xlab='',ylab='',xlim=c(24500,27500),main='',breaks=6)
abline(v=HDobs95,col='red')
title(main=getP(HDperms95,HDobs95),font.main=1,cex.main=0.8,adj=0.08,line=-1)
hist(HEperms95,col='grey',xlab='',ylab='',xlim=c(24500,27500),main='',breaks=6)
abline(v=HEobs95,col='red')
title(main=getP(HEperms95,HEobs95),font.main=1,cex.main=0.8,adj=0.1,line=-1)
hist(HEperms95_2,col='grey',xlab='',ylab='Counts',xlim=c(24500,27500),main='',breaks=6)
abline(v=HEobs95_2,col='red')
title(main=getP(HEperms95_2,HEobs95_2),font.main=1,cex.main=0.8,adj=0.1,line=-1)

hist(HAperms99,col='grey',xlab='Missense outliers',ylab='Counts',xlim=c(4500,6000),main='',breaks=6)
abline(v=HAobs99,col='red')
title(main=getP(HAperms99,HAobs99),font.main=1,cex.main=0.8,adj=0.1,line=-1)
mtext('top 1%',side=2,las=3,line=2.5)
hist(HCperms99,col='grey',xlab='Missense outliers',ylab='',xlim=c(4500,6000),main='',breaks=6)
abline(v=HCobs99,col='red')
title(main=getP(HCperms99,HCobs99),font.main=1,cex.main=0.8,adj=0.9,line=-1)
hist(HDperms99,col='grey',xlab='Missense outliers',ylab='',xlim=c(4500,6000),main='',breaks=6)
abline(v=HDobs99,col='red')
title(main=getP(HDperms99,HDobs99),font.main=1,cex.main=0.8,adj=0.08,line=-1)
hist(HEperms99,col='grey',xlab='Missense outliers',ylab='',xlim=c(4500,6000),main='',breaks=6)
abline(v=HEobs99,col='red')
title(main=getP(HEperms99,HEobs99),font.main=1,cex.main=0.8,adj=0.8,line=-1)
hist(HEperms99_2,col='grey',xlab='Missense outliers',ylab='Counts',xlim=c(4500,6000),main='',breaks=6)
abline(v=HEobs99_2,col='red')
title(main=getP(HEperms99_2,HEobs99_2),font.main=1,cex.main=0.8,adj=0.8,line=-1)
dev.off()

g1<-sns[sns$gene=='Amillepora12599',]
fisher.test(g1$HA.HE>0.99&g1$HC.HE>0.99&g1$HD.HE>0.99,g1$missense)
g2<-sns[sns$gene=='Amillepora12602',]
fisher.test(g2$HA.HE>0.99&g2$HC.HE>0.99&g2$HD.HE>0.99,g2$missense)
l1<-sns[sns$loc1,]
fisher.test(g2$HA.HE>0.99&g2$HC.HE>0.99&g2$HD.HE>0.99,g2$missense)

#########################
#########################

getTopo<-function(tree){
	tip<-which(tree$tip.label=='HE')
	tree$tip.label[tree$edge[(tree$edge[,1]==tree$edge[tree$edge[,2]==tip,1]) & tree$edge[,2]<5 & tree$edge[,2]!=tip,2]]
}

drawTree<-function(x,y,lab,h,v,t){
	segments(x-h,y,x+h,y)
	segments(x-h,y,x-(2*h),y+v)
	segments(x-h,y,x-(2*h),y-v)
	segments(x+h,y,x+(2*h),y+v)
	segments(x+h,y,x+(2*h),y-v)
	text(c(x-(2*h)-t,x-(2*h)-t,x+(2*h)+t,x+(2*h)+t),c(y-v,y+v,y-v,y+v),labels=lab)
}

topos<-sapply(trees,getTopo)
png('../Figures/FigureS1.png',width=4,height=4,res=300,units='in')
par(mar=c(3,4,2,1),mgp=c(2.5,0.1,0),tck=-0.01)
barplot(as.matrix(table(topos)),ylab='Counts',xlim=c(0,3))
drawTree(2,43000,c('HA','HC','HD','HE'),0.1,2000,0.2)
drawTree(2,35000,c('HA','HD','HC','HE'),0.1,2000,0.2)
drawTree(2,18000,c('HD','HC','HA','HE'),0.1,2000,0.2)
dev.off()
###################
####################

############
#Figure 1
#############

png('../Figures/Figure1.png',width=6,height=6,res=300,units='in')
layout(rbind(c(1,1,1),c(2,3,4),c(5,6,7)))
par(mar=c(3,3,2,1),mgp=c(1.7,0.1,0),las=1,tck=0.01,cex.axis=1,lwd=1,cex=0.5,bty='l',cex.main=1.5)

crypcol= rainbow(100,v=0.8)[c(18,5,50,80)]
crypcols= crypcol[factor(meta$cryp)]

#A
bp<-barplot(admix,col=crypcol,names.arg=rep('',ncol(admix)),las=1,ylab='Proportion ancestry')
mtext('HA',side=1,line=1,at=bp[7])
mtext('HC',side=1,line=1,at=bp[25])
mtext('HD',side=1,line=1,at=bp[50])
mtext('HE',side=1,line=1,at=bp[90])
title('A',adj=0)

#B

covmat<-as.matrix(read.delim('pcangsd.cov',header=F))
colnames(covmat)<-samps
rownames(covmat)<-samps
covmat<-covmat[match(meta$filename,rownames(covmat)),match(meta$filename,colnames(covmat))]
e<-eigen(covmat)
plot(e$vectors[,1],e$vectors[,2],col=crypcols,xlab='PC 1',ylab='PC 2')
title('B',adj=0)

#C

globalfst<-read.delim('global_fst.txt',header=F,sep=' ')
globaltree<-getTrees(globalfst[,2])
globaltree$tip.label<-c('HA','HC','HD','HE')
plot(globaltree,type='u',cex=1, edge.width=2,lab4ut='axial') 
segments(.05,0,.1,0,lwd=2)
text(0.075,0.01,label=expression(italic(T) == 0.05))
title('C',adj=0)

#D
plotperm<-data.frame(perms=c(HAperms90,HCperms90,HDperms90,HEperms90),ind=c(rep(1,100),rep(2,100),rep(3,100),rep(4,100)))
plot(jitter(plotperm$ind,0.3),plotperm$perms,cex=1,ylim=c(50000,55500),lwd=0.5,col='grey',xlab='Species',ylab=expression('Missense outliers x'~10^3),axes=F,xlim=c(0.4,4.5))
axis(1,at=1:4,labels=c('HA','HC','HD','HE'))
axis(2,at=c(50000,52000,54000),labels=c(50,52,54))
points(1:4,c(HAobs90,HCobs90,HDobs90,HEobs90),cex=1,col='red',pch=19)
box(bty='l')
legend('topright',fill=c('grey','red'),legend=c('permuted','observed'),bty='n')
title('D',adj=0)

#E
seq<-ape::read.dna('ahy_mitosnps.fa',format='fasta')
meta2<-meta[match(labels(seq),meta$filename),]
table(labels(seq)==meta2$filename)
d<-dist.dna(seq)
hc<-(hclust(d))
cut<-cutree(hc,k=2)
meta$mitotype<-cut[match(meta$filename,meta2$filename)]

h<-pegas::haplotype(seq)
h<-sort(h,what='label')
labs=labels(h)
l<-rep(NA,nrow(seq))
for(i in 1:length(labs)){
	l[attr(h,'index')[[i]]]<-labs[i]
}
sp=table(l,meta2$cryp)
net<-pegas::haploNet(h)
table(attr(net,'labels')==rownames(sp))
attr(net,'alter.links')<-NULL

plot(net,size=attr(net,'freq')**.5,pie=sp,show.mutation=2,scale.ratio=1,fast=F,labels=F,bg=crypcol)
legend('topleft',legend=c('HA','HC','HD','HE'),fill=crypcol,bty='n')
title('E',adj=0)
par(xpd=F)
text(0,-10,label='HG1')
text(30,-10,label='HG2')


#F
plot(jitter(as.numeric(factor(meta$cryp))),meta$propd,axes=F,xlab='Species',ylab=expression(paste("Proportion",italic(" Durusdinium")," symbionts")),type='n',cex=1,ylim=c(0,1),las=1,cex.axis=1)
axis(1,at=1:4,labels=c('HA','HC','HD','HE'))
axis(2)
box(bty='l')
abline(h=0.5,lwd=2,lty=2)
points(jitter(as.numeric(factor(meta$cryp))),meta$propd,col=crypcols,cex=1)

title('F',adj=0)


#G
bl<-table(meta$binbl,interaction(meta$cryp,meta$propd>0.5))
bltab<-matrix(0,5,8)
bltab[5,]<-bl[2,]
bltab[1,c(1,5)]<-bl[1,c(1,5)]
bltab[2,c(2,6)]<-bl[1,c(2,6)]
bltab[3,c(3,7)]<-bl[1,c(3,7)]
bltab[4,c(4,8)]<-bl[1,c(4,8)]
bp<-barplot(bltab,col=c(crypcol,grey(0.95)),names.arg=c('HA','HC','HD','HE','HA','HC','HD','HE'),cex.names=1,ylab='Count',las=2)
mtext('Sym. C',at=mean(bp[2:3]),side=1,las=1,line=1.5,cex=0.7)
mtext('Sym. D',at=mean(bp[6:7]),side=1,las=1,line=1.5,cex=0.7)
legend('topleft',fill=grey(0.95),legend='bleached',bty='n')
title('G',adj=0)


dev.off()


#############

goodcov<-read.delim('good_coverage.txt',na.strings='.',header=F)
rownames(goodcov)<-paste(goodcov[,1],goodcov[,2],sep=':')
goodcov<-goodcov[winpos,]
repeatcov<-read.delim('repeat_coverage.bed',na.strings='.',header=F)
rownames(repeatcov)<-paste(repeatcov[,1],repeatcov[,2],sep=':')
repeatcov<-repeatcov[winpos,]
cor.test(repeatcov[,7],goodcov[,7])

piHA<-read.delim('HA.pestPG')
piHA$winpos<-paste(piHA$Chr,piHA$WinCenter-5000,sep=':')
piHA<-piHA[match(winpos,piHA$winpos),]
piHC<-read.delim('HC.pestPG')
piHC$winpos<-paste(piHC$Chr,piHC$WinCenter-5000,sep=':')
piHC<-piHC[match(winpos,piHC$winpos),]
piHD<-read.delim('HD.pestPG')
piHD$winpos<-paste(piHD$Chr,piHD$WinCenter-5000,sep=':')
piHD<-piHD[match(winpos,piHD$winpos),]
piHE<-read.delim('HE.pestPG')
piHE$winpos<-paste(piHE$Chr,piHE$WinCenter-5000,sep=':')
piHE<-piHE[match(winpos,piHE$winpos),]

meanPi<-(piHA$tP/goodcov[,4]+piHC$tP/goodcov[,4]+piHD$tP/goodcov[,4]+piHE$tP/goodcov[,4])/4

edges<-c()
for(chr in unique(HA.HE$chr)){
	curr<-which(HA.HE$chr==chr)
	if(length(curr)<10){
		edges<-c(edges,curr)
		next
	}
	edges<-c(edges,curr[1:5],curr[(length(curr)-4):length(curr)])
}

getPBSwin<-function(focal.ref,focal.out,ref.out,k=10,edges){
	T1=-log(1-rollapply(focal.ref$num,FUN=mean,width=k,na.rm=T,fill=NA)/rollapply(focal.ref$den,FUN=mean,width=k,na.rm=T, fill=NA))
	T2=-log(1-rollapply(focal.out $num,FUN=mean,width=k,na.rm=T,fill=NA)/rollapply(focal.out $den,FUN=mean,width=k,na.rm=T, fill=NA))
	T3=-log(1-rollapply(ref.out $num,FUN=mean,width=k,na.rm=T,fill=NA)/rollapply(ref.out $den,FUN=mean,width=k,na.rm=T, fill=NA))
	PBS=(T1+T2-T3)/2
	PBS[edges]<-NA
	return(PBS)
}


# P1P2=HA.HC
# P1P3=HA.HD
# P1P4=HA.HE
# P2P3=HC.HD
# P2P4=HC.HE
# P3P4=HD.HE

getMeanPBS<-function(P1P2,P1P3,P1P4,P2P3,P2P4,P3P4,k=10,edges,remove.edges=T){
	TP1P2=-log(1-rollapply(P1P2$num,FUN=mean,width=k,na.rm=T,fill=NA)/rollapply(P1P2$den,FUN=mean,width=k,na.rm=T, fill=NA))
	TP1P3=-log(1-rollapply(P1P3$num,FUN=mean,width=k,na.rm=T,fill=NA)/rollapply(P1P3$den,FUN=mean,width=k,na.rm=T, fill=NA))
	TP1P4=-log(1-rollapply(P1P4$num,FUN=mean,width=k,na.rm=T,fill=NA)/rollapply(P1P4$den,FUN=mean,width=k,na.rm=T, fill=NA))
	TP2P3=-log(1-rollapply(P2P3$num,FUN=mean,width=k,na.rm=T,fill=NA)/rollapply(P2P3$den,FUN=mean,width=k,na.rm=T, fill=NA))
	TP2P4=-log(1-rollapply(P2P4$num,FUN=mean,width=k,na.rm=T,fill=NA)/rollapply(P2P4$den,FUN=mean,width=k,na.rm=T, fill=NA))
	TP3P4=-log(1-rollapply(P3P4$num,FUN=mean,width=k,na.rm=T,fill=NA)/rollapply(P3P4$den,FUN=mean,width=k,na.rm=T, fill=NA))
	P1=((TP1P2+TP1P3-TP2P3)+(TP1P2+TP1P4-TP2P4)+(TP1P3+TP1P4-TP3P4))/6
	P2=((TP1P2+TP2P3-TP1P3)+(TP1P2+TP2P4-TP1P4)+(TP2P3+TP2P4-TP3P4))/6
	P3=((TP1P3+TP2P3-TP1P2)+(TP1P3+TP3P4-TP1P4)+(TP2P3+TP3P4-TP2P4))/6
	P4=((TP1P4+TP2P4-TP1P2)+(TP1P4+TP3P4-TP1P3)+(TP2P4+TP3P4-TP2P3))/6
	res<-as.data.frame(cbind(P1,P2,P3,P4))
	if(remove.edges) res[edges,]<-NA
	return(res)
}
meanPBS<-getMeanPBS(HA.HC,HA.HD,HA.HE,HC.HD,HC.HE,HD.HE,edges=edges)

getPerm<-function(P1P2,P1P3,P1P4,P2P3,P2P4,P3P4,pis,edges,reps=10,k=10){
	perms<-replicate(reps,{
		shuf<-shuffle(length(pis),control=how(blocks=cut(pis,quantile(pis,seq(0,1,0.1),na.rm=T))))
		getMeanPBS(P1P2[shuf,],P1P3[shuf,],P1P4[shuf,],P2P3[shuf,],P2P4[shuf,],P3P4[shuf,],k=k,edges=edges)
	})
	return(perms)
}

PBS<-data.frame(row.names=paste(HA.HC$chr,HA.HC$start,sep=':'),chr=HA.HC$chr)
PBS$chr<-factor(PBS$chr,levels=c(paste('chr',1:14,sep=''),'unplaced'))
PBS$chr[is.na(PBS$chr)]<-'unplaced'
PBS$HA<-meanPBS[,1]
PBS$HC<-meanPBS[,2]
PBS$HD<-meanPBS[,3]
PBS$HE<-meanPBS[,4]
PBS$chrcol<-c(grey(0.2),grey(0.4))[1+as.numeric(PBS$chr)%%2]
PBS$chrcol[PBS$chr=='unplaced']<-grey(0.7)
PBS$chrcolsave<-PBS$chrcol

# perms<-getPerm(HA.HC,HA.HD,HA.HE,HC.HD,HC.HE,HD.HE, pis=meanPi,edges=edges,reps=10,k=10)
# permHA<-perms[,1,]
# permHC<-perms[,2,]
# permHD<-perms[,3,]
# permHE<-perms[,4,]
# ecdf.HA<-ecdf(permHA)
# ecdf.HC<-ecdf(permHC)
# ecdf.HD<-ecdf(permHD)
# ecdf.HE<-ecdf(permHE)
# save(perms,permHA,permHC,permHD,permHE,ecdf.HA,ecdf.HC,ecdf.HD,ecdf.HE,file='perms_meanPBS.Rdata')

load('perms_meanPBS.Rdata')

outlpadjHA<-p.adjust(2*pmin((1-ecdf.HA(PBS$HA)),ecdf.HA(PBS$HA)),method='BH')
outlpadjHC<-p.adjust(2*pmin((1-ecdf.HC(PBS$HC)),ecdf.HC(PBS$HC)),method='BH')
outlpadjHD<-p.adjust(2*pmin((1-ecdf.HD(PBS$HD)),ecdf.HD(PBS$HD)),method='BH')
outlpadjHE<-p.adjust(2*pmin((1-ecdf.HE(PBS$HE)),ecdf.HE(PBS$HE)),method='BH')

PBS$colchrom<-rainbow(11,v=0.8)[1+as.numeric(factor(HA.HC$chr))%%11]
overplot<-which(outlpadjHE<0.05&PBS$HE>mean(PBS$HE,na.rm=T))
PBS$chrcol[overplot]<-PBS$colchrom[overplot]
PBS$chrnum<-(as.numeric(gsub('chr','',PBS$chr)))
PBS$chrnum[is.na(PBS$chrnum)]<-15
PBS$start<-HA.HC$start

PBS$HAcol=PBS$chrcolsave
HAover<-which(outlpadjHA<0.05&PBS$HA>mean(PBS$HA,na.rm=T))
PBS$HAcol[HAover]<-PBS$colchrom[HAover]
PBS$HCcol=PBS$chrcolsave
HCover<-which(outlpadjHC<0.05&PBS$HC>mean(PBS$HC,na.rm=T))
PBS$HCcol[HCover]<-PBS$colchrom[HCover]
PBS$HDcol=PBS$chrcolsave
HDover<-which(outlpadjHD<0.05&PBS$HD>mean(PBS$HD,na.rm=T))
PBS$HDcol[HDover]<-PBS$colchrom[HDover]
PBS$HAoutl<-outlpadjHA<0.05&PBS$HA>mean(PBS$HA,na.rm=T)
PBS$HCoutl<-outlpadjHC<0.05&PBS$HC>mean(PBS$HC,na.rm=T)
PBS$HDoutl<-outlpadjHD<0.05&PBS$HD>mean(PBS$HD,na.rm=T)
PBS$HEoutl<-outlpadjHE<0.05&PBS$HE>mean(PBS$HE,na.rm=T)

HAwin<-t(sapply(winpos[HAover],function(v) c(strsplit(v,split=':')[[1]])))
HAbed<-cbind(HAwin[,1],as.numeric(HAwin[,2])-50000,as.numeric(HAwin[,2])+50000)
HCwin<-t(sapply(winpos[HCover],function(v) c(strsplit(v,split=':')[[1]])))
HCbed<-cbind(HCwin[,1],as.numeric(HCwin[,2])-50000,as.numeric(HCwin[,2])+50000)
HDwin<-t(sapply(winpos[HDover],function(v) c(strsplit(v,split=':')[[1]])))
HDbed<-cbind(HDwin[,1],as.numeric(HDwin[,2])-50000,as.numeric(HDwin[,2])+50000)
HEwin<-t(sapply(winpos[overplot],function(v) c(strsplit(v,split=':')[[1]])))
HEbed<-cbind(HEwin[,1],as.numeric(HEwin[,2])-50000,as.numeric(HEwin[,2])+50000)
write.table(HAbed,row.names=F,col.names=F,quote=F,file='HA_outlier_regions.bed',sep='\t')
write.table(HCbed,row.names=F,col.names=F,quote=F,file='HC_outlier_regions.bed',sep='\t')
write.table(HDbed,row.names=F,col.names=F,quote=F,file='HD_outlier_regions.bed',sep='\t')
write.table(HEbed,row.names=F,col.names=F,quote=F,file='HE_outlier_regions.bed',sep='\t')

plotPBS<-na.omit(PBS[order(PBS$chrnum),])
plotPBS$HA[plotPBS$HA<0]<-0
plotPBS$HC[plotPBS$HC<0]<-0
plotPBS$HD[plotPBS$HD<0]<-0
plotPBS$HE[plotPBS$HE<0]<-0

HA.HE.dxy<-read.delim('HA.HE.dxy',header=F,na.strings='.')
HA.HE.dxy[is.na(HA.HE.dxy[,4]),4]<-0
rownames(HA.HE.dxy)<-paste(HA.HE.dxy[,1],HA.HE.dxy[,2],sep=':')
HA.HE.dxy<-HA.HE.dxy[winpos,]

HC.HE.dxy<-read.delim('HC.HE.dxy',header=F,na.strings='.')
HC.HE.dxy[is.na(HC.HE.dxy[,4]),4]<-0
rownames(HC.HE.dxy)<-paste(HC.HE.dxy[,1],HC.HE.dxy[,2],sep=':')
HC.HE.dxy<-HC.HE.dxy[winpos,]

HD.HE.dxy<-read.delim('HD.HE.dxy',header=F,na.strings='.')
HD.HE.dxy[is.na(HD.HE.dxy[,4]),4]<-0
rownames(HD.HE.dxy)<-paste(HD.HE.dxy[,1],HD.HE.dxy[,2],sep=':')
HD.HE.dxy<-HD.HE.dxy[winpos,]


####################
#Figure 2          
####################

HA.HE.rolldxy<-rollapply(HA.HE.dxy[,4],FUN=sum,na.rm=T,width=10,fill=NA)/rollapply(goodcov[,4],FUN=sum,na.rm=T,width=10,fill=NA)
HA.rollpi<-rollapply(piHA$tP,width=10,FUN=sum,na.rm=T,fill=NA)/rollapply(goodcov[,4],FUN=sum,na.rm=T,width=10,fill=NA)
HE.rollpi<-rollapply(piHE$tP,width=10,FUN=sum,na.rm=T,fill=NA)/rollapply(goodcov[,4],FUN=sum,na.rm=T,width=10,fill=NA)
nsites<-rollapply(goodcov[,4],FUN=sum,na.rm=T,width=10,fill=NA)
HA.HE.rollpi<-(HA.rollpi+HE.rollpi)/2

HC.HE.rolldxy<-rollapply(HC.HE.dxy[,4],FUN=sum,na.rm=T,width=10,fill=NA)/rollapply(goodcov[,4],FUN=sum,na.rm=T,width=10,fill=NA)
HC.rollpi<-rollapply(piHC$tP,width=10,FUN=sum,na.rm=T,fill=NA)/rollapply(goodcov[,4],FUN=sum,na.rm=T,width=10,fill=NA)
HE.rollpi<-rollapply(piHE$tP,width=10,FUN=sum,na.rm=T,fill=NA)/rollapply(goodcov[,4],FUN=sum,na.rm=T,width=10,fill=NA)
HC.HE.rollpi<-(HC.rollpi+HE.rollpi)/2

HD.HE.rolldxy<-rollapply(HD.HE.dxy[,4],FUN=sum,na.rm=T,width=10,fill=NA)/rollapply(goodcov[,4],FUN=sum,na.rm=T,width=10,fill=NA)
HD.rollpi<-rollapply(piHD$tP,width=10,FUN=sum,na.rm=T,fill=NA)/rollapply(goodcov[,4],FUN=sum,na.rm=T,width=10,fill=NA)
HE.rollpi<-rollapply(piHE$tP,width=10,FUN=sum,na.rm=T,fill=NA)/rollapply(goodcov[,4],FUN=sum,na.rm=T,width=10,fill=NA)
HD.HE.rollpi<-(HD.rollpi+HE.rollpi)/2

png('../Figures/Figure2.png',width=6,height=4,res=300,units='in')
par(mar=c(3,3,2,1),mgp=c(1.7,0.1,0),las=1,tck=0.01,cex.axis=1,lwd=1,bty='l')
layout(rbind(c(1,1,1),c(2,3,4)))
plot(plotPBS$HE,col=plotPBS$chrcol,ylab='HE PBS',xlab='Genomic position (100 kb window, 10 kb step)',axes=F,cex=0.5,ylim=c(0,1.8))
axis(2)
chrlab<-aggregate((1:nrow(plotPBS))~plotPBS$chr,FUN=mean)
axis(1,labels=c(1:14,'unplaced'),at=chrlab[,2],lwd=0)
title('A',adj=0)
text(19600,1.8,labels='HES1')
text(34100,1.18,labels='HES2')

plot(HA.HE.rolldxy[nsites>50000],HA.HE.rollpi[nsites>50000],col=PBS$chrcol[nsites>50000],cex=0.5,xlab=expression(HA~vs.~HE~pi[b]),ylab=expression(pi[w]))
points(HA.HE.rolldxy[overplot],HA.HE.rollpi[overplot],col=PBS$chrcol[overplot],cex=0.5)
title('B',adj=0)

plot(HC.HE.rolldxy[nsites>50000],HC.HE.rollpi[nsites>50000],col=PBS$chrcol[nsites>50000],cex=0.5,xlab=expression(HC~vs.~HE~pi[b]),ylab=expression(pi[w]))
points(HC.HE.rolldxy[overplot],HC.HE.rollpi[overplot],col=PBS$chrcol[overplot],cex=0.5)
title('C',adj=0)

plot(HD.HE.rolldxy[nsites>50000],HD.HE.rollpi[nsites>50000],col=PBS$chrcol[nsites>50000],cex=0.5,xlab=expression(HD~vs.~HE~pi[b]),ylab=expression(pi[w]))
points(HD.HE.rolldxy[overplot],HD.HE.rollpi[overplot],col=PBS$chrcol[overplot],cex=0.5)
title('D',adj=0)

dev.off()


################################
#Figure S1
###############################

png('../Figures/FigureS3.png',width=6,height=6,res=300,units='in')
par(mar=c(3,3,2,1),mgp=c(1.7,0.1,0),las=1,tck=0.01,cex.axis=1,lwd=1,bty='l',mfrow=c(4,1))
plot(plotPBS$HA,col=plotPBS$HAcol,ylab='HA PBS',xlab='',axes=F,cex=0.5,ylim=c(0,1.7))
points(which(plotPBS$HAoutl),plotPBS$HA[plotPBS$HAoutl],col=plotPBS$HAcol[plotPBS$HAoutl],cex=0.5)
axis(2)
chrlab<-aggregate((1:nrow(plotPBS))~plotPBS$chr,FUN=mean)
axis(1,labels=c(1:14,'unplaced'),at=chrlab[,2],lwd=0)
title('A',adj=0)

plot(plotPBS$HC,col=plotPBS$HCcol,ylab='HC PBS',xlab='',axes=F,cex=0.5,ylim=c(0,1.7))
points(which(plotPBS$HCoutl),plotPBS$HC[plotPBS$HCoutl],col=plotPBS$HCcol[plotPBS$HCoutl],cex=0.5)
axis(2)
chrlab<-aggregate((1:nrow(plotPBS))~plotPBS$chr,FUN=mean)
axis(1,labels=c(1:14,'unplaced'),at=chrlab[,2],lwd=0)
title('B',adj=0)

plot(plotPBS$HD,col=plotPBS$HDcol,ylab='HD PBS',xlab='',axes=F,cex=0.5,ylim=c(0,1.7))
points(which(plotPBS$HDoutl),plotPBS$HD[plotPBS$HDoutl],col=plotPBS$HDcol[plotPBS$HDoutl],cex=0.5)
axis(2)
chrlab<-aggregate((1:nrow(plotPBS))~plotPBS$chr,FUN=mean)
axis(1,labels=c(1:14,'unplaced'),at=chrlab[,2],lwd=0)
title('C',adj=0)

plot(plotPBS$HE,col=plotPBS$chrcol,ylab='HE PBS',xlab='Genomic position  (100 kb window, 10 kb step)',axes=F,cex=0.5,ylim=c(0,1.7))
points(which(plotPBS$HEoutl),plotPBS$HE[plotPBS$HEoutl],col=plotPBS$chrcol[plotPBS$HEoutl],cex=0.5)
axis(2)
chrlab<-aggregate((1:nrow(plotPBS))~plotPBS$chr,FUN=mean)
axis(1,labels=c(1:14,'unplaced'),at=chrlab[,2],lwd=0)
title('D',adj=0)

dev.off()

#################

####################
# Figure 3
####################

covmat<-as.matrix(read.delim('chr7:20350000-20570000.cov',header=F))
colnames(covmat)<-samps
rownames(covmat)<-samps
covmat<-covmat[match(meta$filename,rownames(covmat)),match(meta$filename,colnames(covmat))]
e1<-eigen(covmat)
e1$values[1]/114
# plot(e1$vectors[,1],e1$vectors[,2],col=crypcols,xlab='PC1 (59% of variation)',ylab='PC2 (7% of variation)')

covmat<-as.matrix(read.delim('Sc0000015:2400000-2570000.cov',header=F))
colnames(covmat)<-samps
rownames(covmat)<-samps
covmat<-covmat[match(meta$filename,rownames(covmat)),match(meta$filename,colnames(covmat))]
e2<-eigen(covmat)
e2$values[1]/114
# plot(e2$vectors[,1],e2$vectors[,2],col=crypcols,xlab='PC1 (69% of variation)',ylab='PC2 (7% of variation)')


winpos[overplot]
#chr7:20350000-20570000
#Sc0000015:2400000-2570000

ahy2amil<-read.delim('ahy2amil_blast.txt',header=F)
amil2ahy<-read.delim('amil2ahy_blast.txt',header=F)
ahy2amil$recip<-amil2ahy[match(ahy2amil[,2],amil2ahy[,1]),2]
ahy2amil$rbh<-ahy2amil[,1]==ahy2amil$recip
ahy2amil[ahy2amil[,2]=='Amillepora12599-RA','recip']
max(c(ahy2amil[,11]),c(amil2ahy[,11]))
de<-read.delim('HEvsHC_DE_contigs.txt',header=F)
ahy2amil$de<-ahy2amil[,1]%in%de[,1]
annot<-read.delim('~/Downloads/Amil_v2.01/Amillepora_euk.emapper.annotations',comment.char='#',header=F)
genes<-read.delim('outlier_genes.bed',header=F)
genes$gene<-gsub(';Name=.*','',gsub("ID=",'',genes[,9]))
genes$de<-genes$gene%in%gsub('-RA','',ahy2amil[ahy2amil$de,2])
genes$gene2<-gsub('Amillepora','gene ',gsub(';Name=.*','',gsub("ID=",'',genes[,9])))
genes$annot<-annot[match(genes$gene,annot[,1]),5]
genes$annot[is.na(genes$annot)]<-''
genes$description<-annot[match(genes$gene,annot[,1]),13]
genes[genes$gene=='Amillepora27759','annot']<-'NDUFA10'
genes[genes$gene=='Amillepora27763','annot']<-'MULE'
genes$gene2<-paste(genes$gene2,genes$annot)
win10<-getMeanPBS(HA.HC,HA.HD,HA.HE,HC.HD,HC.HE,HD.HE,edges=edges,k=1,remove.edges=F)
PBS$HEwin<-win10[,4]


########

png('../Figures/Figure3.png',width=6,height=6,res=300,units='in')
par(mar=c(3,3,2,1),mgp=c(1.7,0.1,0),las=1,tck=0.01,cex.axis=1,lwd=1,bty='l',mfrow=c(2,2))

L1<-PBS[PBS$chr=='chr7'&PBS$start>=20250000&PBS$start<20670000,]
L1sns<-sns[sns$chr=='chr7'&sns$pos>= 20250000&sns$pos<20670000&sns$missense,]
plot(L1$start,L1$HEwin,type='l',ylab='HE PBS',xlab='Chromosome 7 position (Mb)',ylim=c(0,8.5),xlim=c(20250000, 20670000),axes=F)
points(L1sns$pos,L1sns$HE,cex=0.2,col='darkred')
axis(2,at=c(0,2,4,6))
axis(1,at=seq(20300000,20600000,by=100000),labels=seq(20300000,20600000,by=100000)/1000000)
g7<-genes[genes[,1]=='chr7',]
g7<-g7[order(g7[,4]),]
g7$height<-rep(c(6,6.3,6.6,6.9,7.2,7.5,7.8,8.1),10)[1:nrow(g7)]
segments(g7[,4],g7$height,g7[,5],g7$height)
text((g7[,4]+g7[,5])/2,g7$height+0.13,labels=g7$gene2,cex=0.5,font=1+as.numeric(g7$de))
title('A',adj=0)
legend('topleft',legend='missense SNPs',fill='darkred',bty='n',cex=0.7)

L2<-PBS[HA.HC$chr=='Sc0000015'&PBS$start>=2350000&PBS$start<2650000,]
L2sns<-sns[sns$chr=='Sc0000015'&sns$pos>=2350000&sns$pos<2650000&sns$missense,]
plot(L2$start,L2$HEwin,type='l',ylab='HE PBS',xlab='Sc0000015 position (Mb)',ylim=c(0,4.5),xlim=c(2300000, 2650000),axes=F)
points(L2sns$pos,L2sns$HE,cex=0.2,col='darkred')
axis(2,at=c(0:3))
axis(1,at=seq(2340000, 2650000,by=100000),labels=seq(2340000, 2650000,by=100000)/1000000)
gs<-genes[genes[,1]=='Sc0000015',]
gs<-gs[order(gs[,4]),]
gs$height<-rep(seq(1.1,1.6,.1),10)[1:nrow(gs)]
gs$height<-rep(c(2.4,2.6,2.8,3.2,3.4,3.6,3.8,4,4.2),10)[1:nrow(gs)]
segments(gs[,4],gs$height,gs[,5],gs$height)
text((gs[,4]+gs[,5])/2,gs$height+0.1,labels=gs$gene2,cex=0.5,font=1+as.numeric(gs$de))
title('B',adj=0)
legend('topleft',legend='missense SNPs',fill='darkred',bty='n',cex=0.7)

bleachcol<-crypcols
bleachcol[meta$binbl==1]<-grey(0.9)
plot(jitter(as.numeric(meta$cryp)),e1$vectors[,1],bg=bleachcol,pch=21,axes=F,ylab='HES1 PC1 (62% of variation)',xlab='Cryptic species')
axis(1,at=1:4,labels=c('HA','HC','HD','HE'),lwd=0)
axis(2,at=c(-0.1,0,0.1))
title('C',adj=0)
legend('top',legend='bleached',fill=grey(0.95),bty='n')

plot(jitter(as.numeric(meta$cryp)),e2$vectors[,1],bg=bleachcol,pch=21,axes=F,ylab='HES2 PC1 (70% of variation)',xlab='Cryptic species')
axis(1,at=1:4,labels=c('HA','HC','HD','HE'),lwd=0)
axis(2,at=c(-0.1,0,0.1))
title('D',adj=0)

legend('top',legend='bleached',fill=grey(0.95),bty='n')
arrows(1.6,0.03,1.3,0.05,length=0.05,lwd=1)
arrows(1.6,0.02,1.3,0,length=0.05,lwd=1)
arrows(2.6,-0.03,2.3,-0.03,length=0.05,lwd=1)
arrows(3.6,-0.05,3.3,-0.05,length=0.05,lwd=1)

dev.off()


######
cryp1<-ape::read.dna(gzfile('chr1.fna.gz'),format='fasta')
cryp2<-ape::read.dna(gzfile('chr2.fna.gz'),format='fasta')

png('../Figures/Figure4.png',width=5,height=6,res=300,units='in')
par(mar=c(3,3,2,1),mgp=c(1.7,0.1,0),las=1,tck=0.01,cex.axis=1,lwd=1,bty='l',mfrow=c(3,2),font.main=1)

Amillepora12599<-read.tree('proteins/Amillepora12599.new')
Amillepora12599$tip.label[5]<-'A. mil.'
plot(Amillepora12599,cex=1,main='gene 12599',font.main=1)
title('A',adj=0)
segments(0,1.2,0.01,1.2,lwd=1)
text(0.005,1.4,labels='1%')

Amillepora12602<-read.tree('proteins/Amillepora12602.new')
Amillepora12602 $tip.label[5]<-'A. mil.'
plot(Amillepora12602,cex=1,main='gene 12602',font.main=1)
title('B',adj=0)
segments(0,1.2,0.01,1.2,lwd=1)
text(0.005,1.4,labels='1%')

locus1cryp<-ape::read.dna('Locus1.fna',format='fasta')
locus1dist<-dist.dna(locus1cryp,model='raw')
locus1.nj<-bionj(locus1dist)
locus1.nj$tip.label[1]<-'A. mil.'
plot(locus1.nj,type='u',main='HES1')
segments(0,0.02,0.002,0.02,lwd=1)
text(0.0025,0.023,labels='0.5%')
title('C',adj=0)

locus2cryp<-ape::read.dna('Locus2.fna',format='fasta')
locus2dist<-dist.dna(locus2cryp,model='raw')
locus2.nj<-bionj(locus2dist)
locus2.nj$tip.label[1]<-'A. mil.'
plot(locus2.nj,type='u',main='HES2')
segments(0.005,0.002,0.005,0.007,lwd=1)
text(0.0065,0.0045,labels='0.5%')
title('D',adj=0)

crypdist<-dist.dna(cryp1,model='raw')
cryp.nj<-bionj(crypdist)
cryp.nj$tip.label[1]<-'A. mil.'
plot(cryp.nj,type='u',main='Chromosome 1')
segments(0,0.003,0.005,0.003,lwd=1)
text(0.0025,0.004,labels='0.5%')
title('E',adj=0)

crypdist<-dist.dna(cryp2,model='raw')
cryp.nj<-bionj(crypdist)
cryp.nj$tip.label[1]<-'A. mil.'
plot(cryp.nj,type='u',main='Chromosome 2')
segments(0,0.01,0.005,0.01,lwd=1)
text(0.0025,0.011,labels='0.5%')
title('F',adj=0)

dev.off()


################
# Figure S2
##################

tmf<-read.delim('tmf_counts.txt',header=F)
tmf$loc<-e2$vectors[match(tmf[,1],as.character(meta$tag)),1]
tmf$cryp<-meta$cryp[match(tmf[,1],as.character(meta$tag))]
cor(tmf$loc,tmf[,2],use='p')
anova(lm(tmf[,2]~tmf$cryp+tmf$loc))

png('../Figures/FigureS4.png',width=3,height=3,res=300,units='in')
par(mar=c(3,2,1,1),tck=-0.01,mgp=c(1.1,.1,0),bty='l')
plot(tmf$loc,tmf[,2],bg='black',col=crypcols[match(tmf[,1],as.character(meta$tag))],pch=19,xlab='PC1 HES2',ylab='TMF norm. counts')
abline(lm(tmf[,2]~tmf$loc),lwd=2,lty=2)
legend('topleft',legend=c('HC','HE'),fill=crypcol[c(2,4)],bty='n')
dev.off()

###############
# Figure S3
###############

old<-read.delim('Rose2017meta.txt')
meta$oldsym<-old$propd[match(meta$tag,rownames(old))]

symfst<-read.delim('HE_symD_pool400.HE_symC_pool400.10kb.fst',col.names=c('chr','midpos','snites','fst'),row.names=1)
rownames(symfst)<-paste(symfst$chr,symfst$midpos-5000,sep=':')
symfst$outl<-rownames(symfst)%in%winpos[overplot]
symfst$chrnum<-(as.numeric(gsub('chr','',symfst$chr)))
symfst$chrnum[is.na(symfst$chrnum)]<-15
symfst$chrcol<-c(grey(0.2),grey(0.4))[1+symfst$chrnum%%2]
symfst$chrcol[symfst$chrnum==15]<-grey(0.7)
symfst$colchrom<-rainbow(11,v=0.8)[1+as.numeric(factor(symfst$chr))%%11]
symfst$chrcol[which(symfst$outl)]<-symfst$colchrom[which(symfst$outl)]
symfst<-symfst[order(symfst$chrnum),]

png('../Figures/FigureS5.png',width=6,height=3,res=300,units='in')
par(mar=c(3,3,2,1),mgp=c(1.7,0.1,0),las=1,tck=0.01,cex.axis=1,lwd=1,bty='l')
layout(t(as.matrix(c(1,1,2))))

plot(symfst$fst,col=symfst$chrcol,ylab=expression(F[ST]),xlab='Genomic position  (100 kb window, 10 kb step)')
points(which(symfst$outl),symfst$fst[symfst$outl],col=symfst$chrcol[symfst$outl])
text(43000,0.55,labels=c('xfSc0000366'),cex=0.7)
symfst[which.max(symfst$fst),]
covmat<-1000000*t(as.matrix(idx[,3:116]))/colSums(idx[,3:116])
title('A',adj=0)

plot(covmat[,idx[,1]=='xfSc0000366'],propd,ylab='Proportion sym. D',xlab='xfSc0000366 reads per million')
cor.test(covmat[,idx[,1]=='xfSc0000366'],propd)
title('B',adj=0)

dev.off()

####################

library(car)
lm.out<-lm(propd~pool+cryp,data=meta)
aov.out<-Anova(lm.out,type=2)

cat('###Type II ANOVA of symbiont association explained by pool, cryp###\n')
aov.out

cat('\n\n\n###eta squared proportions of symbiont variance explained by pool, cryp, residual###\n')
aov.out[[1]]/sum(aov.out[[1]])

cat('\n\n\n###linear model of symbiont association explained by pool, cryp###\n')
summary(lm.out)

cat('\n\n\n###Analysis of deviance of simple model of bleaching explained by cryp, symbiont, pool###\n')
glm.out<-step(glm(binbl~cryp+propd+pool,data=meta,family='binomial'))
aov.out<-Anova(glm.out,type=2)
aov.out

simpleAIC<-extractAIC(glm.out)


cat('\n\n\n###eta squared proportions of bleaching variance explained by cryp, propd, pool, residual###\n')
aov.out[[1]]/sum(aov.out[[1]])

cat('\n\n\n###glm model of bleaching explained by cryp, symbiont, pool###\n')
summary(glm.out)

cat('\n\n\n###test effect of L2###\n')
glm.out0<-(glm(binbl~cryp+propd+pool,data=meta,family='binomial'))
glm.out<-(glm(binbl~cryp+propd+pool+e2$vectors[,1],data=meta,family='binomial'))
aov.out<-anova(glm.out0,glm.out,test='Chisq')
aov.out

cat('\n\n\n###glm model of bleaching explained by cryp, symbiont, pool, L1, L2###\n')
glm.out<-step(glm(binbl~cryp+propd+pool+e1$vectors[,1]+e2$vectors[,1]+mitotype,data=meta,family='binomial'))
aov.out<-Anova(glm.out,type=2)
aov.out

cat('\n\n\n###relative likelihood of simple model relative to genomic model###\n')
genomicAIC<-extractAIC(glm.out)
genomicAIC[2]-simpleAIC[2]
exp((genomicAIC[2]-simpleAIC[2])/2)
