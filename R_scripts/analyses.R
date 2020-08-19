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

globalfst<-read.delim('global_fst.txt',header=F,sep=' ')
globaltree<-getTrees(globalfst[,2])
globaltree$tip.label<-c('HA','HC','HD','HE')
plot(globaltree,type='u',cex=1, edge.width=2,lab4ut='axial') 
segments(.05,0,.1,0,lwd=2)
text(0.075,0.01,label=expression(italic(T) == 0.05))
title('C',adj=0)


#D
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
title('D',adj=0)
par(xpd=F)
text(0,-10,label='HG1')
text(30,-10,label='HG2')


#E
plot(jitter(as.numeric(factor(meta$cryp))),meta$propd,axes=F,xlab='Species',ylab=expression(paste("Proportion",italic(" Durusdinium")," symbionts")),type='n',cex=1,ylim=c(0,1),las=1,cex.axis=1)
axis(1,at=1:4,labels=c('HA','HC','HD','HE'))
axis(2)
box(bty='l')
abline(h=0.5,lwd=2,lty=2)
points(jitter(as.numeric(factor(meta$cryp))),meta$propd,col=crypcols,cex=1)

title('E',adj=0)

#F
getPBS<-function(focal.ref,focal.out,ref.out){
	T1=-log(1-focal.ref)
	T2=-log(1-focal.out)
	T3=-log(1-ref.out)
	PBS=(T1+T2-T3)/2
	PBS[PBS<0]<-0
	return(PBS)
}


sns<-fread('snpfst_sns_nowarning.txt.gz')
colnames(sns)<-c('SNP','HA.HC','HA.HD','HA.HE','HC.HD','HC.HE','HD.HE','gene','Eff')
sns$HA<-getPBS(sns$HA.HE,sns$HA.HD,sns$HD.HE)
sns$HC<-getPBS(sns$HC.HD,sns$HC.HE,sns$HD.HE)
sns$HD<-getPBS(sns$HC.HD,sns$HD.HE,sns$HC.HE)
sns$HE<-getPBS(sns$HC.HE,sns$HD.HE,sns$HC.HD)
HA<-fisher.test(sns$HA> quantile(sns$HA,0.95),sns$Eff=='missense')
HC<-fisher.test(sns$HC> quantile(sns$HC,0.95),sns$Eff=='missense')
HD<-fisher.test(sns$HD> quantile(sns$HD,0.95),sns$Eff=='missense')
HE<-fisher.test(sns$HE> quantile(sns$HE,0.95),sns$Eff=='missense')

tables2<-sns[sns$HE> quantile(sns$HE,0.95),c('SNP','HA.HC','HA.HD','HA.HE','HC.HD','HC.HE','HD.HE','HA','HC','HD','HE','gene','Eff')]
tables2<-tables2[order(gsub(':.*','',gsub('chr','',tables2$SNP))),]
colnames(tables2)<-c('SNP','HA.HC.FST','HA.HD.FST','HA.HE.FST','HC.HD.FST','HC.HE.FST','HD.HE.FST','HA.PBS','HC.PBS','HD.PBS','HE.PBS','gene','Eff')
write.table(tables2,file='../Tables/TableS2.txt',quote=F,row.names=F,sep='\t')

ors<-cbind(HA$estimate,HC$estimate,HD$estimate,HE$estimate)
errs<-rbind(HA$conf.int,HC$conf.int,HD$conf.int,HE$conf.int)
ortab<-matrix(0,4,4)
ortab[1,1]<-ors[1]
ortab[2,2]<-ors[2]
ortab[3,3]<-ors[3]
ortab[4,4]<-ors[4]
bp<-barplot(ortab,ylim=c(0,1.2),col=crypcol,border=NA,names.arg=c('HA','HC','HD','HE'),las=1,ylab='NS vs. S outlier odds ratio',cex.names=1)
segments(bp,errs[,1],bp,errs[,2],lwd=2)
abline(h=1,lwd=2,lty=3)
title('F',adj=0)

quantile(sns$HE,0.95)
table(sns$HE>quantile(sns$HE,0.95))
length(unique(unlist(sns[sns$HE>quantile(sns$HE,0.95),'gene'])))

quantile(sns$HE,0.99)
table(sns$HE>quantile(sns$HE,0.99))
length(unique(unlist(sns[sns$HE>quantile(sns$HE,0.99),'gene'])))


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

###########
# Figure S1
############
Fst10kb<-data.frame(row.names=winpos,HA.HC=HA.HC$num/HA.HC$den)
Fst10kb$HA.HD<-HA.HD$num/HA.HD$den
Fst10kb$HA.HE<-HA.HE$num/HA.HE$den
Fst10kb$HC.HD<-HC.HD$num/HC.HD$den
Fst10kb$HC.HE<-HC.HE$num/HC.HE$den
Fst10kb$HD.HE<-HD.HE$num/HD.HE$den
Fst10kb<-na.omit(Fst10kb)
Fst10kb<-as.matrix(Fst10kb)
trees<-apply(Fst10kb,1,getTrees)

getTopo<-function(tree){
	tip<-which(tree$tip.label=='HE')
	tree$tip.label[tree$edge[(tree$edge[,1]==tree$edge[tree$edge[,2]==tip,1]) & tree$edge[,2]<5 & tree$edge[,2]!=tip,2]]
}

drawTree<-function(x,y,lab){
	segments(x-0.15,y,x+0.15,y)
	segments(x-0.15,y,x-0.3,y+0.05)
	segments(x-0.15,y,x-0.3,y-0.05)
	segments(x+0.15,y,x+0.3,y+0.05)
	segments(x+0.15,y,x+0.3,y-0.05)
	text(c(x-0.42,x-0.42,x+0.42,x+0.42),c(y-0.05,y+0.05,y-0.05,y+0.05),labels=lab)
}

topos<-sapply(trees,getTopo)

png('../Figures/FigureS1.png',width=3,height=4,units='in',res=300)
par(mfrow=c(1,2),tck=-0.01,mar=c(1,3,2,1),mgp=c(1.5,0.1,0))
barplot(as.matrix(table(topos)),ylab='Counts')
par(tck=-0.01,mar=c(1,1,2,1),mgp=c(1.5,0.1,0))
plot(1,type='n',axes=F,xlab='',ylab='',xlim=c(0,1),ylim=c(0,1))
drawTree(.5,.95,c('HA','HC','HD','HE'))
drawTree(.5,.74,c('HA','HD','HC','HE'))
drawTree(.5,.3,c('HD','HC','HA','HE'))
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

getPerm<-function(df1,df2,df3,pis,edges,reps=10,k=10){
	perms<-replicate(reps,{
		shuf<-shuffle(length(pis),control=how(blocks=cut(pis,quantile(pis,seq(0,1,0.1),na.rm=T))))
		getPBSwin(df1[shuf,],df2[shuf,],df3[shuf,],k=k,edges=edges)
	})
	return(perms)
}

getPBSwin<-function(focal.ref,focal.out,ref.out,k=10,edges){
	T1=-log(1-rollapply(focal.ref$num,FUN=mean,width=k,na.rm=T,fill=NA)/rollapply(focal.ref$den,FUN=mean,width=k,na.rm=T, fill=NA))
	T2=-log(1-rollapply(focal.out $num,FUN=mean,width=k,na.rm=T,fill=NA)/rollapply(focal.out $den,FUN=mean,width=k,na.rm=T, fill=NA))
	T3=-log(1-rollapply(ref.out $num,FUN=mean,width=k,na.rm=T,fill=NA)/rollapply(ref.out $den,FUN=mean,width=k,na.rm=T, fill=NA))
	PBS=(T1+T2-T3)/2
	PBS[edges]<-NA
	return(PBS)
}

PBS<-data.frame(row.names=paste(HA.HC$chr,HA.HC$start,sep=':'),chr=HA.HC$chr)
PBS$chr<-factor(PBS$chr,levels=c(paste('chr',1:14,sep=''),'unplaced'))
PBS$chr[is.na(PBS$chr)]<-'unplaced'
PBS$HA<-getPBSwin(HA.HE, HA.HD, HD.HE,edges=edges)
PBS$HC<-getPBSwin(HC.HD, HC.HE, HD.HE,edges=edges)
PBS$HD<-getPBSwin(HC.HD, HD.HE, HC.HE,edges=edges)
PBS$HE<-getPBSwin(HA.HE, HD.HE, HA.HD,edges=edges)
PBS$chrcol<-c(grey(0.2),grey(0.4))[1+as.numeric(PBS$chr)%%2]
PBS$chrcol[PBS$chr=='unplaced']<-grey(0.7)
PBS$chrcolsave<-PBS$chrcol

# permHA<-getPerm(HA.HE, HA.HD, HD.HE,pis=meanPi,edges=edges,reps=10,k=10)
# permHC<-getPerm(HC.HD, HC.HE, HD.HE,pis=meanPi,edges=edges,reps=10,k=10)
# permHD<-getPerm(HC.HD, HD.HE, HC.HE,pis=meanPi,edges=edges,reps=10,k=10)
# permHE<-getPerm(HA.HE, HD.HE, HA.HD,pis=meanPi,edges=edges,reps=10,k=10)
# ecdf.HA<-ecdf(permHA)
# ecdf.HC<-ecdf(permHC)
# ecdf.HD<-ecdf(permHD)
# ecdf.HE<-ecdf(permHE)
# save(perms,ecdf.outl,permHA,permHC,permHD,permHE,ecdf.HA,ecdf.HC,ecdf.HD,ecdf.HE,file='perms.Rdata')

load('perms.Rdata')
outlpadjHA<-p.adjust(2*pmin((1-ecdf.outl(PBS$HA)),ecdf.outl(PBS$HA)),method='BH')
outlpadjHC<-p.adjust(2*pmin((1-ecdf.outl(PBS$HC)),ecdf.outl(PBS$HC)),method='BH')
outlpadjHD<-p.adjust(2*pmin((1-ecdf.outl(PBS$HD)),ecdf.outl(PBS$HD)),method='BH')
outlpadjHE<-p.adjust(2*pmin((1-ecdf.outl(PBS$HE)),ecdf.outl(PBS$HE)),method='BH')

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
HEwin<-t(sapply(winpos[overplot],function(v) c(strsplit(v,split=':')[[1]])))
HEbed<-cbind(HEwin[,1],as.numeric(HEwin[,2])-50000,as.numeric(HEwin[,2])+50000)
write.table(HAbed,row.names=F,col.names=F,quote=F,file='HA_outlier_regions.bed',sep='\t')
write.table(HCbed,row.names=F,col.names=F,quote=F,file='HC_outlier_regions.bed',sep='\t')
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
text(19600,1.75,labels='Locus 1')
text(34100,1.02,labels='Locus 2')

plot(HA.HE.rolldxy[nsites>50000],HA.HE.rollpi[nsites>50000],col=PBS$chrcol[nsites>50000],cex=0.5,xlab=expression(pi[b]),ylab=expression(pi[w]))
points(HA.HE.rolldxy[overplot],HA.HE.rollpi[overplot],col=PBS$chrcol[overplot],cex=0.5)
title('B',adj=0)

plot(HC.HE.rolldxy[nsites>50000],HC.HE.rollpi[nsites>50000],col=PBS$chrcol[nsites>50000],cex=0.5,xlab=expression(pi[b]),ylab=expression(pi[w]))
points(HC.HE.rolldxy[overplot],HC.HE.rollpi[overplot],col=PBS$chrcol[overplot],cex=0.5)
title('C',adj=0)

plot(HD.HE.rolldxy[nsites>50000],HD.HE.rollpi[nsites>50000],col=PBS$chrcol[nsites>50000],cex=0.5,xlab=expression(pi[b]),ylab=expression(pi[w]))
points(HD.HE.rolldxy[overplot],HD.HE.rollpi[overplot],col=PBS$chrcol[overplot],cex=0.5)
title('D',adj=0)

dev.off()


################################
#Figure S2
###############################

png('../Figures/FigureS2.png',width=6,height=6,res=300,units='in')
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
PBS$HEwin<-getPBSwin(HA.HE, HD.HE, HA.HD,edges=edges,k=1)

#########
sns2<-fread('locus1_sns.txt')
sns2$Eff<-'synonymous'
sns2$Eff[grepl('missense_variant',sns2[,V9])]<-'missense'
colnames(sns2)<-c('chr','pos','HA.HC','HA.HD','HA.HE','HC.HD','HC.HE','HD.HE','ANN','Eff')
sns2$HE<-getPBS(sns2$HC.HE,sns2$HD.HE,sns2$HC.HD)

table(sns2$Eff)
table(sns2$Eff,sns2$HE>quantile(sns$HE,0.95))
fisher.test(sns2$Eff,sns2$HE>quantile(sns$HE,0.95))

########

png('../Figures/Figure3.png',width=6,height=6,res=300,units='in')
par(mar=c(3,3,2,1),mgp=c(1.7,0.1,0),las=1,tck=0.01,cex.axis=1,lwd=1,bty='l',mfrow=c(2,2))

L1<-PBS[PBS$chr=='chr7'&PBS$start>=20250000&PBS$start<20670000,]
plot(L1$start,L1$HEwin,type='l',ylab='HE PBS',xlab='Chromosome 7 position (Mb)',ylim=c(0,4.5),xlim=c(20250000, 20670000),axes=F)
axis(2,at=c(0:3))
axis(1,at=seq(20300000,20600000,by=100000),labels=seq(20300000,20600000,by=100000)/1000000)
g7<-genes[genes[,1]=='chr7',]
g7<-g7[order(g7[,4]),]
g7$height<-rep(c(2.6,2.8,3.2,3.4,3.6,3.8,4,4.2),10)[1:nrow(g7)]
segments(g7[,4],g7$height,g7[,5],g7$height)
text((g7[,4]+g7[,5])/2,g7$height+0.1,labels=g7$gene2,cex=0.5,font=1+as.numeric(g7$de))
title('A',adj=0)

L2<-PBS[HA.HC$chr=='Sc0000015'&PBS$start>=2350000&PBS$start<=2650000,]
plot(L2$start,L2$HEwin,type='l',ylab='HE PBS',xlab='Sc0000015 position (Mb)',ylim=c(0,4.5),xlim=c(2300000, 2650000),axes=F)
axis(2,at=c(0:3))
axis(1,at=seq(2340000, 2650000,by=100000),labels=seq(2340000, 2650000,by=100000)/1000000)
gs<-genes[genes[,1]=='Sc0000015',]
gs<-gs[order(gs[,4]),]
gs$height<-rep(seq(1.1,1.6,.1),10)[1:nrow(gs)]
gs$height<-rep(c(2.4,2.6,2.8,3.2,3.4,3.6,3.8,4,4.2),10)[1:nrow(gs)]
segments(gs[,4],gs$height,gs[,5],gs$height)
text((gs[,4]+gs[,5])/2,gs$height+0.1,labels=gs$gene2,cex=0.5,font=1+as.numeric(gs$de))
title('B',adj=0)

bleachcol<-crypcols
bleachcol[meta$binbl==1]<-grey(0.9)
plot(jitter(as.numeric(meta$cryp)),e1$vectors[,1],bg=bleachcol,pch=21,axes=F,ylab='Locus 1 PC1 (62% of variation)',xlab='Cryptic species')
axis(1,at=1:4,labels=c('HA','HC','HD','HE'),lwd=0)
axis(2,at=c(-0.1,0,0.1))
title('C',adj=0)
legend('top',legend='bleached',fill=grey(0.95),bty='n')

plot(jitter(as.numeric(meta$cryp)),e2$vectors[,1],bg=bleachcol,pch=21,axes=F,ylab='Locus 2 PC1 (70% of variation)',xlab='Cryptic species')
axis(1,at=1:4,labels=c('HA','HC','HD','HE'),lwd=0)
axis(2,at=c(-0.1,0,0.1))
title('D',adj=0)

legend('top',legend='bleached',fill=grey(0.95),bty='n')
arrows(1.6,0.03,1.3,0.05,length=0.05,lwd=1)
arrows(1.6,0.02,1.3,0,length=0.05,lwd=1)
arrows(2.6,-0.03,2.3,-0.03,length=0.05,lwd=1)
arrows(3.6,-0.05,3.3,-0.05,length=0.05,lwd=1)
dev.off()


################
# Figure S3
##################

tmf<-read.delim('tmf_counts.txt',header=F)
tmf$loc<-e2$vectors[match(tmf[,1],as.character(meta$tag)),1]
tmf$cryp<-meta$cryp[match(tmf[,1],as.character(meta$tag))]
cor(tmf$loc,tmf[,2],use='p')
anova(lm(tmf[,2]~tmf$cryp+tmf$loc))

png('../Figures/FigureS3.png',width=3,height=3,res=300,units='in')
par(mar=c(3,2,1,1),tck=-0.01,mgp=c(1.1,.1,0),bty='l')
plot(tmf$loc,tmf[,2],bg='black',col=crypcols[match(tmf[,1],as.character(meta$tag))],pch=19,xlab='PC1 Locus 2',ylab='TMF norm. counts')
abline(lm(tmf[,2]~tmf$loc),lwd=2,lty=2)
legend('topleft',legend=c('HC','HE'),fill=crypcol[c(2,4)],bty='n')
dev.off()

###############
# Figure S4
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

png('../Figures/FigureS4.png',width=6,height=3,res=300,units='in')
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

sink('../Tables/TableS4.txt')

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

sink()