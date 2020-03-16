v1<-read.delim('~/Dropbox/GradWork/genomes2020/adiv1.picard')
v2<-read.delim('~/Dropbox/GradWork/genomes2020/iter1.picard')

plot(v1$PF_HQ_ALIGNED_READS/v1$TOTAL_READS,v2$PF_HQ_ALIGNED_READS/v2$TOTAL_READS)
abline(0,1)

plot(v1$PF_HQ_ERROR_RATE,v2$PF_HQ_ERROR_RATE,xlim=c(0.015,0.035),ylim=c(0.015,0.035))
abline(0,1)

mean(v1$PF_HQ_ALIGNED_BASES/400e6)