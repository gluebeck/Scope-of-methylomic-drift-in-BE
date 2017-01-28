# Figure S2  finestructure on MGMT

CpGIslName = "chr10:131264948-131265710"
tmp = avgIslPlot_Mval(dat = beM, sqam = sqM, IslName = CpGIslName, ClockSet = `CpGs.hypo.drift1.isl`, includeBar = TRUE, bottomleft = "")
cpgs.fine = tmp$cpgn[17:21]

boxplot(t(sqM.52[cpgs.fine,]),pch=19,cex=.5,xlab='island-CpGs',ylab='M-value',xaxt='n',boxwex=.3,
        notch=T,axes=F,cex.lab=1.3)
text(x=c(1:5),y=rep(par("usr")[3],9)+1.1, labels = cpgs.fine, srt = 45, adj = c(1.2,2), xpd = TRUE, cex=.7)
axis(2,at=c(-8:0))
text(3,-.5,expression(paste(bold("MGMT"))),cex=1.4) #,'\n',"island: chr10:131264948-13126571")))

aux = funM[cpgs.fine,]
points(rep(1,12),aux[1,],pch=19,cex=.4,col='red')
points(rep(2,12),aux[2,],pch=19,cex=.4,col='red')
points(rep(3,12),aux[3,],pch=19,cex=.4,col='red')
points(rep(4,12),aux[4,],pch=19,cex=.4,col='red')
points(rep(5,12),aux[5,],pch=19,cex=.4,col='red')
legend('topright',c("SQ (n=52)","fundus (n=12)"),col=c("black","red"),pch=c(19,19),
       bty='n',y.intersp=1.4,inset=0.1,cex=1.3)

# legend('topright',c("MGMT"),bty='n',inset=.1)

# > ilogit2(median(sqM.52[cpgs.fine[1],]))
# [1] 0.193793
# > ilogit2(median(sqM.52[cpgs.fine[2],]))
# [1] 0.05031279
# > ilogit2(median(sqM.52[cpgs.fine[3],]))
# [1] 0.03537204
# > ilogit2(median(sqM.52[cpgs.fine[4],]))
# [1] 0.01090893
# > ilogit2(median(sqM.52[cpgs.fine[5],]))
# [1] 0.06567929
