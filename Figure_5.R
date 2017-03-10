#### Figure 5

id.MM.BE = c(3, 18, 23, 25, 26, 30, 40, 42, 44, 45, 46, 59, 60, 63)
id.HM.BE = c(7, 10, 11, 15, 16, 19, 20, 21, 22, 28, 29, 35, 37, 38, 39, 51, 52, 53, 54, 55, 57, 61)

id.MM.EAC = c(2, 7, 11, 14, 15)
id.HM.EAC = c(1, 6, 17, 20, 22, 24)

par(mfrow=c(2,1),mar=c(0,5,4,4)+0.1,oma=c(0,0,0,0))

# plot(density(ilogit2(beM.64[cpgs.drift,-c(id.MM.BE,id.HM.BE)])),type='l',xlab='beta-value',ylab='Density', main='',
#     ylim=c(0,8),lwd=1,xlim=c(0,.9),lty=3)
plot(density(ilogit2(beM.64[cpgs.drift,-c(id.MM.BE,id.HM.BE)])),type='l',xaxt='n',
     xlab='',ylab='density', main='',
     ylim=c(0,8),lwd=1,xlim=c(0,.9),lty=3,cex.lab=1.4)
lines(density(ilogit2(beM.64[cpgs.drift,id.MM.BE])),lty=1,lwd=2)
lines(density(ilogit2(beM.64[cpgs.drift,id.HM.BE])),lty=2,lwd=2)
# lines(density(ilogit2(sqM.52[cpgs.drift,])),lwd=2,col='grey')
## BETRNet EAC
lines(density(ilogit2(eacM[cpgs.drift,-c(id.MM.EAC,id.HM.EAC)])),lty=3,lwd=1,col=2)
lines(density(ilogit2(eacM[cpgs.drift,id.MM.EAC])),lty=1,lwd=2,col=2)
lines(density(ilogit2(eacM[cpgs.drift,id.HM.EAC])),lty=2,lwd=2,col=2)
## Krause EAC
# lines(density(dat.EAC47[,-c(id.mm)]),lty=2,lwd=2,col=2)
# lines(density(dat.EAC47[,id.mm]),lty=1,lwd=2,col=2)
# lines(density(dat.EAC47[,id.hm]),lty=3,col=2)
## SQ
# lines(density(ilogit2(sqM.52[cpgs.drift,])),col="grey",lwd=2)
# abline(v=median(ilogit2(sqM.52[cpgs.drift,])),col="grey")

legend('topright',c("BE-L","BE-I","BE-H","","EAC-L","EAC-I","EAC-H"),
       lty=c(1,3,2,1,1,3,2),
       col=c(1,1,1,0,2,2,2),lwd=c(2,1,2,0,2,1,2),bty='n',y.intersp = 1.3,inset=.1)
# legend('topright',c("SQ","BE-L","BE-H","EAC-L","EAC-H"),
#        lty=c(1,1,2,1,2),
#        col=c("grey","black","black","red","red"),lwd=c(2,2,2,2,2),bty='n',y.intersp = 1.3,inset=.1)
text(0.2,7,"A",cex=1.8)

par(mar=c(4,5,0,4)+0.1)
plot(x.beta,out[7,],type='l',xlab='beta-value',ylab='density',lwd=2,
     col=4,xlim=c(0,.9),ylim=c(0,8),main='',cex.lab=1.4)
# 
for (i in 7:30) {
  lines(x.beta,out[i,])
}
lines(x.beta,out[31,],lwd=2,col="orange")
legend('topright',c("beta-value density (t=600)","beta-value density (t=3000)"),lwd=c(2,2),
       col=c("blue","orange"),bty='n',y.intersp=1.5,inset = .15)
text(0.2,7,"B",cex=1.8)


