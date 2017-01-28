## Fig_S3
r12 = c(0.43, 0.22, 0.56, 0.53, 0.22, 0.27, 0.21, 0.20, 0.26, 0.30)
r21 = c(0.04, 0.07, 0.13, 0.04, 0.08, 0.04, 0.04, 0.05, 0.11, 0.04)
dage.ccf = c(4,5,6,3,5,7,4,8,6,6)
# r12 = r12/dage.ccf
# r21 = r21/dage.ccf

R12 = c(0.26, 0.91, 0.34, 0.39, 0.34, 0.52, 0.15, 0.05, 0.22, 0.35)
R21 = c(0.06, 0.00, 0.02, 0.02, 0.06, 0.04, 0.08, 0.24, 0.03, 0.04)
dage.cw = c(3,7,10,7,7,5,4,8,7,9)
# R12 = R12/dage.cw
# R21 = R21/dage.cw

plot(100*r12,100*r21,pch=19,xlim=c(0,20),ylim=c(0,10),
     xlab='% drift-CpGs advancing per year',
     ylab='% drift-CpGs retarding per year',cex.lab=1.3)
points(100*R12,100*R21,pch=19,col=2)
abline(0,1)
legend("topright",c("CCF","CW"),col=c(1,2),pch=c(19,19),bty='n',y.intersp = 2,inset=.1,cex=1.3)

## alternative plot

plot(c(dage.ccf,dage.cw[-8]),100*c(r12,R12[-8])-c(r21,R21[-8]),pch='+',col=1,xlab='time between biopsies',
     xlim=c(0,10), ylim=c(0,90), ylab='net % drift across threshold',cex=1.3,cex.lab=1.3)
# points(c(dage.ccf,dage.cw[-8]),c(r21,R21[-8]),cex=1.3)
abline(lm(100*(c(r12,R12[-8])-c(r21,R21[-8])) ~ 0 + c(dage.ccf,dage.cw[-8])),lwd=2)

