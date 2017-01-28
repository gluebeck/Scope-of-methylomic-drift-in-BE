## plot autocorr and island density  
out=CpGIsland_AutoCorr(dat = beM.64, IslNames = Isl.Name, incr = 10)

## plot pairwise correlation with index at Island
ids = 1:500
plot(ids*10,out.AutoCorr.Island.hypo$corrDensity[ids],pch=19,col=1,cex=.2,ylim=c(-.2,1),
     xlab='distance (bases)',ylab='correlation') # ,main='within Island CpG-CpG correlations')
lines(loess.smooth(ids*10,out.AutoCorr.Island.hypo$corrDensity[ids],span=.2),
      lwd=2,col=1)
points(ids*10,out.AutoCorr.Island.hypo.drift$corrDensity[ids],pch=19,col=2,cex=.2)
lines(loess.smooth(ids*10,out.AutoCorr.Island.hypo.drift$corrDensity[ids],span=.2),
      lwd=2,col=2)
# tmp = out.AutoCorr.Island.hypo$islDensity[ids]; tmp=tmp/1/max(tmp)
# lines(loess.smooth(ids*10,tmp,span=.2),type='l',lty=3,lwd=2)
# tmp = out.AutoCorr.Island.hypo.drift$islDensity[ids]; tmp=tmp/1/max(tmp)
# lines(loess.smooth(ids*10,tmp,span=0.2),type='l',lty=3,col=2,lwd=2)
# legend('topleft',c("B"),bty='n',cex=2)
# abline(v=2500,lty=3)
polygon(x=c(2000,2500,2500,2000),y=c(-.5,-0.5,1.2,1.2),border=NA,density=10)
text(3000,0.9,"Shelf",cex=1.2)
text(1000,0.9,"Island+Shore",cex=1.2)
legend('bottomleft',inset=.05,c("static islands","drift islands"),lwd=c(2,2),col=1:2,bty='n',
       y.intersp = 1.5,cex=.9)
# legend('topright',c("density static","density drift"),col=1:2,lty=c(3,3),
#       lwd=c(2,2),bty = 'n', y.intersp = 1.3,cex=1.2)
abline(h=0)
