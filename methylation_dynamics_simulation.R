## time vector
a  = seq(0,3000,1)
len = length(a)

M = 50 # number of CpGs in region
n = 1000 # number of cells being arrayed
k = 1

## track trajectories every 50 time steps 
out = matrix(0,ncol=512,nrow=31)
traj = matrix(0,nrow=max(a)/10 + 1,ncol=10) # time trajectories of mean methylation levels for 10 cells

## set up starting state matrix
state0  = matrix(0,ncol=n,nrow=M)
for (j in 1:n) { state0[,j] = rbinom(M,size=1,prob=0.06) }  # about 6% for NS tissue 

## heterogeneous rate of methyl prob (gamma distributed) 
mu0 = .0001; var.mu0 =(.0002)^2
# mu0 = .0002; var.mu0 =(.0004)^2
shape = mu0*mu0/var.mu0
scale= mu0/shape
##  shape1 = .5
##  shape2 = shape1*(1-mu0)/mu0
##  mu = rbeta(M,shape1=shape1,shape2=shape2) #rexp(M,rate=1/mu0)
mu = rgamma(M,shape=shape,scale=scale)

## mu = .00005* rlnorm(M,meanlog = 0, sdlog = 1)
plot(density(mu),type='l')

## loop over time s
## mu depends on island methylation state in cell
Mu = Mu0 = matrix(rep(mu,n),nrow=M,ncol=n)
## initialize state
beta = state0

l=1
for (s in a) {
  for (i in 1:M) {
    # flip some zeros with rate mu (per year)
    beta[i,rbinom(n,size=1,prob=Mu[i,])==1]=1 # random flip up
    beta[i,rbinom(n,size=1,prob=0.3*Mu[i,])==1]=0 # random flip down
  }
  # check island level methylation (cell by cell, n cells)
  pm = apply(beta,2,mean) 
  for (j in 1:n) {
    tmp = Mu0[,j]*(1+100/(1+exp(-100*(pm[j]-0.2))))
    tmp[tmp>1]=1; Mu[,j] = tmp
  }
  
  if(s%%10==0) {
    traj[k,] = pm[1:10]
    k = k+1
  }
  
  if(s%%100==0) {
    dum = density(apply(beta,1,mean),from=0,to=1)
    out[l,]=out[l,]+dum$y
    l=l+1
    plot(dum,type='l',xlim=c(0,1),ylim=c(0,8))
  } #; Sys.sleep(.2)
  cat('time:',s,'\n')
}

# 
plot(dum$x,out[7,],type='l',xlab='beta-value',ylab='density',lwd=2,
     col=4,xlim=c(0,.9),ylim=c(0,9),main='')
# 
for (i in 7:30) {
   lines(dum$x,out[i,])
}
lines(dum$x,out[31,],lwd=2,col=2)
# lines(dum$x,out[13,],lwd=2,col='orange')

legend('topright',c("beta-value density (t=600)","beta-value density (t=3000)"),lwd=c(2,2),
       col=c(4,2),bty='n',y.intersp=1.5,inset = .15)

plot(traj[,1],type='l',ylim=c(0,.7),lwd=2,xlab='time/10',ylab='mean beta-value (50 CpGs)',cex.lab=1.3)
for(i in 2:10) lines(traj[,i],lwd=2, col=i)

