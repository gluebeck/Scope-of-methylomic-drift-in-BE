#### Figure 4 (finestructure)
#### subsample SQ, BE, Dysp, and EAC drift using drift CpGs (CpGs.hypo.drift5.isl)

#### BETRNet methylation data: sqM.52, beM.64, lgdM, hgdM, eacM, and eacM.TCGA (TCGA)
#### patient info (age, sex): pheno (phenodata) 

len = length(CpGs.hypo.drift5.isl)

# subsample 1000 out of set CpGs.hypo.drift5.isl
idx.ran = sample(1:len,1000,replace=F)
cpgs.drift.ran1k = CpGs.hypo.drift5.isl[idx.ran]

ids = na.omit(match(cpgs.drift.ran1k,rownames(sqM.52)))
dat.1000.sq52 = sqM.52[ids,]
dat.1000.be64 = beM.64[ids,]
dat.1000.dysp = cbind(lgdM[ids,],hgdM[ids,])
dat.1000.eac  = eacM[ids,]
ids = na.omit(match(cpgs.drift.ran1k,rownames(eacM.TCGA)))
iex = c(1:1000)[is.na(match(cpgs.drift.ran1k,rownames(eacM.TCGA)))]
dat.1000.tcga = eacM.TCGA[ids,]

## nSQ mean ranks for drift Cpgs (reflecting methylation fine-structure)
dat = dat.1000.sq52
# dat = matrix(runif(1000*52),ncol=52,nrow=1000)

aux = matrix(0,nrow=nrow(dat),ncol=ncol(dat))
for(i in 1:ncol(dat)) {# ids = sample(1:2511,2511,replace=F)
  idum = rank(dat[,i])
  aux[,i]=idum # sample(1:2511,2511,replace=F) #idum
} #order(dat.2516.be64[,i])}
tmp = apply(aux,1,mean)/nrow(dat)
rank.sq = tmp
rank.sq.tcga = tmp[-iex]

## BE mean ranks for drift Cpgs (reflecting methylation fine-structure)
dat = dat.1000.be64
aux = matrix(0,nrow=nrow(dat),ncol=ncol(dat))
for(i in 1:ncol(dat)) {# ids = sample(1:2511,2511,replace=F)
  idum = rank(dat[,i])
  aux[,i]=idum # sample(1:2511,2511,replace=F) #idum
} 
tmp = apply(aux,1,mean)/nrow(dat)
rank.BE = tmp

## Dyspl mean ranks for drift Cpgs (reflecting methylation fine-structure)
dat = dat.1000.dysp
aux = matrix(0,nrow=nrow(dat),ncol=ncol(dat))
for(i in 1:ncol(dat)) {# ids = sample(1:2511,2511,replace=F)
  idum = rank(dat[,i])
  aux[,i]=idum # sample(1:2511,2511,replace=F) #idum
} 
tmp = apply(aux,1,mean)/nrow(dat)
rank.dysp = tmp

## TCGA mean ranks for drift Cpgs (reflecting methylation fine-structure)
dat = dat.1000.tcga
aux = matrix(0,nrow=nrow(dat),ncol=ncol(dat))
for(i in 1:ncol(dat)) {# ids = sample(1:2511,2511,replace=F)
  idum = rank(dat[,i])
  aux[,i]=idum # sample(1:2511,2511,replace=F) #idum
} 
tmp = apply(aux,1,mean)/nrow(dat)
rank.tcga = tmp

## EAC (BETRNet) mean ranks for drift Cpgs (reflecting methylation fine-structure)
dat = dat.1000.eac
aux = matrix(0,nrow=nrow(dat),ncol=ncol(dat))
for(i in 1:ncol(dat)) {# ids = sample(1:2511,2511,replace=F)
  idum = rank(dat[,i])
  aux[,i]=idum # sample(1:2511,2511,replace=F) #idum
} 
tmp = apply(aux,1,mean)/nrow(dat)
rank.EAC = tmp


#### plot correlation with NS mean rank order of 1000 randomly selected drift CpGs

dat = dat.1000.tcga  #dat.1000.dysp # dat.1000.be64 
# idum1=which(pheno$Patient_Dx[ids.BE37]=="HGD")
# ind.progr = c(rep("NP",64),rep("P",37))
# col = c(rep(1,64),rep(2,37),rep(3,42)); col[64+idum]=4
aux = pval = numeric(ncol(dat))
for(i in 1:ncol(dat)) {
  tmp = rank(dat[,i])
  aux[i]=cor(tmp,rank.sq[-iex])
  pval[i] = cor.test(tmp,rank.sq[-iex])$p.value
}
cor.rank.tcga = aux
# cor.rank.eac = aux
# cor.rank.dysp = aux
# cor.rank.be64 = aux

## compute mean differential drift for cpgs.drift.ran1k
ids = na.omit(match(cpgs.drift.ran1k,CpGs))
aux = eacM[cpgs.drift,] # cbind(lgdM[cpgs.drift,],hgdM[cpgs.drift,]) # beM.64[cpgs.drift,]

ids = na.omit(match(cpgs.drift.ran1k[-iex],CpGs))
aux = eacM.TCGA[cpgs.drift.ran1k[-iex],]
for(i in 1:997) {
  # aux[i,] = aux[i,]-ancova.SQ$eps[ids[i]]-ancova.SQ$rate[ids[i]]*age.BE64
  # aux[i,] = aux[i,]-ancova.SQ$eps[ids[i]]-ancova.SQ$rate[ids[i]]*c(age.LGD,age.HGD)
  # aux[i,] = aux[i,]-ancova.SQ$eps[ids[i]]-ancova.SQ$rate[ids[i]]*age.EAC
  aux[i,] = aux[i,]-ancova.SQ$eps[ids[i]]-ancova.SQ$rate[ids[i]]*age.EAC.TCGA
}
# mean.drift.be64 = apply(aux,2,mean)
# mean.drift.dysp = apply(aux,2,mean)
# mean.drift.eac = apply(aux,2,mean)
mean.drift.tcga = apply(aux,2,mean)

# cor.test(mean.drift.be64,cor.rank.be64)
# cor.test(mean.drift.dysp,cor.rank.dysp)
# cor.test(mean.drift.eac,cor.rank.eac)
# cor.test(mean.drift.tcga,cor.rank.tcga)

# mean beta-value drift 
mean.drift.be64 = apply(ilogit2(dat.1000.be64),2,mean)
mean.drift.dysp = apply(ilogit2(dat.1000.dysp),2,mean)
mean.drift.eac =  apply(ilogit2(dat.1000.eac),2,mean)
mean.drift.tcga = apply(ilogit2(dat.1000.tcga),2,mean)

par(mar=c(4,5,2,2))

plot(mean.drift.be64,cor.rank.be64,pch=19,cex=.5,xlim=c(0.05,0.65),ylim=c(0,.9),
     xlab='drift (beta-value)',ylab='Pearson correlation (fine-structure)')
lines(smooth.spline(mean.drift.be64,cor.rank.be64,spar=1.1),col=1,lwd=2)
points(mean.drift.dysp,cor.rank.dysp,pch=19,cex=.5,col=2)
lines(smooth.spline(mean.drift.dysp,cor.rank.dysp,spar=1),col=2,lwd=2,lty=1)
# iex2 = which.min(mean.drift.eac)
points(mean.drift.eac,cor.rank.eac,pch=19,cex=.5,col='blue')
lines(smooth.spline(mean.drift.eac,cor.rank.eac,spar=1.1),col='blue',lwd=2)
points(mean.drift.tcga,cor.rank.tcga,pch=17,cex=.5,col=4)
lines(smooth.spline(mean.drift.tcga,cor.rank.tcga,spar=1.1),col=4,lwd=2,lty=3)

legend('topright',c("BE (n=64)","Dyspl. (n=43)","EAC (n=24)","TCGA (n=87)"),
       col=c(1,2,4,4),lty=c(1,3,1,3),pch=c(19,19,19,17),bty='n', y.intersp = 1.4,
       inset=.06)
