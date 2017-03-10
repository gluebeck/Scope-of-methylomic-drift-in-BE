#' Mvalue Ploting
#'
#' This function takes multiple M-values for two datasets, averages the m-value at each cpg for each dataset ,and connects
#' the two  sets of averages with a line at each cpg location.
#'
#' @param dat contains a dataframe with illumina 450k methylation data plots non-black cirlced points
#' @param squam second dataframe with same type as dat used to plot black circled points
#' @param IslName is the name of cpg island
#' @param Clockset dataframe including info on CpG clocks on Islands
#' @keywords methylation
#' @export
#' @examples
#' avgIslPlot_Mval(dat = beM.64, sqam = beM.64, IslName = CpGIslName, ClockSet = `CpGs.hypo.drift10.isl`, includeBar = TRUE, bottomleft = "")
# plots M-values for different people at different CpGs


avgIslPlot_Mval = function(dat, sqam, IslName, ClockSet, includeBar = TRUE, bottomleft = "") {

  #select  a particular island on a gene and display it's CpG's
  gofish = as.character(IslName)
  cpgn.GENE = manifestData[manifestData$Islands_Name==gofish,"Name"]#name of cpg
  gene.GENE = manifestData[manifestData$Islands_Name==gofish,"UCSC_RefGene_Name"]#gene name

  #gene GENE contains list of CPGs with gene PDE3A Pde3A
  pos.GENE = manifestData[manifestData$Islands_Name==gofish,"pos"]
  isl.GENE = manifestData[manifestData$Islands_Name==gofish,"Relation_to_Island"]
  refgrp.GENE = manifestData[manifestData$Islands_Name==gofish,"UCSC_RefGene_Group"]

  np = ncol(dat)
  ido = order(pos.GENE)
  print(ido)
  cpgn = cpgn.GENE[ido]
  labels = refgrp.GENE[ido]

  #number of cpgs
  ncpg = length(cpgn)

  # colors
  cols = rep(1,ncpg)
  cols[isl.GENE[ido]=="Island"]=2
  cols[grepl("Shore",isl.GENE[ido])]= 3
  cols[grepl("Shelf",isl.GENE[ido])]=4 #4
  #ids of CpGs
  ids1 = na.omit(match(cpgn,rownames(dat)))
  ids2 = na.omit(match(cpgn,rownames(sqam)))

  #cpg marker coordinates
  inc=(1:ncpg)[!is.na(match(cpgn,rownames(dat)))]

  plot(1:ncpg,rep(10,ncpg),pch=2,ylim=c(-6,8),xaxt='n',xlab="CpG",ylab="M-value")
  #plot(points(c(1,2),c(1,2),pch = "X"))
  axis(1,at=c(1:ncpg),labels=FALSE) # adds ticks to the side of the graph
  text(x=c(1:ncpg),y=rep(par("usr")[3],9), labels = labels, srt = 45, adj = c(1.2,2), xpd = TRUE, cex=.5)

  # inc
  for(i in 1:length(inc)) {
    #plots average m-values for points in
    mean = mean(dat[ids1[i],])
    sec_mean = mean(sqam[ids2[i],])
    points(inc[i],mean,col=cols[inc[i]],pch=19,cex=1.4)
    points(inc[i],sec_mean,col=cols[inc[i]],pch=19,cex=1.4)
    points(inc[i],sec_mean,col= "black",pch= 1,cex=1.4, lwd = 2)
    segments(inc[i], mean, inc[i] , sec_mean,col=cols[inc[i]],cex = 1.4)
  }

  ## gene memberships; indicate range by CpGs
  #makes top of of diagram red and grey strip
  if(includeBar){
    genes = unique(unlist(strsplit(gene.GENE,";")))
    for(i in 1:length(genes)) {
      i1 = i-1
      idum1 = (1:ncpg)[grepl(genes[i],gene.GENE[ido])]
      idum2 = (1:ncpg)[grepl("Island", isl.GENE[ido])]
      lines(sort(idum1),rep(7.5-0.8*i1,length(idum1)),lwd=10,col='darkgrey')
      lines(sort(idum2),rep(7.5,length(idum2)),lwd=5,col='red')
      text(mean(idum1),6.8-0.8*i1,genes[i],cex=1)
    }

    ## indicate clock CpGs
    idum = na.omit(match(ClockSet,cpgn))
    ids = na.omit(match(cpgn[idum],rownames(dat)))
    #if there is clock plot bar
    points(idum,rep(7.5,length(idum)),pch='|',cex=.7)

    # for(i in 1:length(idum)) {points(rep(idum[i],np),dat[ids[i],],col="orange",pch=19,cex=.2)}
    legend('bottomleft',bottomleft,bty='n',y.intersp=1.25,cex=1.2)
    return(list(cpgn=cpgn,genes=genes))
  }

}
