CpGisl.level.info = function(set = ILS.hypo.drift5, isls = Islands.hypo.drift5, 
                             cpgs = CpGs.hypo.drift5.isl, dat, full = T, prom = T) {
  # input set: island names of interest
  # input isls and cpgs: vectors of the same length with CpG names and associated island names
  # input data: methylation data 
  # input str1 and str2: optional, to constrain selection of CpGs for mean calculation
  
  len = length(set) # number of Islands
  Mvals = sd.Mvals = matrix(0,ncol=ncol(dat),nrow=len)
  bslope = sd.bslope = numeric()
  # requires input: aux1 =  slope.ccf, aux2 = CpGs.ccf
  genes = list()
  islands = vector()
  Prom = rep(T,len)
  
  ii = 1
  for (i in 1:len) {
    
    if(full==T) {
      cpgsi = manifestData$Name[manifestData$Islands_Name == set[i]]  # cpgs on one of the islands
    } else {
      cpgsi = cpgs[isls == set[i]]  # cpgs on one of the islands
    }
    
    # remove shelves !!! already checked
    dum = manifestData[cpgsi,"Relation_to_Island"]
    cpgsi = cpgsi[dum!="N_Shelf" & dum!="S_Shelf"]  # & dum!="N_Shore" & dum!="S_Shore"] 

    # check exclusions (only informative islands etc. )
    dum = manifestData[cpgsi,"UCSC_RefGene_Group"]
    tmp = cpgsi[grepl("TSS",dum)] # & !grepl("Body",dum) & !grepl("3'UTR",dum)]
    if(prom == T) {
      if(length(tmp)==0) {next} else {cpgsi = tmp}
    }
    if(length(tmp)==0) Prom[ii]=FALSE
    
    # see which cpgs are in dat, if any .... 
    idum = na.omit(match(cpgsi,rownames(dat)))
    len.idum = length(idum)
    if(len.idum == 0) {next}
    if(len.idum ==1) {
      Mvals[ii,]=dat[idum,]; sd.Mvals[ii,]=0
    } else {
      Mvals[ii,] = apply(dat[idum,],2,mean,na.rm=T); sd.Mvals[ii,] = apply(dat[idum,],2,sd,na.rm=T)
    }
  
    ## get b-slope from CCF data
    idum = na.omit(match(cpgsi,CpGs.ccf))
    bslope[ii] = mean(slope.ccf[idum])
    # idum = na.omit(match(cpgsi,CpGs))
    # bslope[ii] = mean(ancova.BE$rate[idum])
    # sd.bslope[ii] = sd(aux1[idum])
    
    # get gene name for island
    genes[[ii]] = dum = unique(unlist(strsplit(manifestData[cpgsi,"UCSC_RefGene_Name"],split=';')))
    islands[ii] = set[i]
    
    ii = ii+1
  }
  
  rownames(Mvals) = rownames(sd.Mvals) = set
  colnames(Mvals) = colnames(sd.Mvals) = colnames(dat)
  Mvals = Mvals[1:(ii-1),]; sd.Mvals = sd.Mvals[1:(ii-1),]
  bslope = bslope[1:(ii-1)]
  
  tmp = lapply(genes,strsplit,';')
  tmp = lapply(tmp,unlist)
  genes = lapply(tmp,unique)
  Prom = Prom[1:(ii-1)]

  # return(list(Mvals=Mvals, sd.Mvals=sd.Mvals, bslope=bslope, sd.bslope=sd.bslope, genes=genes))
  return(list(Mvals=Mvals, sd.Mvals=sd.Mvals, genes=genes, islands=islands, prom=Prom, bslope=bslope))
}
