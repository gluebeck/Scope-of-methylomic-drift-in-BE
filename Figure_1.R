## Figure 1 
## island-level: TSS
## requires vector of Island names (e.g. ILS.hypo, ILS.hypo.drift1, ...) and manifest  
## requires vector with CpG probes (here they are island-related excluding shores) 
TSS = grepl("TSS",manifestData[CpGs.hypo.isl,"UCSC_RefGene_Group"])
TSS.drift = grepl("TSS",manifestData[CpGs.hypo.drift1.isl,"UCSC_RefGene_Group"])
Body = grepl("Body",manifestData[CpGs.hypo.isl,"UCSC_RefGene_Group"])
Body.drift = grepl("Body",manifestData[CpGs.hypo.drift1.isl,"UCSC_RefGene_Group"])

i=0
for (isl in ILS.hypo) {
  cpgi = manifestData[manifestData$Islands_Name==isl,"UCSC_RefGene_Group"]
  dum = manifestData[manifestData$Islands_Name==isl,"Relation_to_Island"]
  ind = grepl("Island",dum) | grepl("Shore",dum)
  cpgi = cpgi[ind]
  # if(sum(grepl("TSS",cpgi))>0) i=i+1
  # if(sum(grepl("Body",cpgi))>0 & !sum(grepl("TSS",cpgi))>0) i=i+1
  if(all(cpgi=="")) i=i+1
}
# TSS: 12409/16984 = 73%
# Body: 1650/16984 = 10%
# intergeneic: 2652/16984 = 16%
i=0
for (isl in ILS.hypo.drift1) {
  cpgi = manifestData[manifestData$Islands_Name==isl,"UCSC_RefGene_Group"]
  dum = manifestData[manifestData$Islands_Name==isl,"Relation_to_Island"]
  ind = grepl("Island",dum) | grepl("Shore",dum)
  cpgi = cpgi[ind]
  # if(sum(grepl("TSS",cpgi))>0) i=i+1
  # if(sum(grepl("Body",cpgi))>0 & !sum(grepl("TSS",cpgi))>0) i=i+1
  if(all(cpgi=="")) i=i+1
}
# TSS: 2552/4024 = 63%
# Body: 452/4024 = 11%
# intergenic: 947/4024 = 24%
