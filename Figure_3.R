#### Figure_3 (heatmap)
####
#### library(gplots)
#### library(grid)
# 
#### input: CpGs.hypo.drift5.isl; Islands.hypo.drift5; ILS.hypo.drift5
#### BETRNet methylation data: sqM.52, beM.64, lgdM, hgdM, eacM, and eacM.TCGA (TCGA)
#### patient info (age, sex): pheno (phenodata) 

cpgs.drift = CpGs.hypo.drift5.isl

# use first 10 NS samples (6th is missing age information)
dat = ilogit2(sqM.52[cpgs.drift,c(1:5,7:11)])
out.sq = CpGisl.level.info(set = ILS.hypo.drift5, isls = Islands.hypo.drift5,
  cpgs = CpGs.hypo.drift5.isl, dat=dat, full=T, prom=F)
# use all 64 BE samples 
dat = ilogit2(beM.64[cpgs.drift,])
out.beB64 = CpGisl.level.info(set = ILS.hypo.drift5, isls = Islands.hypo.drift5,
  cpgs = CpGs.hypo.drift5.isl, dat=dat, full=T, prom=F)

tmp = as.character(pheno$Patient_ID[match(colnames(beM.64),pheno$X)])
# truncate long IDs
tmp = substr(tmp,1,7)

# missing age and sex for 6th subject
sex = factor(c(sex.SQ52[c(1:5,7:11)],sex.BE64))
age = c(age.SQ52[c(1:5,7:11)],age.BE64)
age = factor(10*floor(age/10))

id.MM.BE = c(3, 18, 23, 25, 26, 30, 40, 42, 44, 45, 46, 59, 60, 63)
id.HM.BE = c(7, 10, 11, 15, 16, 19, 20, 21, 22, 28, 29, 35, 37, 38, 39, 51, 52, 53, 54, 55, 57, 61)
id.BM.BE = (1:64)[-c(id.MM.BE,id.HM.BE)]

pat = apply(out.beB64$Mvals,2,mean)
isl = apply(out.beB64$Mvals,1,mean)
# isl = out$bslope

pat1 = apply(out.beB64$Mvals[,id.MM.BE],2,mean)
pat2 = apply(out.beB64$Mvals[,c(id.BM.BE,id.HM.BE)],2,mean)
# pat3 = apply(out.beB64$Mvals[,id.HM.BE],2,mean)

type.order = c(1:10,10+order(pat1),10+length(pat1)+order(pat2)) 
# ,10+length(pat1)+length(pat2)+order(pat3))

aux = cbind(out.sq$Mvals,out.beB64$Mvals[,c(id.MM.BE,id.BM.BE,id.HM.BE)]) #ilogit2(out.be64$Mvals)
aux.o = aux[order(isl,decreasing = T),type.order] #c(1:10,10+order(pat))]

labels = c(colnames(sqM)[1:10],tmp[c(id.MM.BE,id.BM.BE,id.HM.BE)])
# labels.o = labels[type.order]

annotation_col = data.frame(Age=age[type.order],Sex=sex[type.order])
rownames(annotation_col) = labels.o = paste0("Subject", type.order)
tmp = RColorBrewer::brewer.pal(7, "Greys")
ann_colors =  list(Sex=c(Female="pink",Male="cyan"), Age=c("20"=tmp[1],"30"=tmp[2],"40"=tmp[3],
                                "50"=tmp[4],"60"=tmp[5],"70"=tmp[6],"80"=tmp[7]))
col.pal <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
colnames(aux.o) = labels.o

pheatmap((aux.o), cluster_cols=F, cluster_rows=F, 
        main = " ",
        color=col.pal,
        show_colnames = T, show_rownames = F,
        annotation_names_col = F,
        annotation_col = annotation_col,
        annotation_colors = ann_colors,
        labels_col = labels.o,
        fontsize_col = 6,
        gaps_col = c(10,24), 
         )

grid.text("SQ",x=0.07,y=.97)
grid.text("BE (unimodal)",x=.21,y=.97)
grid.text("BE (bimodal)",x=.5,y=.97)
grid.text("beta-value",x=.845,y=.36,rot = 90)
