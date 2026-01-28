## DEFIL data
rm(list=ls())
gc()

library(devtools)
load_all("PlasmaToolsHiseq.hg19-main/")
targetfile <- target55
target <- as.data.table(targetfile)
setnames(target, c("seqnames", "gcmed"), c("chr", "target"))

load_all("reproduce_lucas_wflow-master/code/PlasmaTools.lucas/")
load("reproduce_lucas_wflow-master/code/PlasmaTools.lucas/data/filters.hg19.rda")
filters <- as.data.table(filters.hg19)
setnames(filters, "seqnames",  "chr")
setkey(filters, chr, start, end)

bins <- data.table::fread("reproduce_lucas_wflow-master/code/preprocessing/bins_5mb.csv")
bins <- bins[map >= 0.90 & gc >= 0.30]
bins <- bins[,chr:=factor(chr, paste0("chr", c(1:22, "X")))]
bins[,bin:=1:.N]
setkey(bins, chr, start, end)


bedfile=list.files(path = "E:/proj/MEDIP/bam_bed_defil/",pattern = "bed$",full.names = T)
load("medip_brca_samp135_phen.RData")
table(basename(bedfile) %in% paste0(phen_brca$Sample2,"_read_length_GC.bed"))
bedfile=bedfile[basename(bedfile) %in% paste0(phen_brca$Sample2,"_read_length_GC.bed")]
binFrags <- function(fragments, bins, cutoff=85,
                     chromosomes=paste0("chr",c(1:22, "X"))) {
    setkey(bins, chr, start, end)
    fragbins <- foverlaps(fragments[chr %in% chromosomes],
                          bins, type="within", nomatch=NULL)
    bins2 <- fragbins[,.(arm=unique(arm), gc=gc[1], map=map[1],
                         short = sum(w >= 20 & w <= cutoff ),
                         long = sum(w > cutoff & w <= 500),
                         short.cor = sum(weight[w >= 20 & w <= cutoff]),
                         long.cor = sum(weight[w > cutoff & w <= 500]),
                         ultrashort = sum(w < 100),
                         ultrashort.cor = sum(weight[w < 100]),
                         multinucs = sum(w > 250),
                         multinucs.cor = sum(weight[w > 250]),
                         mediansize = as.integer(median(w)),
                         frag.gc = mean(fraggc)),
                      by=.(chr, start, end)]
    
    setkey(bins2, chr, start, end)
    bins2 <- bins2[bins]
    bins2 <- bins2[is.na(i.gc), which(grepl("short|long|multi", colnames(bins2))):=0]
    bins2[,`:=`(gc=i.gc, map=i.map, arm=i.arm)]
    bins2[,which(grepl("^i.", colnames(bins2))):=NULL]
    bins2[, bin:=1:.N]
    setcolorder(bins2, c("chr", "start", "end", "bin"))
    bins2[]
}

a=data.frame()
for (x in bedfile) {
    bed <- fread(x,data.table = F)
    bed = bed[bed$`13_seq_len` <500,]
    bed=bed[,c(1:4,6,13)]
    colnames(bed)=c("chr", "start", "end", "mapq","gc","w")
    bed <- as.data.table(bed)
    bed <- bed[,start:=start+1][mapq >= 30]
    bed <- bed[which(bed$w<1000)]
    
    fragments <- makeGRangesFromDataFrame(bed, keep.extra.columns=TRUE)
    fragments <- keepSeqlevels(fragments, paste0("chr", c(1:22, "X")),pruning.mode="coarse")
    fragments <- as.data.table(fragments)
    setnames(fragments,  "seqnames", "chr")
    fragdt <- foverlaps(fragments, filters, type="any")
    fragdt <- fragdt[is.na(start)][,`:=`(start=NULL, end=NULL, width=NULL,strand=NULL, name=NULL, score=NULL)]
    setnames(fragdt, c("i.start", "i.end", "i.width",  "i.strand"),c("start", "end", "width", "strand"))
    fragdt <- gcCorrectTarget(fragments, target)
    fragdt <- fragdt[,chr:=factor(chr, paste0("chr", c(1:22, "X")))]
    setnames(fragdt, "gc", "fraggc")
    setkey(fragdt, chr, start, end)
    
    bins2 <- binFrags(fragdt, bins,cutoff = 85)
    bins2[,id:=gsub("_read_length_GC.bed","",basename(x))]
    setcolorder(bins2, c("id", "chr", "start", "end", "bin"))
    a=rbind(a,as.data.frame(bins2))
    message(x)
}
data.table::fwrite(a,file = "a_cut85_le500.tsv",quote = F)

a<- fread("a_cut85_le500.tsv")
x<-setDT(a %>% filter(chr != "chrX"))
setkey(x, id)
x[,cov:=short+long]


load("medip_brca_samp135_phen.RData")
load("medip_brca_peak_samp135_valid_samp.RData")
phen_brca=phen_brca[!phen_brca$Sample %in% valid_sample,]
refids=phen_brca[phen_brca$Group %in% "Normal","Sample2"]

load_all("reproduce_lucas_wflow-master/code/PlasmaTools.lucas/R")
binsforzscores<-countAndNormalize(x, measure="cov")
armmeansdt_all <- PlasmaTools.lucas:::getArmMeans(binsforzscores)
armmeansdt_ref <- armmeansdt_all %>% filter(id %in% refids)
armmeansdt_lucas <- armmeansdt_all %>% filter(!id %in% refids)
stats <- armmeansdt_ref %>% group_by(arm) %>% dplyr::summarize(Mean=mean(armmean))
stats1 <- armmeansdt_ref %>% group_by(arm)  %>% dplyr::summarize(STD=sd(armmean))
stats <- inner_join(stats,stats1)

zscores <- rep(NA, length(armmeansdt_lucas$armmean))
for (i in 1:length(armmeansdt_lucas$armmean)){
    m <- armmeansdt_lucas[i]$armmean
    a <- armmeansdt_lucas[i]$arm
    mu <- (stats %>% filter(arm==a))$Mean
    sigma <- (stats %>% filter(arm==a))$STD
    zscores[i] <- (m - mu)/sigma
}
armmeansdt_lucas[,zscore := zscores]

#reframe the data-table shape
armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
armmeansdt_lucas[,arm:=factor(arm, armlevels)]
armmeansdt_lucas[,armvar:=factor(paste0("zscore_", arm),  paste0("zscore_", armlevels))]
features.zscores_lucas <- dcast(armmeansdt_lucas, id  ~ armvar, value.var="zscore")
write_csv(features.zscores_lucas,"zscores_hiseq55_cut85_le500.csv")

a<- fread("a_cut85_le500.tsv")
x_f <- a %>% filter(chr != "chrX")
setkey(x_f, id)
x_f[,cov:=short+long]
x_f[,ratio.cor:=short.cor/long.cor]
x_f[,ratio.scaled:=scale(ratio.cor), by=id]
x_f[,ratiovar:=factor(paste0("ratio_", bin), paste0("ratio_", 1:.N)), by=id]

features.ratios <- dcast(x_f, id ~ ratiovar, value.var="ratio.scaled")

x_f[,cov.cor:=short.cor+long.cor]
x_f[,cov.scaled:=scale(cov.cor), by=id]
x_f[,covvar:=factor(paste0("cov_", bin), paste0("cov_", 1:.N)), by=id]

features.covs <- dcast(x_f, id ~ covvar, value.var="cov.scaled")
setkey(features.ratios, id)
setkey(features.covs, id)
features.full <- features.covs[features.ratios]
write_csv(features.full,"coverage_hiseq55_cut85_le500.csv")

zscores<- read_csv("zscores_hiseq55_cut85_le500.csv")
cov <- read_csv("coverage_hiseq55_cut85_le500.csv")
setDT(zscores)
setDT(cov)
setkey(cov, id)
setkey(zscores, id)
full_features <- zscores[cov]
write_csv(full_features,"allfeatures_cut85_le500.csv")

full_features=as.data.frame(full_features)
rownames(full_features)=full_features$id
full_features$id=NULL
dat=as.data.frame(t(full_features))
dat[is.na(dat)]=0
save(dat,file = "dat_all_cut85_le500.RData")



## DEFIL de 
rm(list=ls())
gc()
options(stringsAsFactors = F)

load("medip_brca_peak_samp135_subsamp.RData")
load("dat_all_cut80_le500.RData")
load("medip_brca_samp135_phen.RData")
table(colnames(dat) %in% phen_brca$Sample2)
colnames(dat)=phen_brca[match(colnames(dat),phen_brca$Sample2),"Sample"]

phen=na.omit(phen_brca[,c(2,3,5)])
phen$Gender=as.numeric(factor(phen$Gender,levels = c("F","M")))-1
phen$Group=factor(phen$Group,levels = c("Normal","Tumor"))

limma_de_all=list()
for(i in 1:100){
    train_samp=subsample[[i]]$train_samp
    train_samp=train_samp[train_samp %in% rownames(phen)]
    coldata=phen[train_samp,]
    psm <- matchit(Group~Age+Gender,data=coldata,method="nearest",distance = "logit",ratio = 1)
    coldata <- match.data(psm)
    design=model.matrix(~0+coldata$subclass+coldata$Group)
    colnames(design)=gsub("coldata[$]","",colnames(design))
    res=lmFit(dat[,rownames(coldata)],design = design,) %>% eBayes(.,trend = T) %>% topTable(., number =Inf,adjust.method = "BH",coef="GroupTumor") %>% .[order(.$adj.P.Val),]  %>% na.omit(.)
    limma_de_all[[i]]=res[abs(res$logFC)>0.5 & res$adj.P.Val<0.05,]
    message(nrow(limma_de_all[[i]]))
}
# save(limma_de_all,file="medip_brca_defil_samp135_norm_voom_10_01_limma_PSM_de_noxy_all.RData")


## de feature
rm(list=ls())
options(stringsAsFactors = F)

load("medip_brca_peak_samp135_subsamp.RData")
load("medip_brca_peak_samp135_valid_samp.RData")
load("dat_all_cut80_le500.RData")
load("medip_brca_samp135_phen.RData")
table(colnames(dat) %in% phen_brca$Sample2)
colnames(dat)=phen_brca[match(colnames(dat),phen_brca$Sample2),"Sample"]

load("medip_brca_defil_samp135_norm_voom_10_01_limma_PSM_de_noxy_all.RData")

all_de2=data.frame()
for (r in 1:100) {
    res=limma_de_all[[r]]
    if(nrow(res)>10){res=res[1:10,]}
    res$regule=ifelse(res$logFC>0,"Hypermethylation","Hypomethylation")
    res$peak=rownames(res)
    if(nrow(res)>0){all_de2=rbind(all_de2,data.frame(round=r,res))}
}
rownames(all_de2)=NULL

imp_de_gene=as.data.frame(table(all_de2$peak,all_de2$regule))
imp_de_gene=reshape2::dcast(imp_de_gene,Var1~Var2,value.var = "Freq")
table(imp_de_gene$Hypermethylation>10)
de_gene=as.character(imp_de_gene[imp_de_gene$Hypomethylation>10 | imp_de_gene$Hypermethylation>10,"Var1"])


## rf model 
phen_brca$Group=factor(phen_brca$Group,levels = c("Normal","Tumor"))
de_exp=as.data.frame(t(dat[de_gene,!colnames(dat) %in% valid_sample]))
de_exp=as.data.frame(t(dat[de_gene,colnames(dat) %in% phen_brca$Sample]))
de_exp$group=phen_brca[rownames(de_exp),"Group"]

valid_data=as.data.frame(t(dat[de_gene,valid_sample]))
valid_data$group=as.factor(phen_brca[rownames(valid_data),"Group"])


ctrl <- trainControl(method = "repeatedcv",number = 10,repeats = 5,summaryFunction = twoClassSummary,classProbs = TRUE)
limma_modlist=list()
for(i in 1:100){
    train_samp=subsample[[i]]$train_samp
    train_samp=train_samp[train_samp %in% rownames(de_exp)]
    test_samp=subsample[[i]]$test_samp
    test_samp=test_samp[test_samp %in% rownames(de_exp)]
    train_data=de_exp[rownames(de_exp) %in% train_samp,c(de_gene,"group")]
    test_data=de_exp[rownames(de_exp) %in% test_samp,c(de_gene,"group")]
    
    rf_mod=randomForest(group~. ,train_data,important=TRUE,proximity=TRUE,ntree=100,mtry=10)
    train_pre=predict(rf_mod,train_data,type="prob")%>%data.frame
    train_pre$group=as.factor(phen_brca[rownames(train_pre),"Group"])
    test_pre=predict(rf_mod,test_data,type="prob")%>%data.frame
    test_pre$group=as.factor(phen_brca[rownames(test_pre),"Group"])
    valid_pre=predict(rf_mod,valid_data,type="prob")%>%data.frame
    valid_pre$group=as.factor(phen_brca[rownames(valid_pre),"Group"])
    
    limma_pre=rbind(data.frame(type="train",sample=rownames(train_pre),group=train_pre$group,pre=train_pre$Tumor),
                    data.frame(type="test",sample=rownames(test_pre),group=test_pre$group,pre=test_pre$Tumor),
                    data.frame(type="valid",sample=rownames(valid_pre),group=valid_pre$group,pre=valid_pre$Tumor))
    limma_modlist[[i]]= list(limma_pre=limma_pre,rf_mod=rf_mod)
}
# save(limma_modlist,file="medip_brca_defil_samp135_norm_voom_10_01_rf_modlist.RData") # 0.9877 0.9802


