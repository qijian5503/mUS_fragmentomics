## bed 20k
rm(list=ls())
gc()

peak <- read.delim("hg19_20k_window.bed", header=FALSE)
colnames(peak)=c("chr","start","end")
peak$start=as.numeric(peak$start)
peak$end=as.numeric(peak$end)
peak <- as.data.table(peak)
setkey(peak, chr, start, end)

bedfile=list.files(path = "E:/proj/MEDIP/bam_bed_defil/",pattern = "bed$",full.names = T)
bed <- fread(bedfile[1],data.table = F)
bed = bed[bed$`13_seq_len` <500,]
colnames(bed)[1:3]=c("chr", "start", "end")
bed <- as.data.table(bed)
setkey(bed, chr, start, end)
fragdt <- foverlaps(bed,peak, type="any")
fragdt=na.omit(fragdt)


load("medip_brca_samp135_phen.RData")
bedfile=list.files(path = "E:/proj/MEDIP/bam_bed_defil/",pattern = "bed$",full.names = T)
table(basename(bedfile) %in% paste0(phen_brca$Sample2,"_read_length_GC.bed"))
bedfile=bedfile[basename(bedfile) %in% paste0(phen_brca$Sample2,"_read_length_GC.bed")]
all_data=data.frame()
for (i in bedfile) {
    bed <- fread(i,data.table = F)
    bed = bed[bed$`13_seq_len` <500,]
    colnames(bed)[1:3]=c("chr", "start", "end")
    bed <- as.data.table(bed)
    setkey(bed, chr, start, end)
    fragdt <- foverlaps(bed,peak, type="any")
    fragdt=na.omit(fragdt)
    fragdt$Sample=gsub("_read_length_GC.bed","",basename(i))
    all_data=rbind(all_data,fragdt[,c(1:3,8,15:16)])
}
# save(all_data,file = "medip_brca_samp135_bed_20k.RData")


load("medip_brca_samp135_bed_20k.RData")
all_data$pos=paste(all_data$chr,as.character(all_data$start),as.character(all_data$end),sep = "_")
all_data=all_data[,c("pos","13_seq_len","Sample")]

library(future.apply)
plan(multisession(workers = 4))
options(future.globals.maxSize= 89128960000)

all_dat2=data.frame()
p=unique(all_data$pos)
m=1
for( i in unique(all_data$Sample)){
    dat= as.data.frame(all_data[all_data$Sample %in% i,])
    dat2=do.call(rbind,future_lapply(p,function(x){
        len= na.omit(as.numeric(dat[dat$pos %in% x,"13_seq_len"]))
        if(!is_empty(len)){
            data.frame(pos=x,Sample=i,ratio=sum(len<85)/length(len))
        }else{data.frame(pos=x,Sample=i,ratio=NA)}
    }))
    all_dat2=rbind(all_dat2,dat2)
    message(m)
    m=m+1
    gc()
}
all_dat2$Group=phen_brca[match(all_dat2$Sample,phen_brca$Sample2),"Group"]
save(all_dat2,file = "medip_brca_samp135_bed_20k_le85_ratio.RData")


## bed 20k de

load("medip_brca_peak_samp135_subsamp.RData")
table(rownames(phen_brca) %in% colnames(all_dat2))
colnames(all_dat2)=phen_brca[match(colnames(all_dat2),phen_brca$Sample2),"Sample"]
all_dat2[is.na(all_dat2)]=0

load("medip_brca_samp135_phen.RData")
phen=na.omit(phen_brca[,c(2,3,5)])
phen$Gender=as.numeric(factor(phen$Gender,levels = c("F","M")))-1
phen$Group=factor(phen$Group,levels = c("Normal","Tumor"))

limma_de_all=list()
for(i in 1:100){
    train_samp=subsample[[i]]$train_samp
    train_samp=train_samp[train_samp %in% rownames(phen)]
    coldata=phen[train_samp,]
    psm <- matchit(Group~Age+Gender,data=coldata,method="nearest",distance = "logit",ratio = 3)
    coldata <- match.data(psm)
    design=model.matrix(~0+coldata$subclass+coldata$Group)
    colnames(design)=gsub("coldata[$]","",colnames(design))
    res=lmFit(all_dat2[,rownames(coldata)],design = design,) %>% eBayes(.,trend = T) %>% topTable(., number =Inf,adjust.method = "BH",coef="GroupTumor") %>% .[order(.$adj.P.Val),]  %>% na.omit(.)
    limma_de_all[[i]]=res[res$adj.P.Val<0.05,]
    message(nrow(limma_de_all[[i]]))
}
# save(limma_de_all,file="medip_brca_count_samp135_bed_20k_limma_PSM_de.RData")



## bed 20k
dat=reshape2::dcast(all_dat2,pos~Sample,value.var = "ratio")
table(rownames(phen_brca) %in% colnames(dat))
colnames(dat)=phen_brca[match(colnames(dat),phen_brca$Sample2),"Sample"]

load("medip_brca_count_samp135_bed_20k_limma_PSM_de.RData")
all_de2=data.frame()
for (r in 1:100) {
    res=limma_de_all[[r]]
    if(nrow(res)>100){res=res[1:100,]}
    res$regule=ifelse(res$logFC>0,"Hypermethylation","Hypomethylation")
    res$peak=rownames(res)
    if(nrow(res)>0){all_de2=rbind(all_de2,data.frame(round=r,res))}
}
rownames(all_de2)=NULL

imp_de_gene=as.data.frame(table(all_de2$peak,all_de2$regule))
imp_de_gene=reshape2::dcast(imp_de_gene,Var1~Var2,value.var = "Freq")
table(imp_de_gene$Hypomethylation>0,imp_de_gene$Hypermethylation>0 )
table(imp_de_gene$Hypomethylation>95)
de_gene=as.character(imp_de_gene[imp_de_gene$Hypomethylation>95,"Var1"])
bed_20k_gene=de_gene[!grepl("chrX|chrY",de_gene)]

bed_20k_dat=as.data.frame(t(dat[bed_20k_gene,]))
bed_20k_dat$group=phen_brca[rownames(bed_20k_dat),"Group"]
bed_20k_dat$group=factor(bed_20k_dat$group,levels = c("Normal","Tumor"))
# save(bed_20k_dat,file = "bed_20k_dat.RData")


library(glmnet)
bed_20k_de_exp=as.data.frame(t(dat[bed_20k_gene,!colnames(dat) %in% valid_sample]))
bed_20k_de_exp$Group=phen_brca[rownames(bed_20k_de_exp),"Group"]
cv_fit = cv.glmnet(as.matrix(bed_20k_de_exp[,bed_20k_gene]), bed_20k_de_exp$Group, family = "binomial",type.measure = "class")
coef.min = coef(cv_fit, s = "lambda.min")
bed_20k_gene2=coef.min@Dimnames[[1]][coef.min@i+1][-1] 
# save(cv_fit,file = "bed_20k_lasso_cv_fit.RData")

bed_20k_de_exp=as.data.frame(t(dat[bed_20k_gene2,]))
bed_20k_de_exp$group=phen_brca[rownames(bed_20k_de_exp),"Group"]
bed_20k_de_exp$group=factor(bed_20k_de_exp$group,levels = c("Normal","Tumor"))

limma_modlist=list()
for(i in 1:100){
    train_samp=subsample[[i]]$train_samp
    train_samp=train_samp[train_samp %in% rownames(bed_20k_de_exp)]
    test_samp=subsample[[i]]$test_samp
    test_samp=test_samp[test_samp %in% rownames(bed_20k_de_exp)]
    train_data=bed_20k_de_exp[rownames(bed_20k_de_exp) %in% train_samp,c(bed_20k_gene2,"group")]
    test_data=bed_20k_de_exp[rownames(bed_20k_de_exp) %in% test_samp,c(bed_20k_gene2,"group")]
    valid_data=bed_20k_de_exp[rownames(bed_20k_de_exp) %in% valid_sample,c(bed_20k_gene2,"group")]
    
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
# save(limma_modlist,file="medip_brca_bed_20k_samp135_rf_modlist.RData") 

