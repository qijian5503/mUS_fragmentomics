## peak data

library(data.table)
library(tidyverse)
dat <- data.table::fread("peaks_exp.bed", data.table = F)
table(dat$V1)
dat=dat[!grepl("chrUn|hap|random",dat$V1),]
dat=dat[!grepl("chrM|chrX|chrY",dat$V1),]
dat$peak=paste(dat$V1,dat$V2,dat$V3,sep = "_")
dat=reshape2::dcast(dat,peak~V9,fill = 0,fun.aggregate = median,value.var = "V8")
rownames(dat)=dat$peak
dat$peak=NULL

load("medip_brca_samp135_phen.RData")
table(phen$Sample2 %in% colnames(dat))
rownames(phen)=phen$Sample2
table(colnames(dat) %in% phen$Sample2)
peak=dat[,colnames(dat) %in% phen_brca$Sample2]
colnames(peak)=phen_brca[match(colnames(peak),phen_brca$Sample2),"Sample"]
save(peak,file="medip_brca_peak_samp135_exp.RData")

phen_brca=phen_brca[phen_brca$Sample2 %in% colnames(dat),]
rownames(phen_brca)=phen_brca$Sample
save(phen_brca,file="medip_brca_samp135_phen.RData")


library(caret)
load("medip_brca_samp135_phen.RData")
table(phen_brca$Group)
valid_sample=createDataPartition(phen_brca$Group, p = 0.2,times = 1)
valid_sample=phen_brca$Sample[valid_sample[[1]]]
save(valid_sample,file="medip_brca_peak_samp135_valid_samp.RData")

load("medip_brca_peak_samp135_valid_samp.RData")
train_test=phen_brca[!phen_brca$Sample %in% valid_sample,]
train_test_samp=train_test$Sample
subsample=createDataPartition(train_test$Group, p = 0.8,times = 100)
for (i in 1:100) {
    train_samp=train_test_samp[subsample[[i]]]
    test_samp=train_test_samp[!train_test_samp %in% train_samp]
    subsample[[i]]=list(train_samp=train_samp,test_samp=test_samp)
}
save(subsample,file="medip_brca_peak_samp135_subsamp.RData")


## peak selset
library(limma)
library(edgeR)

load("medip_brca_peak_samp135_exp.RData")
load("medip_brca_samp135_phen.RData")
load("medip_brca_peak_samp135_valid_samp.RData")
peak=peak[!grepl("chrX|chrY",rownames(peak)),]
peak=peak[!grepl("chrUn|hap|random",rownames(peak)),]

table(phen_brca$Group)
train_test=phen_brca[!phen_brca$Sample %in% valid_sample,]
train_test_samp=train_test$Sample

train_test=peak[,train_test_samp]
keep <- rowSums(train_test > 10) >= ncol(train_test)*0.1
table(keep)
save(keep,file = "medip_brca_peak_samp135_keep_10_01.RData")

peak_nam=rownames(train_test)[keep]
peak_nam=str_split(peak_nam,"_",simplify = T)
write.table(peak_nam,"peak_keep_name.bed",quote = F,row.names = F,col.names = F,sep = "\t")

train_test=train_test[keep,]
f75 <- rep_len(1,ncol(train_test))
for (j in seq_len(ncol(train_test))) f75[j] <- quantile(train_test[,j], probs=0.75)
lib.size=colSums(train_test)
f75=f75 / lib.size
refColumn <- colnames(train_test)[which.min(abs(f75-mean(f75)))]
refColumn
save(refColumn,file="medip_brca_peak_samp135_refColumn_10_01.RData")

valid=peak[keep,c(valid_sample,refColumn)]

## peak norm
d=DGEList(counts = train_test) %>% calcNormFactors(. ,method = "TMM" ,refColumn = refColumn) %>% voom()
norm=as.data.frame(d$E)
norm[1:3,1:3]
save(norm,file = "medip_brca_peak_samp135_traintest_norm_voom_10_01_noxy.RData")

d=DGEList(counts = valid) %>% calcNormFactors(. ,method = "TMM" ,refColumn = refColumn) %>% voom()
norm_v=as.data.frame(d$E)
table(norm_v[,refColumn]==norm[,refColumn])
norm_v[,refColumn]=NULL
save(norm_v,file = "medip_brca_peak_samp135_norm_voom_10_01_valid_noxy.RData")

ref_exp=data.frame(V1=rownames(valid),V2=valid[,refColumn])
save(ref_exp,file = "peak_ref_exp.RData")


## de analysis   
rm(list=ls())
gc()
options(stringsAsFactors = F)

load("medip_brca_peak_samp135_subsamp.RData")
load("medip_brca_peak_samp135_traintest_norm_voom_10_01_noxy.RData")
load("medip_brca_samp135_phen.RData")

phen_brca=na.omit(phen_brca[,c(2,3,5)])
phen_brca$Gender=as.numeric(factor(phen_brca$Gender,levels = c("F","M")))-1
phen_brca$Group=factor(phen_brca$Group,levels = c("Normal","Tumor"))

limma_de_all=list()
for(i in 1:100){
    train_samp=subsample[[i]]$train_samp
    train_samp=train_samp[train_samp %in% rownames(phen_brca)]
    coldata=phen_brca[train_samp,]
    psm <- matchit(Group~Age+Gender,data=coldata,method="nearest",distance = "logit",ratio = 2)
    coldata <- match.data(psm)
    design=model.matrix(~0+coldata$subclass+coldata$Group)
    colnames(design)=gsub("coldata[$]","",colnames(design))
    res=lmFit(norm[,rownames(coldata)],design = design,) %>% eBayes(.,trend = T) %>% topTable(., number =Inf,adjust.method = "BH",coef="GroupTumor") %>% .[order(.$adj.P.Val),]  %>% na.omit(.)
    limma_de_all[[i]]=res[abs(res$logFC)>1 & res$adj.P.Val<0.01,]
    message(nrow(limma_de_all[[i]]))
}
# save(limma_de_all,file="medip_brca_peak_samp135_norm_voom_10_01_limma_PSM_de_noxy_all.RData")


## de region
load("medip_brca_samp135_phen.RData")
load("medip_brca_peak_samp135_subsamp.RData")
load("medip_brca_peak_samp135_valid_samp.RData")
load("medip_brca_peak_samp135_norm_voom_10_01_valid_noxy.RData")
load("medip_brca_peak_samp135_traintest_norm_voom_10_01_noxy.RData")
methy_data=cbind(norm,norm_v)

load("medip_brca_peak_samp135_norm_voom_10_01_limma_PSM_ratio3_de_noxy_all.RData")
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
table(imp_de_gene$Hypermethylation>95)
de_gene=as.character(imp_de_gene[imp_de_gene$Hypomethylation>95 | imp_de_gene$Hypermethylation>95,"Var1"])
methy_gene=de_gene[!grepl("chrX|chrY",de_gene)]
methy_data=methy_data[methy_gene,]
# save(methy_data,file = "methy_data.RData")



library(glmnet)
methy_de_exp=as.data.frame(t(norm[methy_gene,]))
methy_de_exp$Group=phen_brca[rownames(methy_de_exp),"Group"]
cv_fit = cv.glmnet(as.matrix(methy_de_exp[,methy_gene]), methy_de_exp$Group, family = "binomial",type.measure = "class")
coef.min = coef(cv_fit, s = "lambda.min")
methy_gene2=coef.min@Dimnames[[1]][coef.min@i+1][-1]  ## 筛选得到的基因
# save(cv_fit,file = "methy_lasso_cv_fit.RData")

methy_de_exp=as.data.frame(t(methy_data[methy_gene2,]))
methy_de_exp$group=phen_brca[rownames(methy_de_exp),"Group"]
methy_de_exp$group=factor(methy_de_exp$group,levels = c("Normal","Tumor"))

limma_modlist=list()
for(i in 1:100){
    train_samp=subsample[[i]]$train_samp
    train_samp=train_samp[train_samp %in% rownames(methy_de_exp)]
    test_samp=subsample[[i]]$test_samp
    test_samp=test_samp[test_samp %in% rownames(methy_de_exp)]
    train_data=methy_de_exp[rownames(methy_de_exp) %in% train_samp,c(methy_gene2,"group")]
    test_data=methy_de_exp[rownames(methy_de_exp) %in% test_samp,c(methy_gene2,"group")]
    valid_data=methy_de_exp[rownames(methy_de_exp) %in% valid_sample,c(methy_gene2,"group")]
    
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
# save(limma_modlist,file="medip_brca_peak_samp135_norm_voom_10_01_rf_modlist.RData") 


