## FSD
rm(list=ls())
gc()

load("medip_brca_samp135_phen.RData")
load("medip_brca_peak_samp135_subsamp.RData")
load("medip_brca_peak_samp135_valid_samp.RData")
load("dat_all_cut85_le500.RData")
dat=dat[!grepl("cov_",rownames(dat)),]
table(rownames(phen_brca) %in% colnames(dat))
colnames(dat)=phen_brca[match(colnames(dat),phen_brca$Sample2),"Sample"]

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
    res=lmFit(dat[,rownames(coldata)],design = design,) %>% eBayes(.,trend = T) %>% topTable(., number =Inf,adjust.method = "BH",coef="GroupTumor") %>% .[order(.$adj.P.Val),]  %>% na.omit(.)
    limma_de_all[[i]]=res[res$adj.P.Val<0.05,]
    message(nrow(limma_de_all[[i]]))
}
# save(limma_de_all,file="medip_brca_FSD_samp135_limma_PSM_de_all.RData")


load("medip_brca_FSD_samp135_limma_PSM_de_all.RData")
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
FSD_gene=de_gene[!grepl("chrX|chrY",de_gene)]

library(glmnet)
FSD_de_exp=as.data.frame(t(dat[FSD_gene,!colnames(dat) %in% valid_sample]))
FSD_de_exp$Group=phen_brca[rownames(FSD_de_exp),"Group"]
cv_fit = cv.glmnet(as.matrix(FSD_de_exp[,FSD_gene]), FSD_de_exp$Group, family = "binomial",type.measure = "class")
coef.min = coef(cv_fit, s = "lambda.min")
FSD_gene2=coef.min@Dimnames[[1]][coef.min@i+1][-1]  ## 筛选得到的基因
# save(cv_fit,file = "FSD_lasso_cv_fit.RData")

FSD_de_exp=as.data.frame(t(dat[FSD_gene2,]))
FSD_de_exp$group=phen_brca[rownames(FSD_de_exp),"Group"]
FSD_de_exp$group=factor(FSD_de_exp$group,levels = c("Normal","Tumor"))

limma_modlist=list()
for(i in 1:100){
    train_samp=subsample[[i]]$train_samp
    train_samp=train_samp[train_samp %in% rownames(FSD_de_exp)]
    test_samp=subsample[[i]]$test_samp
    test_samp=test_samp[test_samp %in% rownames(FSD_de_exp)]
    train_data=FSD_de_exp[rownames(FSD_de_exp) %in% train_samp,c(FSD_gene2,"group")]
    test_data=FSD_de_exp[rownames(FSD_de_exp) %in% test_samp,c(FSD_gene2,"group")]
    valid_data=FSD_de_exp[rownames(FSD_de_exp) %in% valid_sample,c(FSD_gene2,"group")]
    
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
# save(limma_modlist,file="medip_brca_FSD_samp135_rf_modlist.RData") 


