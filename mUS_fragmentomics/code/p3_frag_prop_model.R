## fragment length
rm(list=ls())
gc()

find_peaks <- function (x, m = 3){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
        z <- i - m + 1
        z <- ifelse(z > 0, z, 1)
        w <- i + m + 1
        w <- ifelse(w < length(x), w, length(x))
        if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
    })
    pks <- unlist(pks)
    pks
}

load("medip_brca_samp135_phen.RData")
length_all_samp <- read.table("length_all_samp.txt", quote="\"", comment.char="")
table(phen_brca$Sample2 %in% length_all_samp$V3)
phen_brca$Sample2[!phen_brca$Sample2 %in% length_all_samp$V3]

length_all_samp$Group=phen_brca[match(length_all_samp$V3,phen_brca$Sample2),"Group"]
length_all_samp=na.omit(length_all_samp)
length_all_samp$a=paste0(length_all_samp$Group,length_all_samp$V3)

a=reshape2::dcast(length_all_samp,V3+Group~V2,value.var = "V1",fill = 0)
a$sum=rowSums(a[,3:ncol(a)])
length_all_samp$sum=a[match(length_all_samp$V3,a$V3),"sum"]
length_all_samp$freq=length_all_samp$V1/length_all_samp$sum
length_all_samp$Group=factor(length_all_samp$Group,levels = c("Normal","Tumor"))
save(length_all_samp,file = "length_all_samp135.RData")


## fragment prop
load("medip_brca_samp135_phen.RData")
load("length_all_samp135.RData")
table(length_all_samp$V3 %in% phen_brca$Sample2)
length_all_samp=length_all_samp[length_all_samp$V3 %in% phen_brca$Sample2,]
a=reshape2::dcast(length_all_samp,V2~Group,value.var = "V1",fill = 0,fun.aggregate = median)
a$ALL=a$Normal+a$Tumor
a$Normal=a$Normal/sum(a$Normal)
a$Tumor=a$Tumor/sum(a$Tumor)
a$ALL=a$ALL/sum(a$ALL)

peak=unique(c(20,find_peaks(a$ALL)+19,153,find_peaks(-a$ALL)+19))
peak=peak[peak<500]

head(length_all_samp)
a=reshape2::dcast(length_all_samp,V2~V3,value.var = "freq",fill = 0,fun.aggregate = sum)
rownames(a)=a$V2

all_exp=data.frame()
for (s in peak) {
    for (e in peak) {
        if(s+5 < e){
            prop=colSums(a[a$V2>=s & a$V2<e,])
            dat=data.frame(sample=names(prop),prop=as.numeric(prop),length=paste(s,e,sep = "_"))
            all_exp=rbind(all_exp,dat)
        }
    }
}
all_exp=reshape2::dcast(all_exp,length~sample,fill = 0,value.var = "prop",fun.aggregate = sum)
rownames(all_exp)=paste0("length_",all_exp$length)
all_exp$length=NULL
all_exp=all_exp[,colnames(all_exp) %in% phen_brca$Sample2]
colnames(all_exp)=phen_brca[match(colnames(all_exp),phen_brca$Sample2),"Sample"]
# save(all_exp,file = "frag_prop_all_exp_samp135.RData")


## de feature
table(colnames(all_exp) %in% phen_brca$Sample)
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
    res=lmFit(all_exp[,rownames(coldata)],design = design,) %>% eBayes(.,trend = T) %>% topTable(., number =Inf,adjust.method = "BH",coef="GroupTumor") %>% .[order(.$adj.P.Val),]  %>% na.omit(.)
    limma_de_all[[i]]=res[res$adj.P.Val<0.05,]
    message(nrow(limma_de_all[[i]]))
}
# save(limma_de_all,file="medip_brca_frag_prop_samp135_limma_PSM_ratio3_de_all.RData")


load("medip_brca_frag_prop_samp135_limma_PSM_ratio3_de_all.RData")
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
frag_prop_gene=de_gene[!grepl("chrX|chrY",de_gene)]

library(glmnet)
frag_prop_de_exp=as.data.frame(t(all_exp[frag_prop_gene,!colnames(all_exp) %in% valid_sample]))
frag_prop_de_exp$Group=phen_brca[rownames(frag_prop_de_exp),"Group"]
cv_fit = cv.glmnet(as.matrix(frag_prop_de_exp[,frag_prop_gene]), frag_prop_de_exp$Group, family = "binomial",type.measure = "class")
coef.min = coef(cv_fit, s = "lambda.min")
frag_prop_gene2=coef.min@Dimnames[[1]][coef.min@i+1][-1]  ## 筛选得到的基因
# save(cv_fit,file = "frag_prop_lasso_cv_fit.RData")

frag_prop_de_exp=as.data.frame(t(all_exp[frag_prop_gene2,]))
frag_prop_de_exp$group=phen_brca[rownames(frag_prop_de_exp),"Group"]
frag_prop_de_exp$group=factor(frag_prop_de_exp$group,levels = c("Normal","Tumor"))

limma_modlist=list()
for(i in 1:100){
    train_samp=subsample[[i]]$train_samp
    train_samp=train_samp[train_samp %in% rownames(frag_prop_de_exp)]
    test_samp=subsample[[i]]$test_samp
    test_samp=test_samp[test_samp %in% rownames(frag_prop_de_exp)]
    train_data=frag_prop_de_exp[rownames(frag_prop_de_exp) %in% train_samp,c(frag_prop_gene2,"group")]
    test_data=frag_prop_de_exp[rownames(frag_prop_de_exp) %in% test_samp,c(frag_prop_gene2,"group")]
    valid_data=frag_prop_de_exp[rownames(frag_prop_de_exp) %in% valid_sample,c(frag_prop_gene2,"group")]
    
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
# save(limma_modlist,file="medip_brca_frag_prop_samp135_rf_modlist.RData") 

