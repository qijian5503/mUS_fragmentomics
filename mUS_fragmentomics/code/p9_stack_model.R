## stack
rm(list=ls())
gc()

load("medip_brca_peak_samp135_norm_voom_10_01_rf_modlist2.RData")
methy_modlist=limma_modlist
load("medip_brca_bed_20k_samp135_rf_modlist2.RData")
bed_20k_modlist=limma_modlist
load("medip_brca_frag_prop_samp135_rf_modlist2.RData")
frag_prop_modlist=limma_modlist
load("medip_brca_FSD_samp135_rf_modlist2.RData")
FSD_modlist=limma_modlist
load("medip_brca_FSR_samp135_rf_modlist2.RData")
FSR_modlist=limma_modlist
dim(FSD_modlist[[1]]$rf_mod$importance)

all_pre=data.frame()
all_imp=data.frame()
for (i in 1:100) {
    for (n in c("bed_20k","frag_prop","FSD","FSR","methy")) {
        modlist=get(paste0(n,"_modlist"))
        all_pre=rbind(all_pre,data.frame(feature=n,round=i,modlist[[i]]$limma_pre))
        all_imp=rbind(all_imp,data.frame(feature=n,round=i,feature2=rownames(modlist[[i]]$rf_mod$importance),modlist[[i]]$rf_mod$importance))
    }
}

unique(all_pre$feature)
all_pre$feature=factor(all_pre$feature,levels = c("bed_20k","frag_prop","FSD","FSR","methy"),labels = c("mUS_FR","mUS_FLF","mUS_FSD","mUS_FSR","Methy"))
dat=reshape2::dcast(all_pre,round+type+sample+group~feature,value.var = "pre")
anyNA(dat)
colnames(dat)

load("medip_brca_samp135_phen.RData")
load("medip_brca_peak_samp135_de_methy_frag.RData")
table(all_data$Sample %in% phen_brca$Sample2)
unique(all_data$Group)
all_data$Group=factor(all_data$Group,levels = c("hyp","no_sig"),labels = c("DMR","Non-DMR"))
all_data$Group2=phen_brca[match(all_data$Sample,phen_brca$Sample2),"Group"]
dat3=do.call(rbind,apply(unique(all_data[,c(6:8)]),1,function(x){
    d=all_data[all_data$Group %in% x[1] & all_data$Sample %in% x[2] & all_data$Group2 %in% x[3], ]
    data.frame(Sample=x[2],Group=x[1],Group2=x[3],Freq85=nrow(d[d$`13_seq_len`<=85,])/nrow(d))
}))
dat3=reshape2::dcast(dat3,Sample+Group2~Group,value.var ="Freq85")
table(dat3$Sample %in% dat$sample)
table(dat3$Sample %in% phen_brca$Sample)
dat3$Sample=phen_brca[match(dat3$Sample,phen_brca$Sample2),"Sample"]
dat$mUS_FR_DMR=dat3[match(dat$sample,dat3$Sample),"DMR"]
dat$mUS_FR_non_DMR=dat3[match(dat$sample,dat3$Sample),"Non-DMR"]
dat$mUS_FR_DMR=(dat$mUS_FR_DMR-min(dat$mUS_FR_DMR))/(max(dat$mUS_FR_DMR)-min(dat$mUS_FR_DMR))
dat$mUS_FR_non_DMR=(dat$mUS_FR_non_DMR-min(dat$mUS_FR_non_DMR))/(max(dat$mUS_FR_non_DMR)-min(dat$mUS_FR_non_DMR))
colnames(dat)
# save(dat,file="medip_brca_stack_samp135_dat.RData")

limma_modlist=list()
for(i in 1:100){
    dat2=dat[dat$round %in% i ,]
    rownames(dat2)=dat2$sample
    dat2$sample=NULL
    train_data=dat2[dat2$type %in% "train",c("mUS_FR","mUS_FLF","mUS_FSD","mUS_FSR","mUS_FR_DMR","mUS_FR_non_DMR","Methy","group")]
    test_data=dat2[dat2$type %in% "test",c("mUS_FR","mUS_FLF","mUS_FSD","mUS_FSR","mUS_FR_DMR","mUS_FR_non_DMR","Methy","group")]
    valid_data=dat2[dat2$type %in% "valid",c("mUS_FR","mUS_FLF","mUS_FSD","mUS_FSR","mUS_FR_DMR","mUS_FR_non_DMR","Methy","group")]
    
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
# save(limma_modlist,file="medip_brca_stack_samp135_rf_modlist4.RData")

