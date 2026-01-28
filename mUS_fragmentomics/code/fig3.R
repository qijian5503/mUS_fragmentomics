rm(list=ls())
gc()

list.files(pattern = "cv_fit.RData")
# "bed_20k_lasso_cv_fit.RData"   "frag_prop_lasso_cv_fit.RData"
# "FSD_lasso_cv_fit.RData"       "FSR_lasso_cv_fit.RData"       

list.files(pattern = "modlist")
load("medip_brca_bed_20k_samp135_rf_modlist.RData")
bed_20k_modlist=limma_modlist
load("medip_brca_frag_prop_samp135_rf_modlist.RData")
frag_prop_modlist=limma_modlist
load("medip_brca_FSD_samp135_rf_modlist.RData")
FSD_modlist=limma_modlist
load("medip_brca_FSR_samp135_rf_modlist.RData")
FSR_modlist=limma_modlist
dim(FSD_modlist[[1]]$rf_mod$importance)


all_pre=data.frame()
all_imp=data.frame()
for (i in 1:100) {
    for (n in c("bed_20k","frag_prop","FSD","FSR")) {
        modlist=get(paste0(n,"_modlist"))
        all_pre=rbind(all_pre,data.frame(feature=n,round=i,modlist[[i]]$limma_pre))
        all_imp=rbind(all_imp,data.frame(feature=n,round=i,feature2=rownames(modlist[[i]]$rf_mod$importance),modlist[[i]]$rf_mod$importance))
    }
}

c("bed_20k","frag_prop","FSD","FSR")
# "#00468BFF" "#ED0000FF" "#42B540FF" "#0099B4FF" "#925E9FFF" "#FDAF91FF" "#AD002AFF" "#ADB6B6FF" "#1B1919FF"

unique(all_pre$feature)
all_pre$feature=factor(all_pre$feature,
                       levels = c("bed_20k","frag_prop","FSD","FSR"),
                       labels = c("mus_FR","FLF","mus_FSD","mus_FSR"))
dat=reshape2::dcast(all_pre,round+type+sample+group~feature,value.var = "pre")
dat[is.na(dat)]=0
colnames(dat)
# colnames(dat)[5:8]=c("mus_FR","FLF","mus_FSD","mus_FSR")

for (n in c("mus_FR","FLF","mus_FSD","mus_FSR")) {
    a=all_pre[all_pre$feature %in% n & all_pre$type %in% "test",]
    b=all_pre[all_pre$feature %in% n & all_pre$type %in% "valid",]
    a1=roc(a$group,as.numeric(a$pre))
    auc1=roc(as.numeric(a1$response), as.numeric(a1$predictor))$auc %>% round(4)
    a2=roc(b$group,as.numeric(b$pre))
    auc2=roc(as.numeric(a2$response), as.numeric(a2$predictor))$auc %>% round(4)
    
    p=pROC::ggroc(list(test=a1,valid=a2))+
        ggsci::scale_color_lancet()+geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+
        annotate("text", y = 0.6, x=0.5, label = n,colour ="#925E9FFF")+
        annotate("text", y = 0.5, x=0.5, label = paste0("Test set AUC   = ",auc1),colour ="#00468BFF")+
        annotate("text", y = 0.4, x=0.5, label = paste0("Validation set AUC   = ",auc2),colour ="#ED0000FF")+
        theme_classic() + theme(legend.position = "")
    assign(paste0(n,"_ggroc"),p)
}


load("medip_brca_samp135_phen.RData")
table(phen_brca$Sample %in% all_pre$sample)
all_pre2=merge(all_pre,phen_brca,by.x="sample",by.y="Sample")
unique(all_pre2$feature)

p5=ggplot(all_pre2, aes(Group, pre)) +
    geom_violin(aes(fill=Group),scale = "width",trim=T)+
    labs(title = "")+xlab(label = "")+ylab("Probability of Cancer")+
    facet_grid(.~feature,scales="free")+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+theme_bw()+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),legend.position = "")
p5


p6=ggplot(all_pre2[!is.na(all_pre2$Stage),],aes(Stage, pre))+
    geom_violin(aes(fill=Stage),scale = "width",trim=T)+
    facet_grid(.~feature,scales = "free")+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+
    scale_x_discrete("")+scale_y_continuous("")+theme_bw()+ylab("Probability of Cancer")+
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),legend.position = "")
p6


