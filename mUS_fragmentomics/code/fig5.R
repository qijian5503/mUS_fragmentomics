rm(list=ls())
gc()

load("medip_brca_stack_samp135_rf_modlist.RData")
all_pre2=data.frame()
for (i in 1:100) {
    all_pre2=rbind(all_pre2,data.frame(round=i,limma_modlist[[i]]$limma_pre))
}

a=all_pre2[all_pre2$type %in% "test",]
b=all_pre2[all_pre2$type %in% "valid",]
a1=roc(a$group,as.numeric(a$pre))
auc1=roc(as.numeric(a1$response), as.numeric(a1$predictor))$auc %>% round(4)
a2=roc(b$group,as.numeric(b$pre))
auc2=roc(as.numeric(a2$response), as.numeric(a2$predictor))$auc %>% round(4)

p1=pROC::ggroc(list(test=a1,valid=a2))+
    ggsci::scale_color_lancet()+geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+
    annotate("text", y = 0.6, x=0.5, label = "Stack Model",colour ="#925E9FFF")+
    annotate("text", y = 0.5, x=0.5, label = paste0("Test set AUC   = ",auc1),colour ="#00468BFF")+
    annotate("text", y = 0.4, x=0.5, label = paste0("Validation set AUC   = ",auc2),colour ="#ED0000FF")+
    theme_classic() + theme(legend.position = "")
p1

library(fastshap)
library(magrittr)
library(tidyverse)
load("medip_brca_stack_samp135_dat.RData")
X=dat[,c("mUS_FR","mUS_FLF","mUS_FSD","mUS_FSR","mUS_FR_DMR","mUS_FR_non_DMR","Methy","group")]
X$group=as.numeric(X$group)-1
rf_mod=randomForest(group~. ,X)

library(fastshap)
shap_values <- explain(rf_mod,X=X[,-8],nsim=10,pred_wrapper = function(model,newdata){predict(rf_mod, newdata = newdata, type = "class")})
shap_handle <- shap_values %>% as.data.frame() %>% mutate(id=1:n()) %>% pivot_longer(cols = -(ncol(X[,-8])+1),values_to="shap") # 长宽数据转换
shap_handle

data4 <- X %>% mutate(id=1:n()) %>% pivot_longer(cols = -(ncol(X[,-8])+1))

shap_scale <- shap_handle %>% dplyr::rename("feature"="name")%>%
    group_by(feature)%>%mutate(shap=(shap-min(shap))/(max(shap)-min(shap)))
shap_scale <- shap_handle %>% left_join(data4)%>%
    dplyr::rename("feature"="name")%>% group_by(feature)%>%
    mutate(value=(value-min(value))/(max(value)-min(value))) %>% sample_n(200)
shap_scale=as.data.frame(shap_scale)
do.call(rbind,lapply(unique(shap_scale$feature), function(x){
    data.frame(feature=x,shap_imp=mean(abs(shap_scale[shap_scale$feature %in% x,"shap"]),na.rm=T))
})) %>% .[order(.$shap_imp,decreasing = T),]

# 5 mUS_FR_non_DMR 0.125660000
# 4     mUS_FR_DMR 0.060033267
# 6        mUS_FSD 0.040705314
# 2        mUS_FLF 0.038673467
# 3         mUS_FR 0.017014383
# 1          Methy 0.012681605
# 7        mUS_FSR 0.009994095

shap_scale$feature=factor(shap_scale$feature,
                          levels = c("mUS_FR_non_DMR","mUS_FR_DMR","mUS_FSD","mUS_FLF","mUS_FR","Methy","mUS_FSR"),
                          labels = c("mus_FR_non_DMR","mus_FR_DMR","mus_FSD","FLF","mus_FR","Methylation","mus_FSR"))
shap_scale$feature=factor(shap_scale$feature,levels = rev(c("mus_FR_non_DMR","mus_FR_DMR","mus_FSD","FLF","mus_FR","Methylation","mus_FSR")))
p2=ggplot(data=shap_scale, aes(x=shap, y=feature, color=shap)) +
    geom_jitter(size=2, height=0.1, width=0) +
    scale_color_gradient(low="#FFCC33", high="#6600CC", name="SHAP value", aesthetics = c("color")) +
    theme_bw()+ylab("")+xlab("SHAP")
p2



load("medip_brca_samp135_phen.RData")
table(all_pre2$sample %in% phen_brca$Sample)
unique(phen_brca$Stage)
phen_brca$Stage=ifelse(is.na(phen_brca$Stage) & phen_brca$Group %in% "Normal","Normal",phen_brca$Stage)
phen_brca$Risk=lapply(phen_brca$Sample, function(x){median(all_pre2[all_pre2$sample %in% x,"pre"],na.rm = T)})
phen_brca$Risk=as.numeric(phen_brca$Risk)

phen_brca2=phen_brca
median(phen_brca2$Age,na.rm = T)
phen_brca2$Age=ifelse(phen_brca2$Age >= 60,">=60",ifelse(phen_brca2$Age < 60,"<60",NA))

median(phen_brca2$Concentration,na.rm = T)
phen_brca2$Concentration=ifelse(phen_brca2$Concentration >= 4.2125,">=4.2125",ifelse(phen_brca2$Concentration < 4.2125,"<4.2125",NA))

median(phen_brca2$Raw_Reads,na.rm = T)
phen_brca2$Raw_Reads=ifelse(phen_brca2$Raw_Reads >= 29210767,">=29.2M",ifelse(phen_brca2$Raw_Reads < 29210767,"<29.2M",NA))

phen_brca3=reshape2::melt(phen_brca2[,c(1:2,4:5,7,11)],id.vars = c("Sample","Group","Risk"))
phen_brca3$Feature=paste0(phen_brca3$variable," ",phen_brca3$value)

p3=ggplot(phen_brca3[!is.na(phen_brca3$value),], aes(Feature, Risk,fill=Feature)) + 
    geom_boxplot()+geom_jitter(width = 0.1, alpha = 0.5, color = 'black')+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+
    stat_compare_means(aes(group=Feature),label.y = 0.5,label = "p.format")+
    facet_grid(.~Group)+theme_bw()+ylab("Probability of Cancer")+xlab("")+
    theme(panel.grid = element_blank(),legend.position = "")+
    coord_flip()
p3

unique(phen_brca$Stage)
my_comparisons=list(c("StageI","StageII"),c("StageI","StageIII"),c("StageII","StageIII"))
p4=ggplot(phen_brca[!is.na(phen_brca$Stage),],aes(Stage, Risk)) + 
    geom_boxplot(aes(fill = Stage))+xlab(label = "")+
    stat_compare_means(comparisons = my_comparisons,aes(group=Stage),label = "p.format")+
    geom_point(aes(fill=Stage),position=position_jitterdodge(dodge.width=0.85), alpha=0.2)+theme_bw()+
    scale_fill_manual(values = c("#0571B0","#42B540FF","#CA0020","#0099B4FF","#925E9FFF","#FDAF91FF"))+ylab("Probability of Cancer")
p4


g1=pROC::ggroc(list(Age=roc(as.factor(phen_brca$Group), as.numeric(phen_brca$Age)),
                    Concentration=roc(as.factor(phen_brca$Group), as.numeric(phen_brca$Concentration)),
                    Risk=roc(as.factor(phen_brca$Group), as.numeric(phen_brca$Risk))))
auc1=roc(as.factor(phen_brca$Group), as.numeric(phen_brca$Age))$auc[[1]]
auc2=roc(as.factor(phen_brca$Group), as.numeric(phen_brca$Concentration))$auc[[1]]
auc3=roc(as.factor(phen_brca$Group), as.numeric(phen_brca$Risk))$auc[[1]]

library(ggsci)
pal_lancet("lanonc")(9)
# "#00468BFF" "#ED0000FF" "#42B540FF" "#0099B4FF" "#925E9FFF" "#FDAF91FF" "#AD002AFF" "#ADB6B6FF" "#1B1919FF"
p5=g1+ggsci::scale_color_lancet() +
    theme_classic() + theme(legend.title = element_blank())+
    geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+
    annotate("text", y = 0.6, x=0.5, label = paste0("Age AUC = ",round(auc1,4)),colour ="#00468BFF")+
    annotate("text", y = 0.5, x=0.5, label = paste0("Concentration AUC = ",round(auc2,4)),colour ="#ED0000FF")+
    annotate("text", y = 0.4, x=0.5, label = paste0("Stack Model AUC = ",round(auc3,4)),colour ="#42B540FF")+
    theme_bw()+theme(legend.position = "")
p5



