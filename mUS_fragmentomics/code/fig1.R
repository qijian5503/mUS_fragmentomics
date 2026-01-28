rm(list=ls())
options(stringsAsFactors = F)

load("medip_brca_samp135_phen.RData")
load("medip_brca_peak_samp135_valid_samp.RData")
table(phen_brca$Sample %in% valid_sample)
phen_brca$group2=ifelse(phen_brca$Sample %in% valid_sample,"Validation","Discovery")
phen_brca$Raw_Reads=round(phen_brca$Raw_Reads/1000000,2)
colnames(phen_brca)[7]=c("Reads Number(Millions)")
phen_brca$Stage=ifelse(phen_brca$Group %in% "Normal","Normal",phen_brca$Stage)

library(gtsummary)
t1=tbl_summary(phen_brca[phen_brca$group2 %in% "Discovery",c(2,4,5,6,7)],by = Group,
               statistic = list(all_continuous() ~ "{mean} ({sd})",all_categorical() ~ "{n} / {N} ({p}%)"),#确定分类变量格式
               digits = all_continuous() ~ 2,#确定小数点数
               missing_text = "(Missing)")
t2=tbl_summary(phen_brca[phen_brca$group2 %in% "Validation",c(2,4,5,6,7)],by = Group,
               statistic = list(all_continuous() ~ "{mean} ({sd})",all_categorical() ~ "{n} / {N} ({p}%)"),#确定分类变量格式
               digits = all_continuous() ~ 2,#确定小数点数
               missing_text = "(Missing)")
tbl_merge(tbls = list(t1, t2),tab_spanner = c("**Discovery Set**", "**Validation Set**"))



load("medip_brca_samp135_phen.RData")
table(phen_brca$Stage)
phen_brca$Stage=ifelse(phen_brca$Group %in% "Normal","Normal",phen_brca$Stage)

library(ggpubr)
phen_brca$Group=factor(phen_brca$Group,levels = c("Normal","Tumor"))
phen_brca$Stage=factor(phen_brca$Stage,levels = c("Normal","StageI","StageII","StageIII"),labels = c("Normal","Stage I","Stage II","Stage III"))
median(phen_brca$Age,na.rm = T)
phen_brca$Age=ifelse(phen_brca$Age>=49,">=49",ifelse(phen_brca$Age<49,"<49",NA))
phen_brca$Age=factor(phen_brca$Age,levels = c("<49",">=49"))

p2=ggplot(phen_brca[!is.na(phen_brca$Group) ,], aes(Group, Concentration)) +
    geom_boxplot(aes(fill = Group))+
    geom_jitter(width = 0.1, alpha = 0.5, color = 'black') + 
    stat_compare_means(aes(group=Group),hide.ns = T,label = "p")+labs(title = "")+xlab(label = "")+ylab("cfDNA Concentration(ng/ml)")+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+theme_bw()+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),legend.position = "")

p2

Normal_BRCA=roc(as.factor(phen_brca$Group), as.numeric(phen_brca$Concentration))
g=pROC::ggroc(list(`Normal VS Tumor`=Normal_BRCA))
p3=g+ggsci::scale_color_lancet()+geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+theme_classic() +
    annotate("text", y = 0.4, x=0.5, label = paste0("cfDNA Concentration\nAUC = ",round(Normal_BRCA$auc,4)),colour ="#00468BFF")+
    theme(legend.position = "")
p3




table(phen_brca$Stage)
my_comparisons=list(c("Stage II","Stage I"),c("Stage III","Stage I"),c("Stage III","Stage II"),c("Stage I","Normal"),c("Stage II","Normal"),c("Stage III","Normal"))
p5=ggplot(phen_brca[!is.na(phen_brca$Stage) ,], aes(Stage, Concentration)) +
    geom_boxplot(aes(fill = Stage))+
    geom_jitter(width = 0.1, alpha = 0.5, color = 'black') + 
    stat_compare_means(comparisons = my_comparisons,aes(group=Stage),hide.ns = T)+labs(title = "")+xlab(label = "")+ylab("cfDNA Concentration(ng/ml)")+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+theme_bw()+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),legend.position = "")
p5





library(gtsummary)
load("medip_brca_samp135_phen.RData")
phen_brca$Raw_Reads=round(phen_brca$Raw_Reads/1000000,2)
colnames(phen_brca)[7]=c("Reads Number(Millions)")
phen_brca$Stage=ifelse(phen_brca$Group %in% "Normal","Normal",phen_brca$Stage)

colnames(phen_brca)
t1=tbl_summary(phen_brca[,c(2:7)],by = Group,
               statistic = list(all_continuous() ~ "{mean} ({sd})",all_categorical() ~ "{n} / {N} ({p}%)"),#确定分类变量格式
               digits = all_continuous() ~ 2,#确定小数点数
               missing_text = "(Missing)")
t1

