rm(list=ls())
gc()

load("medip_brca_peak_samp135_norm_voom_10_01_limma_PSM_de_noxy_all.RData")
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
table(imp_de_gene$Hypomethylation>90)
table(imp_de_gene$Hypermethylation>90)
de_gene=as.character(imp_de_gene[imp_de_gene$Hypomethylation>90 | imp_de_gene$Hypermethylation>90,"Var1"])
de_gene=de_gene[!grepl("chrX|chrY",de_gene)]

imp_de_gene$Group=ifelse(imp_de_gene$Hypermethylation>0,"Hypermethylation",
                         ifelse(imp_de_gene$Hypomethylation>0,"Hypomethylation",NA))
imp_de_gene$Var1=as.character(imp_de_gene$Var1)

load("medip_brca_peak_samp135_peakAnno.RData")
de_gene=peak_ana[peak_ana$pos %in% de_gene,]
table(de_gene$pos %in% imp_de_gene$Var1)
de_gene$regule=imp_de_gene[match(de_gene$pos,imp_de_gene$Var1),"Group"]
colnames(de_gene)

load("medip_brca_samp135_phen.RData")
load("medip_brca_peak_samp135_norm_voom_10_01_valid_noxy.RData")
load("medip_brca_peak_samp135_traintest_norm_voom_10_01_noxy.RData")
dat=cbind(norm,norm_v)
dat2=as.data.frame(dat[unique(de_gene$pos),])
table(phen_brca$Sample %in% colnames(dat2))
colnames(phen_brca)
phen_brca$Group=factor(phen_brca$Group,levels = c("Normal","Tumor"))
ann_colors = list(Group=c(Tumor="#CA0020", Normal="#0571B0"))
phen_brca=phen_brca[,c(4,5)]
phen_brca=phen_brca[order(phen_brca$Group),]
phen_brca$Concentration=NULL

p1=ComplexHeatmap::pheatmap(as.matrix(dat2[,rownames(phen_brca)]), scale = "row",cluster_cols = F,
                            show_rownames = F,show_colnames = F,annotation_colors = ann_colors,
                            col = colorRampPalette(c("navy","white","firebrick3"))(100),
                            annotation_col = phen_brca,annotation_names_row = T,
                            annotation_names_col = T ,column_title = NULL,row_title = NULL)
p1



colnames(de_gene)
de_gene$type=ifelse(grepl("^LOC",de_gene$SYMBOL),"non-coding RNA",de_gene$type)
de_gene=unique(de_gene[,c(1,15,6,7:9,11)])
de_gene$pos=NULL
de_gene$enhancer_type=NULL
colnames(de_gene)=c("Regule","Gene_location","Gene_type","CpG_type","Repeat_Region")

library(reshape2)
de_gene$link=1
de_gene$CpG_type=NULL
de_gene=melt(de_gene,id.vars = "link")
variable <- summary(de_gene$variable)
de_gene$flow <- rep(1:variable[1], length(variable))

mycol <- c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462','#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8',
           '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#66C2A5',
           '#6181BD', '#F34800', '#64A10E', '#FF00FF', '#c7475b', '#049a0b', '#BEAED4',
           '#FDC086', '#FFFF99', '#386CB0', '#F0027F', '#4253ff', '#ff4308', '#D8D155','#64495D', '#7CC767')

library(ggalluvial)
p2=ggplot(de_gene, aes(x = variable, y = 1,stratum = value, alluvium = flow, fill = value)) +
    geom_stratum() + geom_flow(aes.flow = 'forward') + scale_fill_manual(values = mycol) + 
    geom_text(stat = 'stratum', infer.label = TRUE, size = 2.5) +
    labs(x = '', y = '') + theme(legend.position = 'none', panel.background = element_blank(),line = element_blank(), axis.text.y = element_blank())
p2


load("medip_brca_peak_samp135_norm_voom_10_01_rf_modlist.RData")
all_pre=data.frame()
for (i in 1:100) {
    all_pre=rbind(all_pre,data.frame(round=i,limma_modlist[[i]]$limma_pre))
}

a=all_pre[!all_pre$type %in% "valid",]
b=all_pre[all_pre$type %in% "valid",]
Normal_BRCA1=roc(as.factor(a$group), as.numeric(a$pre))
Normal_BRCA2=roc(as.factor(b$group), as.numeric(b$pre))

g=pROC::ggroc(list(Test_set=roc(as.factor(a$group), as.numeric(a$pre)),
                   Valid_set=roc(as.factor(b$group), as.numeric(b$pre))))
library(ggsci)
pal_lancet("lanonc")(9)
# "#00468BFF" "#ED0000FF" "#42B540FF" "#0099B4FF" "#925E9FFF" "#FDAF91FF" "#AD002AFF" "#ADB6B6FF" "#1B1919FF"
p3=g+ggsci::scale_color_lancet()+geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+theme_classic() +
    annotate("text", y = 0.6, x=0.5, label = "Methylation Model",colour ="#925E9FFF")+
    annotate("text", y = 0.5, x=0.5, label = paste0("        Test set AUC = ",round(Normal_BRCA1$auc,4)),colour ="#00468BFF")+
    annotate("text", y = 0.4, x=0.5, label = paste0("Validation set AUC = ",round(Normal_BRCA2$auc,4)),colour ="#ED0000FF")+
    theme(legend.position = "")
p3

load("methy_lasso_cv_fit.RData")
plot(cv_fit)
p4 = recordPlot() 

load("medip_brca_samp135_phen.RData")
load("medip_brca_peak_samp135_de_methy_frag.RData")
table(all_data$Sample %in% phen_brca$Sample2)
unique(all_data$Group)
all_data$Group=factor(all_data$Group,levels = c("hyp","no_sig"),labels = c("DMR","Non-DMR"))

all_data$Group2=phen_brca[match(all_data$Sample,phen_brca$Sample2),"Group"]
a= data.frame(Group=all_data$Group,Group2=all_data$Group2,length=all_data$`13_seq_len`) %>% dplyr::group_by(Group2,Group,length) %>% 
    dplyr::count() %>% group_by(Group,Group2) %>% mutate(Freq=n/sum(n)*100)
head(a)
sum(a[a$Group %in% "DMR","Freq"])
a$Group=paste0(a$Group,"_",a$Group2)
p5=ggplot(a[a$length<250,], aes(length, Freq,fill=Group,color=Group))+
    geom_line()+scale_color_manual(values=c("#00468BFF","#CA0020","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+ylab("Fraction")+xlab("cfDNA length")+theme_bw()+
    theme(legend.position = "inside")
p5

dat3=do.call(rbind,apply(unique(all_data[,c(6:8)]),1,function(x){
    d=all_data[all_data$Group %in% x[1] & all_data$Sample %in% x[2] & all_data$Group2 %in% x[3], ]
    data.frame(Sample=x[2],Group=x[1],Group2=x[3],Freq85=nrow(d[d$`13_seq_len`<=85,])/nrow(d))
}))

p6=ggplot(dat3, aes(x = Group, y = Freq85)) +
    geom_boxplot(aes(fill = Group), show.legend = FALSE, width = 0.6) + 
    scale_fill_manual(values = c("#00468BFF", "#CA0020")) + 
    geom_point() + stat_compare_means(label = "p.format", paired = TRUE,hide.ns = T)+
    geom_line(aes(group = Sample), color = '#FDAF91FF', lwd = 0.5, alpha = 1)+
    facet_wrap(~Group2,scales = "free_y",shrink = TRUE)+theme_bw()+xlab("")+ylab("Percentage of fragment length below 85 bp")
p6

c=dat3[dat3$Group %in% "DMR",]
d=dat3[!dat3$Group %in% "DMR",]
Normal_BRCA3=roc(as.factor(c$Group2), as.numeric(c$Freq85))
Normal_BRCA4=roc(as.factor(d$Group2), as.numeric(d$Freq85))

g2=pROC::ggroc(list(DMR=roc(as.factor(c$Group2), as.numeric(c$Freq85)),
                    `Non-DMR`=roc(as.factor(d$Group2), as.numeric(d$Freq85))))
library(ggsci)
pal_lancet("lanonc")(9)
# "#00468BFF" "#ED0000FF" "#42B540FF" "#0099B4FF" "#925E9FFF" "#FDAF91FF" "#AD002AFF" "#ADB6B6FF" "#1B1919FF"
p7=g2+ggsci::scale_color_lancet()+geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+theme_classic() +
    # annotate("text", y = 0.6, x=0.5, label = "Percentage of fragment length below 85 bp",colour ="#925E9FFF")+
    annotate("text", y = 0.5, x=0.5, label = paste0("Fragment in DMR AUC = ",round(Normal_BRCA3$auc,4)),colour ="#00468BFF")+
    annotate("text", y = 0.4, x=0.5, label = paste0("Fragment in Non-DMR AUC = ",round(Normal_BRCA4$auc,4)),colour ="#ED0000FF")+
    theme(legend.position = "")
p7



