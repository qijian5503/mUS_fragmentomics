rm(list=ls())
gc()

load("medip_brca_samp135_phen.RData")
load("length_all_samp135.RData")
length_all_samp$Group=factor(length_all_samp$Group,levels = c("Normal","Tumor"))
ggplot(length_all_samp[length_all_samp$V2<250,], aes(V2, freq,fill=V3,color=Group))+
    geom_line(alpha=0.3)+xlab("cfDNA length")+ylab("Fraction")+theme_bw()+scale_color_manual(values=c("#00468BFF","#CA0020","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))

a=reshape2::dcast(length_all_samp,V2~Group,value.var = "V1",fill = 0,fun.aggregate = median)
a$ALL=a$Normal+a$Tumor
a$Normal=a$Normal/sum(a$Normal)
a$Tumor=a$Tumor/sum(a$Tumor)
a$ALL=a$ALL/sum(a$ALL)
a=reshape2::dcast(length_all_samp,V2~Group,value.var = "V1",fill = 0,fun.aggregate = sum)
a$Tumor=a$Tumor/sum(a$Tumor)
a$Normal=a$Normal/sum(a$Normal)

a=reshape2::melt(a,id.var="V2")
colnames(a)=c("cfDNA length","Group","Freq")
p1=ggplot(a[a$`cfDNA length`<250,], aes(`cfDNA length`, Freq,fill=`Group`,color=`Group`))+geom_line()+scale_color_manual(values=c("#00468BFF","#CA0020","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+ylab("Fraction")+theme_bw()
p1
p2=ggplot(a[a$`cfDNA length`<100 & !a$`Group` %in% "ALL",], aes(`cfDNA length`, Freq,fill=`Group`,color=`Group`))+geom_line()+scale_color_manual(values=c("#00468BFF","#CA0020","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+ylab("Fraction")+theme_bw()
p2

a=reshape2::dcast(length_all_samp,V2~V3,value.var = "freq",fill = 0,fun.aggregate = sum)
rownames(a)=a$V2
dat=data.frame(Group=phen_brca[match(colnames(a)[2:ncol(a)],phen_brca$Sample2),"Group"],
               prop=colSums(a[a$V2>=20 & a$V2<=85,2:ncol(a)]))
dat$Group=factor(dat$Group,levels = c("Normal","Tumor"))
library(ggplot2)
library(ggpubr)
p3=ggplot(dat, aes(Group, prop)) + 
    geom_boxplot(aes(fill=Group))+
    stat_compare_means(hide.ns = T)+
    xlab("fragment length below 85bp")+ ylab("Fraction")+geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+theme(legend.position = "")
p3


Normal_BRCA=roc(as.factor(dat$Group), as.numeric(dat$prop))
g=pROC::ggroc(list(`Normal VS Tumor`=Normal_BRCA))
p4=g+ggsci::scale_color_lancet()+geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+theme_classic() +
    annotate("text", y = 0.4, x=0.5, label = paste0("AUC = ",round(Normal_BRCA$auc,4)),colour ="#00468BFF")+
    theme(legend.position = "")
p4


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
table(imp_de_gene$Hypomethylation>90)
de_gene=as.character(imp_de_gene[imp_de_gene$Hypomethylation>90,"Var1"])
de_gene=de_gene[!grepl("chrX|chrY",de_gene)]


load("medip_brca_samp135_bed_20k_le85_ratio.RData")
dat=reshape2::dcast(all_dat2,pos~Sample,value.var = "ratio")
nam=str_split(dat$pos,pattern = "_",simplify = T)
nam=as.data.frame(nam)
nam$V2=ifelse(grepl("e[+]",nam$V2),format(as.numeric(nam$V2), scientific = F),nam$V2)
nam$V3=ifelse(grepl("e[+]",nam$V3),format(as.numeric(nam$V3), scientific = F),nam$V3)
nam$V4=ifelse(grepl("e[+]",nam$V4),format(as.numeric(nam$V4), scientific = F),nam$V4)
nam$V5=ifelse(grepl("e[+]",nam$V5),format(as.numeric(nam$V5), scientific = F),nam$V5)
dat$pos=ifelse(nam$V4 %in% "",paste(nam$V1,nam$V2,nam$V3,sep = "_"),
               ifelse(nam$V5 %in% "",paste(nam$V1,nam$V2,nam$V3,nam$V4,sep = "_"),paste(nam$V1,nam$V2,nam$V3,nam$V4,nam$V5,sep = "_")))
rownames(dat)=dat$pos
dat$pos=NULL

load("medip_brca_samp135_phen.RData")
colnames(dat)=phen_brca[match(colnames(dat),phen_brca$Sample2),"Sample"]
dat[is.na(dat)]=0
dat=as.data.frame(t(dat))
dat$Group=phen_brca[match(rownames(dat),phen_brca$Sample),"Group"]

p5=ggplot(dat, aes(Group, chr16_33860000_33880000)) + 
    geom_boxplot(aes(fill=Group))+
    stat_compare_means(hide.ns = T)+
    xlab("chr16_33860000_33880000")+ ylab("Short fragment ratio")+geom_jitter(width = 0.3, alpha = 0.7, color = 'black')+
    scale_fill_manual(values = c("#0571B0","#CA0020","#00468BFF","#42B540FF","#0099B4FF","#925E9FFF","#FDAF91FF"))+theme(legend.position = "")
p5

Normal_BRCA2=roc(as.factor(dat$Group), as.numeric(dat$chr16_33860000_33880000))
g2=pROC::ggroc(list(`Normal VS Tumor`=Normal_BRCA2))
p6=g2+ggsci::scale_color_lancet()+geom_segment(aes(x=0,xend=1,y=1,yend=0),color="darkgrey",linetype=4)+theme_classic() +
    annotate("text", y = 0.4, x=0.5, label = paste0("AUC = ",round(Normal_BRCA2$auc,4)),colour ="#00468BFF")+
    theme(legend.position = "")
p6


dat2=as.data.frame(t(dat[,de_gene]))
table(phen_brca$Sample %in% colnames(dat2))
colnames(phen_brca)
phen_brca=phen_brca[,c(5,6,2,4)]
phen_brca=phen_brca[order(phen_brca$Group),]

phen_brca$Group=factor(phen_brca$Group,levels = c("Normal","Tumor"))
phen_brca$Stage=factor(phen_brca$Stage,labels = c("Normal","Stage I","Stage II","Stage III"),levels = c("Normal","StageI","StageII","StageIII"))

median(phen_brca$Age,na.rm = T)
phen_brca$Age=ifelse(phen_brca$Age>=49,">=49",ifelse(phen_brca$Age<49,"<49",NA))
phen_brca$Age=factor(phen_brca$Age,levels = c("<49",">=49"))
colnames(phen_brca)
colnames(phen_brca)[4]="cfDNA Concentration"

ann_colors = list(
    Group=c(Tumor="#CA0020", Normal="#0571B0"),
    Age=c(`>=49`="#4DAF4A",`<49`="#984EA3"),
    Stage=c(`Normal`="#0571B0",`Stage I`="#F5F5F5",`Stage II`="#80CDC1",`Stage III`="#01665E")
)

p7=ComplexHeatmap::pheatmap(as.matrix(dat2[,rownames(phen_brca)]), scale = "row",cluster_cols = F,
                            show_rownames = F,show_colnames = F,annotation_colors = ann_colors,
                            col = colorRampPalette(c("navy","white","firebrick3"))(10),
                            annotation_col = phen_brca,annotation_names_row = T,
                            annotation_names_col = T ,column_title = NULL,row_title = NULL)
p7





library(ggplot2)
library(magrittr)
library(tidyverse)
library(fs)
library(grid)
library(data.table)
library(cowplot)
library(devtools)
library(here)

a<- fread("a_cut85_le500.tsv")
bins5mb <- as_tibble(a)
bins5mb <- bins5mb %>% mutate(ratio.cor = short.cor/long.cor)
bins5mb <- bins5mb %>% group_by(id) %>%
    mutate(ratio.centered = scale(ratio.cor, scale=FALSE)[,1])


load("medip_brca_samp135_phen.RData")
table(phen_brca$Sample2 %in% bins5mb$id)
phen_brca$id=phen_brca$Sample2
bins5mb=bins5mb[bins5mb$id %in% phen_brca$id,]
fp2 <- inner_join(bins5mb, phen_brca, by="id")

fp2 <-  arrange(fp2, id, bin) %>%
    mutate(bin=factor(bin),arm=factor(arm, levels=unique(arm))) %>%
    mutate(dx=factor(Group, levels=c("Normal","Tumor")))

panel.labels <- fp2 %>%
    group_by(dx) %>%
    summarize(n=length(unique(id)),.groups="drop") %>%
    mutate(labels=paste0(c("Normal (n=","Tumor (n="),n, ")"),arm="1p") %>%
    mutate(x=rep(1, 2), y=rep(0.2, 2))


arm <- fp2 %>% group_by(arm) %>%
    summarize(n=n(), .groups="drop") %>%
    mutate(arm = as.character(arm))
small.arms <- setNames(c("", "", "12q", "", "","", "", "", "","", "", "", "","", ""),
                       c("10p", "12p", "12q", "16p", "16q","17p", "17q", "18p", "18q","19p", "19q", "20p", "20q","21q", "22q"))
arm.labels <- setNames(arm$arm, arm$arm)
arm.labels[names(small.arms)] <- small.arms

p8 <- fp2 %>% group_by(dx) %>%
    ggplot(aes(x = bin, y = ratio.centered, group = id)) +
    geom_line(size = 0.1, color="gray60") +
    labs(x = "",y = "Fragmentation profile\n", color = "") +
    facet_grid(dx~arm, space="free_x", scales="free_x",switch="x") +    
    theme_classic() +theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),strip.background=element_blank())
p8


load("medip_brca_samp135_phen.RData")
load("medip_brca_peak_samp135_subsamp.RData")
load("medip_brca_peak_samp135_valid_samp.RData")
load("dat_all_cut85_le500.RData")
dat=dat[!grepl("cov_",rownames(dat)),]
table(rownames(phen_brca) %in% colnames(dat))
colnames(dat)=phen_brca[match(colnames(dat),phen_brca$Sample2),"Sample"]

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

phen_brca$Group=factor(phen_brca$Group,levels = c("Normal","Tumor"))
ann_colors = list(Group=c(Tumor="#CA0020", Normal="#0571B0"))
phen_brca=phen_brca[,c(4,5)]
phen_brca=phen_brca[order(phen_brca$Group),]
phen_brca$Concentration=NULL

p9=ComplexHeatmap::pheatmap(as.matrix(dat[FSD_gene,rownames(phen_brca)]), scale = "row",cluster_cols = F,
                            show_rownames = F,show_colnames = F,annotation_colors = ann_colors,
                            col = colorRampPalette(c("navy","white","firebrick3"))(100),
                            annotation_col = phen_brca,annotation_names_row = T,
                            annotation_names_col = T ,column_title = NULL,row_title = NULL)
p9




