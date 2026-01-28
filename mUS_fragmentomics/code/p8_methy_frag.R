## DMR non-DMR frag
rm(list=ls())
load("medip_brca_count_samp135_norm_voom_10_01_limma_PSM_de_noxy_all.RData")
all_de2=data.frame()
for (r in 1:100) {
    res=limma_de_all[[r]]
    res=res[abs(res$logFC)> 0.5 & res$adj.P.Val<0.001,]
    res$regule=ifelse(res$logFC>0,"Hypermethylation","Hypomethylation")
    res$peak=rownames(res)
    if(nrow(res)>0){all_de2=rbind(all_de2,data.frame(round=r,res))}
}
rownames(all_de2)=NULL

imp_de_gene=as.data.frame(table(all_de2$peak,all_de2$regule))
imp_de_gene=reshape2::dcast(imp_de_gene,Var1~Var2,value.var = "Freq")
table(imp_de_gene$Hypomethylation>10)
table(imp_de_gene$Hypermethylation>10)
de_gene=as.character(imp_de_gene[imp_de_gene$Hypomethylation>10 | imp_de_gene$Hypermethylation>10,"Var1"])
de_gene=de_gene[!grepl("chrX|chrY",de_gene)]
no_sig_gene=rownames(limma_de_all[[r]])[!rownames(limma_de_all[[r]]) %in% de_gene]

peak=as.data.frame(str_split(de_gene,"_",simplify = T))
colnames(peak)=c("chr","start","end")
peak$start=as.numeric(peak$start)
peak$end=as.numeric(peak$end)
peak$start=peak$start-1000
peak$end=peak$end+1000
peak <- as.data.table(peak)
setkey(peak, chr, start, end)

no_sig_peak=as.data.frame(str_split(no_sig_gene,"_",simplify = T))
colnames(no_sig_peak)=c("chr","start","end")
no_sig_peak$start=as.numeric(no_sig_peak$start)
no_sig_peak$end=as.numeric(no_sig_peak$end)
no_sig_peak$start=no_sig_peak$start-1000
no_sig_peak$end=no_sig_peak$end+1000
no_sig_peak <- as.data.table(no_sig_peak)
setkey(no_sig_peak, chr, start, end)


bedfile=list.files(path = "E:/proj/MEDIP/bam_bed_defil/",pattern = "bed$",full.names = T)
bed <- fread(bedfile[1],data.table = F)
bed = bed[bed$`13_seq_len` <500,]
colnames(bed)[1:3]=c("chr", "start", "end")
bed <- as.data.table(bed)
setkey(bed, chr, start, end)
fragdt <- foverlaps(bed,peak, type="any")
fragdt=na.omit(fragdt)
fragdt$Group="hyp"

fragdt_no_sig <- foverlaps(bed,no_sig_peak, type="any")
fragdt_no_sig=na.omit(fragdt_no_sig)
fragdt_no_sig$Group="no_sig"
fragdt=rbind(fragdt,fragdt_no_sig)
table(fragdt$Group)


dat=as.data.frame(table(fragdt$Group,fragdt$`13_seq_len`))
dat$Var2=as.numeric(as.character(dat$Var2))
ggplot(dat, aes(Var2, log2(Freq),color=Var1))+geom_line()+theme_classic()

a= data.frame(Group=fragdt$Group,length=fragdt$`13_seq_len`) %>% dplyr::group_by(Group,length) %>% 
    dplyr::count() %>% group_by(Group) %>% mutate(Freq=n/sum(n)*100)
head(a)
sum(a[a$Group %in% "hyp","Freq"])
ggplot(a, aes(length, Freq,color=Group))+geom_line()+theme_classic()

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
    fragdt$Group="hyp"
    
    fragdt_no_sig <- foverlaps(bed,no_sig_peak, type="any")
    fragdt_no_sig=na.omit(fragdt_no_sig)
    fragdt_no_sig$Group="no_sig"
    fragdt=rbind(fragdt,fragdt_no_sig)
    fragdt$Sample=i
    all_data=rbind(all_data,fragdt[,c(1:3,8,15:17)])
}
all_data$Sample=gsub("_read_length_GC.bed","",basename(all_data$Sample))
save(all_data,file = "medip_brca_count_samp135_de_methy_frag.RData")

