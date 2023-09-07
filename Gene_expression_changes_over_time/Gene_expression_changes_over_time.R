library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
good.shapes = c(1:25,33:127)

#Immune response genes
g<-c( "Tnf", "Il1b",  "Ccl2", "Ccl5",  "Ccl20", "Cxcl10", "Irf7", "Irf9", "Isg15", "Gbp2", "Gbp3", "Gbp5", "Mx1", "Mx2", "Tlr3", "Tlr7", "Ifnar1", "Ifnar2")
##"Ifg","Il-6","Ccl19" not available


##For OE WT: 
filenames<-c("WT-UnInf-OE-vs-WT-07dpi-OE", "WT-UnInf-OE-vs-WT-14dpi-OE", "WT-UnInf-OE-vs-WT-30dpi-OE")

d<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select( log2FoldChange , Gene.name)%>%
    filter(Gene.name%in%g)
  de$time<-str_sub(filenames[i], -8, -4)
  d<-bind_rows(d, de)
}
d$time2<-as.numeric(gsub("dpi", "",d$time))

p<-ggplot(d, aes(time2,log2FoldChange, color= Gene.name , shape=Gene.name , group=Gene.name))+geom_point(size=3)+geom_line()+theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="", breaks=c(7, 14, 30), labels=c("7dpi", "14dpi", "30dpi"))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.title = element_blank())+
  scale_shape_manual(values=good.shapes[1:18]) 

ggsave(filename ="/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/OE_WT_Immune_response_genes.tiff",p, 
       width =16, height =14, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")

##For OE KO: 
filenames<-c("IRF3-KO-UnInf-OE-vs-IRF3-KO-07dpi-OE", "IRF3-KO-UnInf-OE-vs-IRF3-KO-14dpi-OE","IRF3-KO-UnInf-OE-vs-IRF3-KO-30dpi-OE")

d<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select( log2FoldChange , Gene.name)%>%
    filter(Gene.name%in%g)
  de$time<-str_sub(filenames[i], -8, -4)
  d<-bind_rows(d, de)
}
d$time2<-as.numeric(gsub("dpi", "",d$time))

p<-ggplot(d, aes(time2,log2FoldChange, color= Gene.name , shape=Gene.name , group=Gene.name))+geom_point(size=3)+geom_line()+theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="", breaks=c(7, 14, 30), labels=c("7dpi", "14dpi", "30dpi"))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.title = element_blank())+
  scale_shape_manual(values=good.shapes[1:18]) 

ggsave(filename ="/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/OE_KO_Immune_response_genes.tiff",p, 
       width =16, height =14, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")


#For OB WT:
filenames<-c("WT-UnInf-OB-vs-WT-07dpi-OB",  "WT-UnInf-OB-vs-WT-30dpi-OB")

d<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select( log2FoldChange , Gene.name)%>%
    filter(Gene.name%in%g)
  de$time<-str_sub(filenames[i], -8, -4)
  d<-bind_rows(d, de)
}
d$time2<-as.numeric(gsub("dpi", "",d$time))

p<-ggplot(d, aes(time2,log2FoldChange, color= Gene.name , shape=Gene.name , group=Gene.name))+geom_point(size=3)+geom_line()+theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="", breaks=c(7,  30), labels=c("7dpi",  "30dpi"))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.title = element_blank())+
  scale_shape_manual(values=good.shapes[1:18]) 

ggsave(filename ="/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/OB_WT_Immune_response_genes.tiff",p, 
       width =16, height =14, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")

##For OB KO: 
filenames<-c("IRF3-KO-UnInf-OB-vs-IRF3-KO-07dpi-OB", "IRF3-KO-UnInf-OB-vs-IRF3-KO-30dpi-OB")

d<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select( log2FoldChange , Gene.name)%>%
    filter(Gene.name%in%g)
  de$time<-str_sub(filenames[i], -8, -4)
  d<-bind_rows(d, de)
}
d$time2<-as.numeric(gsub("dpi", "",d$time))

p<-ggplot(d, aes(time2,log2FoldChange, color= Gene.name , shape=Gene.name , group=Gene.name))+geom_point(size=3)+geom_line()+theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="", breaks=c(7,  30), labels=c("7dpi", "30dpi"))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.title = element_blank())+
  scale_shape_manual(values=good.shapes[1:18]) 

ggsave(filename ="/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/OB_KO_Immune_response_genes.tiff",p, 
       width =16, height =14, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")



####For OE only: olfactory receptor genes: ORs - top 5 upregulated and 5 downregulated ORs
filenames<-c("WT-UnInf-OE-vs-WT-07dpi-OE", "WT-UnInf-OE-vs-WT-14dpi-OE", "WT-UnInf-OE-vs-WT-30dpi-OE")
de <- read.csv("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/WT-UnInf-OE-vs-WT-07dpi-OE/Differential_expression_analysis_table.csv")
OR<-de$Gene.name[grepl("Olfr",de$Gene.name)]
OR<-OR[!grepl("-p", OR)] #remove pseudogenes


up5<-de%>%
  filter(Gene.name%in%OR)%>%
  filter(log2FoldChange>1 &padj<0.05)%>%
  filter(!(Gene.name %in%c("Olfr66", "Olfr447", "Olfr355", "Olfr372", "Olfr978")))%>%
  top_n(5, wt=log2FoldChange)

dw5<-de%>%
  filter(Gene.name%in%OR)%>%
  filter(log2FoldChange< -1&padj<0.3)%>%
  top_n(5, wt=abs(log2FoldChange))

d<-data.frame()
for (i in 1:length(filenames)){
  dx <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select( log2FoldChange ,padj, Gene.name)%>%
    filter(Gene.name%in%OR)
  dx$time<-str_sub(filenames[i], -8, -4)
  d<-bind_rows(d, dx)
}

d$time2<-as.numeric(gsub("dpi", "",d$time))

p<-ggplot(d%>%filter(Gene.name%in%c(up5$Gene.name, dw5$Gene.name)), aes(time2,log2FoldChange, color= Gene.name , shape=Gene.name , group=Gene.name))+geom_point(size=3)+geom_line()+theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="", breaks=c(7, 14,30), labels=c("7dpi", "14dpi", "30dpi"))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.title = element_blank())+
  scale_shape_manual(values=good.shapes[1:18]) 

ggsave(filename ="/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/OE_WT_OR_UP_DW.tiff",p, 
       width =16, height =12, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")

###
filenames<-c("IRF3-KO-UnInf-OE-vs-IRF3-KO-07dpi-OE", "IRF3-KO-UnInf-OE-vs-IRF3-KO-14dpi-OE","IRF3-KO-UnInf-OE-vs-IRF3-KO-30dpi-OE")
de <- read.csv("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/IRF3-KO-UnInf-OE-vs-IRF3-KO-07dpi-OE/Differential_expression_analysis_table.csv")
OR<-de$Gene.name[grepl("Olfr",de$Gene.name)]
OR<-OR[!grepl("-p", OR)] #remove pseudogenes

up5<-de%>%
  filter(Gene.name%in%OR)%>%
  filter(log2FoldChange>1 &padj<0.05)%>%
  top_n(5, wt=log2FoldChange)
dw5<-de%>%
  filter(Gene.name%in%OR)%>%
  filter(log2FoldChange< -1&padj<0.05)%>%
  top_n(5, wt=abs(log2FoldChange))

d<-data.frame()
for (i in 1:length(filenames)){
  dx <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select( log2FoldChange ,padj, Gene.name)%>%
    filter(Gene.name%in%OR)
  dx$time<-str_sub(filenames[i], -8, -4)
  d<-bind_rows(d, dx)
}

d$time2<-as.numeric(gsub("dpi", "",d$time))

p<-ggplot(d%>%filter(Gene.name%in%c(up5$Gene.name, dw5$Gene.name)), aes(time2,log2FoldChange, color= Gene.name , shape=Gene.name , group=Gene.name))+geom_point(size=3)+geom_line()+theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="", breaks=c(7, 14,30), labels=c("7dpi", "14dpi", "30dpi"))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.title = element_blank())+
  scale_shape_manual(values=good.shapes[1:18]) 

ggsave(filename ="/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/OE_KO_OR_UP_DW.tiff",p, 
       width =16, height =12, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")





###For OE only: Olfactory transduction genes: OMP, CNGA2, Ano2, adcy3, gnal, gng13, Scn9a, Adcy3, Gfy
g<-c("Omp", "Cnga2", "Ano2", "Adcy3", "Gnal", "Gng13", "Scn9a", "Gfy")

##For OE WT: 
filenames<-c("WT-UnInf-OE-vs-WT-07dpi-OE", "WT-UnInf-OE-vs-WT-14dpi-OE", "WT-UnInf-OE-vs-WT-30dpi-OE")

d<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select( log2FoldChange , Gene.name)%>%
    filter(Gene.name%in%g)
  de$time<-str_sub(filenames[i], -8, -4)
  d<-bind_rows(d, de)
}
d$time2<-as.numeric(gsub("dpi", "",d$time))

p<-ggplot(d, aes(time2,log2FoldChange, color= Gene.name , shape=Gene.name , group=Gene.name))+geom_point(size=3)+geom_line()+theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="", breaks=c(7, 14, 30), labels=c("7dpi", "14dpi", "30dpi"))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.title = element_blank())+
  scale_shape_manual(values=good.shapes[1:18]) 

ggsave(filename ="/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/OE_WT_Olfactory_transduction_genes.tiff",p, 
       width =16, height =12, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")

##For OE KO: 
filenames<-c("IRF3-KO-UnInf-OE-vs-IRF3-KO-07dpi-OE", "IRF3-KO-UnInf-OE-vs-IRF3-KO-14dpi-OE","IRF3-KO-UnInf-OE-vs-IRF3-KO-30dpi-OE")

d<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select( log2FoldChange , Gene.name)%>%
    filter(Gene.name%in%g)
  de$time<-str_sub(filenames[i], -8, -4)
  d<-bind_rows(d, de)
}
d$time2<-as.numeric(gsub("dpi", "",d$time))

p<-ggplot(d, aes(time2,log2FoldChange, color= Gene.name , shape=Gene.name , group=Gene.name))+geom_point(size=3)+geom_line()+theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="", breaks=c(7, 14, 30), labels=c("7dpi", "14dpi", "30dpi"))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.title = element_blank())+
  scale_shape_manual(values=good.shapes[1:18]) 

ggsave(filename ="/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/OE_KO_Olfactory_transduction_genes.tiff",p, 
       width =16, height =12, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")


##Supporting cells markers: ezr, krt8,

g<-c("Ezr", "Krt8")

##For OE WT: 
filenames<-c("WT-UnInf-OE-vs-WT-07dpi-OE", "WT-UnInf-OE-vs-WT-14dpi-OE", "WT-UnInf-OE-vs-WT-30dpi-OE")

d<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select( log2FoldChange , Gene.name)%>%
    filter(Gene.name%in%g)
  de$time<-str_sub(filenames[i], -8, -4)
  d<-bind_rows(d, de)
}
d$time2<-as.numeric(gsub("dpi", "",d$time))

p<-ggplot(d, aes(time2,log2FoldChange, color= Gene.name , shape=Gene.name , group=Gene.name))+geom_point(size=3)+geom_line()+theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="", breaks=c(7, 14, 30), labels=c("7dpi", "14dpi", "30dpi"))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.title = element_blank())+
  scale_shape_manual(values=good.shapes[1:18]) 

ggsave(filename ="/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/OE_WT_Supporting_cells_markers.tiff",p, 
       width =16, height =12, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")

##For OE KO: 
filenames<-c("IRF3-KO-UnInf-OE-vs-IRF3-KO-07dpi-OE", "IRF3-KO-UnInf-OE-vs-IRF3-KO-14dpi-OE","IRF3-KO-UnInf-OE-vs-IRF3-KO-30dpi-OE")

d<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select( log2FoldChange , Gene.name)%>%
    filter(Gene.name%in%g)
  de$time<-str_sub(filenames[i], -8, -4)
  d<-bind_rows(d, de)
}
d$time2<-as.numeric(gsub("dpi", "",d$time))

p<-ggplot(d, aes(time2,log2FoldChange, color= Gene.name , shape=Gene.name , group=Gene.name))+geom_point(size=3)+geom_line()+theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="", breaks=c(7, 14, 30), labels=c("7dpi", "14dpi", "30dpi"))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.title = element_blank())+
  scale_shape_manual(values=good.shapes[1:18]) 

ggsave(filename ="/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/OE_KO_Supporting_cells_markers.tiff",p, 
       width =16, height =12, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")


##Neuronal targeting markers: Nrp1, Npr2, kirrel2, kirrel3, Robo2, Sema7A

g<-c("Nrp1", "Nrp2","Kirrel2", "Kirrel3", "Robo2","Sema7a")

##For OE WT: 
filenames<-c("WT-UnInf-OE-vs-WT-07dpi-OE", "WT-UnInf-OE-vs-WT-14dpi-OE", "WT-UnInf-OE-vs-WT-30dpi-OE")

d<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select( log2FoldChange , Gene.name)%>%
    filter(Gene.name%in%g)
  de$time<-str_sub(filenames[i], -8, -4)
  d<-bind_rows(d, de)
}
d$time2<-as.numeric(gsub("dpi", "",d$time))

p<-ggplot(d, aes(time2,log2FoldChange, color= Gene.name , shape=Gene.name , group=Gene.name))+geom_point(size=3)+geom_line()+theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="", breaks=c(7, 14, 30), labels=c("7dpi", "14dpi", "30dpi"))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.title = element_blank())+
  scale_shape_manual(values=good.shapes[1:18]) 

ggsave(filename ="/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/OE_WT_Neuronal_targeting_markers.tiff",p, 
       width =16, height =12, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")

##For OE KO: 
filenames<-c("IRF3-KO-UnInf-OE-vs-IRF3-KO-07dpi-OE", "IRF3-KO-UnInf-OE-vs-IRF3-KO-14dpi-OE","IRF3-KO-UnInf-OE-vs-IRF3-KO-30dpi-OE")

d<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select( log2FoldChange , Gene.name)%>%
    filter(Gene.name%in%g)
  de$time<-str_sub(filenames[i], -8, -4)
  d<-bind_rows(d, de)
}
d$time2<-as.numeric(gsub("dpi", "",d$time))

p<-ggplot(d, aes(time2,log2FoldChange, color= Gene.name , shape=Gene.name , group=Gene.name))+geom_point(size=3)+geom_line()+theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="", breaks=c(7, 14, 30), labels=c("7dpi", "14dpi", "30dpi"))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.title = element_blank())+
  scale_shape_manual(values=good.shapes[1:18]) 

ggsave(filename ="/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/OE_KO_Neuronal_targeting_markers.tiff",p, 
       width =16, height =12, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")


###Stem cells/developmental markers: Sox2, Krt5, Ascl1 (Mash1), Gap43, TP63, Insr

g<-c("Sox2", "Krt5","Ascl1", "Gap43", "Trp63","Insr")

##For OE WT: 
filenames<-c("WT-UnInf-OE-vs-WT-07dpi-OE", "WT-UnInf-OE-vs-WT-14dpi-OE", "WT-UnInf-OE-vs-WT-30dpi-OE")

d<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select( log2FoldChange , Gene.name)%>%
    filter(Gene.name%in%g)
  de$time<-str_sub(filenames[i], -8, -4)
  d<-bind_rows(d, de)
}
d$time2<-as.numeric(gsub("dpi", "",d$time))

p<-ggplot(d, aes(time2,log2FoldChange, color= Gene.name , shape=Gene.name , group=Gene.name))+geom_point(size=3)+geom_line()+theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="", breaks=c(7, 14, 30), labels=c("7dpi", "14dpi", "30dpi"))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.title = element_blank())+
  scale_shape_manual(values=good.shapes[1:18]) 

ggsave(filename ="/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/OE_WT_Stem_cells.developmental_markers.tiff",p, 
       width =16, height =12, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")

##For OE KO: 
filenames<-c("IRF3-KO-UnInf-OE-vs-IRF3-KO-07dpi-OE", "IRF3-KO-UnInf-OE-vs-IRF3-KO-14dpi-OE","IRF3-KO-UnInf-OE-vs-IRF3-KO-30dpi-OE")

d<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select( log2FoldChange , Gene.name)%>%
    filter(Gene.name%in%g)
  de$time<-str_sub(filenames[i], -8, -4)
  d<-bind_rows(d, de)
}
d$time2<-as.numeric(gsub("dpi", "",d$time))

p<-ggplot(d, aes(time2,log2FoldChange, color= Gene.name , shape=Gene.name , group=Gene.name))+geom_point(size=3)+geom_line()+theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="", breaks=c(7, 14, 30), labels=c("7dpi", "14dpi", "30dpi"))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.title = element_blank())+
  scale_shape_manual(values=good.shapes[1:18]) 

ggsave(filename ="/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/OE_KO_Stem_cells.developmental_markers.tiff",p, 
       width =16, height =12, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")


###Zonal markers: nqo1, Ncam2 (ocam)

g<-c("Nqo1", "Ncam2")

##For OE WT: 
filenames<-c("WT-UnInf-OE-vs-WT-07dpi-OE", "WT-UnInf-OE-vs-WT-14dpi-OE", "WT-UnInf-OE-vs-WT-30dpi-OE")

d<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select( log2FoldChange , Gene.name)%>%
    filter(Gene.name%in%g)
  de$time<-str_sub(filenames[i], -8, -4)
  d<-bind_rows(d, de)
}
d$time2<-as.numeric(gsub("dpi", "",d$time))

p<-ggplot(d, aes(time2,log2FoldChange, color= Gene.name , shape=Gene.name , group=Gene.name))+geom_point(size=3)+geom_line()+theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="", breaks=c(7, 14, 30), labels=c("7dpi", "14dpi", "30dpi"))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.title = element_blank())+
  scale_shape_manual(values=good.shapes[1:18]) 

ggsave(filename ="/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/OE_WT_Zonal_markers.tiff",p, 
       width =16, height =12, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")

##For OE KO: 
filenames<-c("IRF3-KO-UnInf-OE-vs-IRF3-KO-07dpi-OE", "IRF3-KO-UnInf-OE-vs-IRF3-KO-14dpi-OE","IRF3-KO-UnInf-OE-vs-IRF3-KO-30dpi-OE")

d<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select( log2FoldChange , Gene.name)%>%
    filter(Gene.name%in%g)
  de$time<-str_sub(filenames[i], -8, -4)
  d<-bind_rows(d, de)
}
d$time2<-as.numeric(gsub("dpi", "",d$time))

p<-ggplot(d, aes(time2,log2FoldChange, color= Gene.name , shape=Gene.name , group=Gene.name))+geom_point(size=3)+geom_line()+theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name="", breaks=c(7, 14, 30), labels=c("7dpi", "14dpi", "30dpi"))+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16), legend.title = element_blank())+
  scale_shape_manual(values=good.shapes[1:18]) 

ggsave(filename ="/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/OE_KO_Zonal_markers.tiff",p, 
       width =16, height =12, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")
