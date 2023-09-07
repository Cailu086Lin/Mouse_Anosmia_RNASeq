library(DESeq2)
library(tximeta)
library(pheatmap)
library(dplyr)
library(stringr)


library(biomaRt)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
listDatasets(ensembl)
ensembl_ids <- d$Geneid # Replace with your list of Ensembl Gene IDs.
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                   filters = "ensembl_gene_id", 
                   values = ensembl_ids,
                   mart = ensembl)

##extract row counts 
filename <- read.table("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/hit-counts/filenames.txt", quote="\"", comment.char="")[,1]
for (i in 1:length(filename)){
  d<-read.delim(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/hit-counts/", filename[i]))[,c(1,7)]
  write.table(d,paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/hit-counts/test/sorted.", filename[i]), quote=FALSE, row.names = FALSE, col.names = FALSE)
}

####top50 up_down genes
filenames<-c("WT-UnInf-OE-vs-WT-07dpi-OE", "WT-UnInf-OE-vs-WT-14dpi-OE", "WT-UnInf-OE-vs-WT-30dpi-OE")

dup<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select(ID, log2FoldChange,padj, Gene.name)%>%
    filter(padj<0.05&log2FoldChange>1)%>%
    top_n(20, wt=log2FoldChange)
  dup<-bind_rows(dup, de)
}

dwn<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select(ID, log2FoldChange,padj, Gene.name)%>%
    filter(padj<0.05&log2FoldChange< -1)%>%
    top_n(27, wt=abs(log2FoldChange))
  dwn<-bind_rows(dwn, de)
}


##htseq-count input

directory <-"/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/hit-counts/test/WTOE"

sampleFiles <- grep("txt",list.files(directory),value=TRUE)


sampleTable <- data.frame(sampleName =sampleFiles,
                          fileName = sampleFiles,
                          type =c("UnInf","UnInf","UnInf", "7dpi", "7dpi", "7dpi", "14dpi", "14dpi", "14dpi", "30dpi", "30dpi", "30dpi"),
                          group = c("UnInf","UnInf","UnInf","Inf","Inf","Inf","Inf","Inf","Inf","Inf","Inf","Inf"))

sampleTable$group<-factor(sampleTable$group)
sampleTable$type<-factor(sampleTable$type)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ group)
ddsHTSeq

##
dds<-collapseReplicates(ddsHTSeq, groupby=ddsHTSeq@colData@listData[["type"]], renameCols = TRUE)

##normalize

vsd <- vst(dds, blind=FALSE)

head(assay(vsd), 3)

##Interaction
dup<-dup%>%distinct(Gene.name, .keep_all = TRUE)
###
dwn<-dwn%>%distinct(Gene.name, .keep_all = TRUE)
select<-dwn$ID
mt<-data.frame(t(scale(t(assay(vsd))))[select,])[, c(4,3,1,2)]
rownames(mt)<-dwn$Gene.name
colnames(mt)<-c("control", "7dpi", "14dpi", "30dpi")
###

df<-data.frame(Group=c("control", "infected","infected","infected"))
rownames(df)<-colnames(mt)

pheatmap(mt, cluster_rows=FALSE, show_rownames=TRUE, 
         cluster_cols=FALSE, annotation_col=df)

###KO_OE
filenames<-c("IRF3-KO-UnInf-OE-vs-IRF3-KO-07dpi-OE", "IRF3-KO-UnInf-OE-vs-IRF3-KO-14dpi-OE","IRF3-KO-UnInf-OE-vs-IRF3-KO-30dpi-OE")
dup<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select(ID, log2FoldChange,padj, Gene.name)%>%
    filter(padj<0.05&log2FoldChange>1)%>%
    top_n(20, wt=log2FoldChange)
  dup<-bind_rows(dup, de)
}

dwn<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select(ID, log2FoldChange,padj, Gene.name)%>%
    filter(padj<0.05&log2FoldChange< -1)%>%
    top_n(18, wt=abs(log2FoldChange))
  dwn<-bind_rows(dwn, de)
}


##htseq-count input

directory <-"/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/hit-counts/test/KOOE"

sampleFiles <- grep("txt",list.files(directory),value=TRUE)


sampleTable <- data.frame(sampleName =sampleFiles,
                          fileName = sampleFiles,
                          type =c("UnInf","UnInf","UnInf", "7dpi", "7dpi", "7dpi", "14dpi", "14dpi", "14dpi", "30dpi", "30dpi", "30dpi"),
                          group = c("UnInf","UnInf","UnInf","Inf","Inf","Inf","Inf","Inf","Inf","Inf","Inf","Inf"))

sampleTable$group<-factor(sampleTable$group)
sampleTable$type<-factor(sampleTable$type)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ group)
ddsHTSeq

##
dds<-collapseReplicates(ddsHTSeq, groupby=ddsHTSeq@colData@listData[["type"]], renameCols = TRUE)

##normalize

vsd <- vst(dds, blind=FALSE)

head(assay(vsd), 3)

##Interaction
dup<-dup%>%distinct(Gene.name, .keep_all = TRUE)
###
dwn<-dwn%>%distinct(Gene.name, .keep_all = TRUE)
select<-dup$ID
mt<-data.frame(t(scale(t(assay(vsd))))[select,])[, c(4,3,1,2)]
rownames(mt)<-dup$Gene.name
colnames(mt)<-c("control", "7dpi", "14dpi", "30dpi")
###

df<-data.frame(Group=c("control", "infected","infected","infected"))
rownames(df)<-colnames(mt)

pheatmap(mt, cluster_rows=FALSE, show_rownames=TRUE, 
         cluster_cols=FALSE, annotation_col=df)

##WT_OB
filenames<-c("WT-UnInf-OB-vs-WT-07dpi-OB",  "WT-UnInf-OB-vs-WT-30dpi-OB")
dup<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select(ID, log2FoldChange,padj, Gene.name)%>%
    filter(padj<0.05&log2FoldChange>1)%>%
    top_n(25, wt=log2FoldChange)
  dup<-bind_rows(dup, de)
}

dwn<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select(ID, log2FoldChange,padj, Gene.name)%>%
    filter(padj<0.05&log2FoldChange< -1)%>%
    top_n(29, wt=abs(log2FoldChange))
  dwn<-bind_rows(dwn, de)
}


##htseq-count input

directory <-"/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/hit-counts/test/WTOB"

sampleFiles <- grep("txt",list.files(directory),value=TRUE)


sampleTable <- data.frame(sampleName =sampleFiles,
                          fileName = sampleFiles,
                          type =c("UnInf","UnInf","UnInf", "7dpi", "7dpi", "7dpi",  "30dpi", "30dpi", "30dpi"),
                          group = c("UnInf","UnInf","UnInf","Inf","Inf","Inf","Inf","Inf","Inf"))

sampleTable$group<-factor(sampleTable$group)
sampleTable$type<-factor(sampleTable$type)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ group)
ddsHTSeq

##
dds<-collapseReplicates(ddsHTSeq, groupby=ddsHTSeq@colData@listData[["type"]], renameCols = TRUE)

##normalize

vsd <- vst(dds, blind=FALSE)

head(assay(vsd), 3)

##Interaction
dup<-dup%>%distinct(Gene.name, .keep_all = TRUE)
###
dwn<-dwn%>%distinct(Gene.name, .keep_all = TRUE)
select<-dwn$ID
mt<-data.frame(t(scale(t(assay(vsd))))[select,])[, c(3,2,1)]
rownames(mt)<-dwn$Gene.name
colnames(mt)<-c("control", "7dpi","30dpi")
###

df<-data.frame(Group=c("control", "infected","infected"))
rownames(df)<-colnames(mt)

pheatmap(mt, cluster_rows=FALSE, show_rownames=TRUE, 
         cluster_cols=FALSE, annotation_col=df)

##KO_OB
filenames<-c("IRF3-KO-UnInf-OB-vs-IRF3-KO-07dpi-OB", "IRF3-KO-UnInf-OB-vs-IRF3-KO-30dpi-OB")

dup<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select(ID, log2FoldChange,padj, Gene.name)%>%
    filter(padj<0.05&log2FoldChange>1)%>%
    top_n(25, wt=log2FoldChange)
  dup<-bind_rows(dup, de)
}

dwn<-data.frame()
for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))%>%
    select(ID, log2FoldChange,padj, Gene.name)%>%
    filter(padj<0.05&log2FoldChange< -1)%>%
    top_n(29, wt=abs(log2FoldChange))
  dwn<-bind_rows(dwn, de)
}


##htseq-count input

directory <-"/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/hit-counts/test/KOOB"

sampleFiles <- grep("txt",list.files(directory),value=TRUE)


sampleTable <- data.frame(sampleName =sampleFiles,
                          fileName = sampleFiles,
                          type =c("UnInf","UnInf","UnInf", "7dpi", "7dpi", "7dpi",  "30dpi", "30dpi", "30dpi"),
                          group = c("UnInf","UnInf","UnInf","Inf","Inf","Inf","Inf","Inf","Inf"))

sampleTable$group<-factor(sampleTable$group)
sampleTable$type<-factor(sampleTable$type)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ group)
ddsHTSeq

##
dds<-collapseReplicates(ddsHTSeq, groupby=ddsHTSeq@colData@listData[["type"]], renameCols = TRUE)

##normalize

vsd <- vst(dds, blind=FALSE)

head(assay(vsd), 3)

##Interaction
dup<-dup%>%distinct(Gene.name, .keep_all = TRUE)
###
dwn<-dwn%>%distinct(Gene.name, .keep_all = TRUE)
select<-dup$ID
mt<-data.frame(t(scale(t(assay(vsd))))[select,])[, c(3,2,1)]
rownames(mt)<-dup$Gene.name
colnames(mt)<-c("control", "7dpi","30dpi")
###

df<-data.frame(Group=c("control", "infected","infected"))
rownames(df)<-colnames(mt)

pheatmap(mt, cluster_rows=FALSE, show_rownames=TRUE, 
         cluster_cols=FALSE, annotation_col=df)
