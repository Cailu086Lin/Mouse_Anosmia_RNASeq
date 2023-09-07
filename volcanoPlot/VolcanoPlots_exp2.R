library(ggplot2)
library(gglabeller)
library(ggrepel)
library(tidyr)

filenames<-c("WT-UnInf-OE-vs-WT-07dpi-OE", "WT-UnInf-OE-vs-WT-14dpi-OE", "WT-UnInf-OE-vs-WT-30dpi-OE", "IRF3-KO-UnInf-OE-vs-IRF3-KO-07dpi-OE", "IRF3-KO-UnInf-OE-vs-IRF3-KO-14dpi-OE",
             "IRF3-KO-UnInf-OE-vs-IRF3-KO-30dpi-OE", "WT-UnInf-OE-vs-IRF3-KO-UnInf-OE", "WT-07dpi-OE-vs-IRF3-KO-07dpi-OE", "WT-14dpi-OE-vs-IRF3-KO-14dpi-OE", "WT-30dpi-OE-vs-IRF3-KO-30dpi-OE")

for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))

  p<- ggplot()+geom_point(de%>% drop_na()%>%filter(padj>0.05|abs(log2FoldChange) < 1), mapping=aes(log2FoldChange, -log10(padj)), alpha=1) +
                                geom_point(de%>% drop_na()%>%filter(padj<0.05&abs(log2FoldChange) > 1), mapping=aes(log2FoldChange, -log10(padj)), color="darkgrey", alpha=0.8)+ 
                                             ylab("Adjusted P value, -Log10")+xlab("Log2FoldChange")+
    geom_text_repel(subset(de, padj<0.05&log2FoldChange < -1), mapping=aes(log2FoldChange, -log10(padj),label = Gene.name), color="blue")+
    geom_text_repel(subset(de, padj<0.05&log2FoldChange >1), mapping=aes(log2FoldChange, -log10(padj),label = Gene.name), color="red")+
    theme_bw()+geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red",linewidth=0.2)+
    geom_vline(xintercept = 1, linetype="dashed", color="red",linewidth=0.2)+geom_vline(xintercept = -1, linetype="dashed", color="red", linewidth=0.2)+
    theme(legend.position="none", axis.title=element_text(size=14, face="bold"), axis.text = element_text(size=12))
  dev.new()
  ggsave(filename =paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/vlc_",filenames[i], ".tiff"),p, width =16, height =14, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")
  dev.off()
}
  