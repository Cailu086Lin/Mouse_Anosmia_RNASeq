library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr")
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

filenames<-c("WT-UnInf-OE-vs-WT-07dpi-OE", "WT-UnInf-OE-vs-WT-14dpi-OE", "WT-UnInf-OE-vs-WT-30dpi-OE", "IRF3-KO-UnInf-OE-vs-IRF3-KO-07dpi-OE", "IRF3-KO-UnInf-OE-vs-IRF3-KO-14dpi-OE",
             "IRF3-KO-UnInf-OE-vs-IRF3-KO-30dpi-OE", "WT-UnInf-OE-vs-IRF3-KO-UnInf-OE", "WT-07dpi-OE-vs-IRF3-KO-07dpi-OE", "WT-14dpi-OE-vs-IRF3-KO-14dpi-OE", "WT-30dpi-OE-vs-IRF3-KO-30dpi-OE")

for (i in 1:length(filenames)){
  de <- read.csv(paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/analysis_mouse/DEG/deseq2/", filenames[i],"/Differential_expression_analysis_table.csv"))
  gup<-subset(de, padj<0.05&log2FoldChange > 1)$Gene.name
  gdw<-subset(de, padj<0.05&log2FoldChange < -1)$Gene.name
  
  dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")
  if (websiteLive) {
    enriched <- enrichr(gup, dbs)
  }
  
  if (websiteLive) enriched[["GO_Cellular_Component_2023"]]
  if (websiteLive) enriched[["GO_Molecular_Function_2023"]]
  if (websiteLive) enriched[["GO_Biological_Process_2023"]]
  
  
  ggsave(filename =paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/pathwayUp_",filenames[i], ".tiff"),
         if (websiteLive) plotEnrich(enriched[[3]], showTerms = 20, numChar = 100, y = "Count", orderBy =  "P.value"), 
         width =24, height =14, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")
  
  dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")
  if (websiteLive) {
    enriched <- enrichr(gdw, dbs)
  }
  
  if (websiteLive) enriched[["GO_Cellular_Component_2023"]]
  if (websiteLive) enriched[["GO_Molecular_Function_2023"]]
  if (websiteLive) enriched[["GO_Biological_Process_2023"]]
  
  
  ggsave(filename =paste0("/media/cailu/Seagate Backup Plus Drive1/data/Wang_anosmiaRNASeq/pathwayDown_",filenames[i], ".tiff"),
         if (websiteLive) plotEnrich(enriched[[3]], showTerms = 20, numChar = 100, y = "Count", orderBy =  "P.value"), 
         width =24, height =14, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")
  
}








