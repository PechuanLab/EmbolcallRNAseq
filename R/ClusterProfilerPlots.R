#' Title
#'
#' @param ewp2
#' @param Prefix
#' @param contraste
#' @param DataBase
#' @param geneList
#'
#' @return GSEA plot
#' @export
#'
#' @examples ClusterProfilerPlots()

ClusterProfilerPlots <- function(ewp2,Prefix,contraste,DataBase,geneList){
  # Aes
  aescnet <- function(g) {
    g$layers[[3]]$aes_params$size = 2
    return(g)
  }
  # Save and prepare for basic plot
  saveRDS(ewp2,paste(DataBase,Prefix,contraste,"GSEA.Rdata",sep="_"))
  fgseaRes = ewp2@result
  write.csv(fgseaRes,paste(Prefix,contraste,DataBase,"GSEATable.csv",sep="_"))
  fgseaRes= fgseaRes %>% dplyr::arrange(pvalue) %>% tibble::as_tibble() %>%
    dplyr::filter(p.adjust<0.2) %>% 
    dplyr::mutate(Sign = paste("Sign",sign(NES))) %>% 
    dplyr::arrange(desc(NES))

  # Get the top pathways
  fgseaResTidy = tidyGSEA(fgseaRes)
  # Barplot
  pdf(paste(DataBase,Prefix,contraste,"GSEA.pdf",sep="_"), width = 12)
  print(
    ggpubr::ggbarplot(fgseaResTidy, x = "Description", y = "NES",
              fill = "Sign",           # change fill color by mpg_level
              color = "white",            # Set bar border colors to white
              palette = c("navy","firebrick"),  
              sort.val = "asc",          # Sort the value in descending order
              sort.by.groups = FALSE,     # Don't sort inside each group
              x.text.angle = 90,          # Rotate vertically x axis texts
              ylab = "Normalized Enrichment Score",
              xlab = "Description",
              rotate = TRUE,
              ggtheme = theme_minimal()
    )+
      theme(legend.position = "none")
  )
  dev.off()

  # Order things a bit
  y1 = fgseaRes %>% filter(qvalue < 0.2) %>% filter(NES < 0) %>%  arrange(desc(abs(NES))) 
  y2 = fgseaRes %>% filter(qvalue < 0.2) %>% filter(NES > 0) %>%  arrange(desc(abs(NES))) 
  
  #Top GSEA
  pdf(paste(DataBase,Prefix,contraste,"NegtopGSEA.pdf",sep="_"))
  print(
    try(enrichplot::gseaplot2(ewp2, geneSetID = y1$ID[1:4],subplots=1:2,pvalue_table = TRUE,color = c( "#ED574EFF","firebrick","#FF7F00","gold")),silent =T)
  )
  dev.off()

  #Top GSEA
  pdf(paste(DataBase,Prefix,contraste,"PostopGSEA.pdf",sep="_"))
  print(
    try(enrichplot::gseaplot2(ewp2, geneSetID = y2$ID[1:4],subplots=1:2,pvalue_table = TRUE,color = c("royalblue","gray","#6F9DC4FF","navy")),silent =T)
  )
  dev.off()

  # Network
  g =  try(cnetplot(ewp2, foldChange=geneList, colorEdge = TRUE,showCategory = 4),silent =T)

  pdf(paste(DataBase,Prefix,contraste,"EnrichNetWork.pdf",sep="_"), width = 18, height = 18)
  print(
    g + theme(legend.position="none") 
  )
  dev.off()

  # HeatMap
  pdf(paste(DataBase,Prefix,contraste,"heatmap.pdf",sep="_"),width = 30)
  print(
  enrichplot::heatplot(ewp2, foldChange=geneList, showCategory=5)
    )
  dev.off()

}
