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
    dplyr::filter(p.adjust<0.2) %>%dplyr::mutate(Sign = paste("Sign",sign(NES))) %>%dplyr::arrange(desc(NES))

  # Get the top pathways
  fgseaResTidy = tidyGSEA(fgseaRes)
  # Barplot
  pdf(paste(DataBase,Prefix,contraste,"GSEA.pdf",sep="_"), width = 12)
  print(
    ggpubr::ggbarplot(fgseaResTidy, x = "Description", y = "NES",
              fill = "Sign",           # change fill color by mpg_level
              color = "white",            # Set bar border colors to white
              palette = c("navy","firebrick"),            # jco journal color palett. see ?ggpar
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
}
