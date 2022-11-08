#' Title
#'
#' @param Contrastes
#' @param lfc
#'
#'
#' @return Plot for lfc DE genes
#' @export
#'
#' @examples 

DEHeatMap <- function(Contrastes,lfc = 2,ExprMat) {
  # Focus on the lfc > 
  FitTreat = treat(fit2,lfc = lfc)
  signature = topTreat(FitTreat, coef= Contrastes, n=Inf,adjust.method="fdr",sort.by = "P",p.value=0.05) 
  topdat = ExprMat[rownames(signature),] 
  scaled_mat = t(scale(t(topdat)))

  # Annotate by Go Terms Kegg
  g = GOHeatMap(scaled_mat = scaled_mat,Gomet = "Kegg",nk = 8,
                organismRef = "Mus musculus",title = paste0("Top ",lfc,Contrastes," Genes"),
                nHL = 5,Csplit=2,top_annotation = top_bar)
  # PDF
  pdf(paste0(Contrastes,"TopDE_genes.pdf"),width =15)
  print(g)
  dev.off()
}

