#' Title
#'
#' @param PCAToolsObject
#' @param npcs
#' @param species
#' @return GSEA database
#' @export
#'
#' @examples  GseaPCA( )

GseaPCA <- function(PCAToolsObject, npcs = 3, species){
  
  # Get the PCA loadings matrix
  loadings_mat = PCAToolsObject$loadings
  
  for (i in 1:npcs) {
    # Prepare stats
    pc = data.frame(symbol = rownames(loadings_mat),logFC = loadings_mat[,i],adj.P.Val = 0)
    stats = GSEAPrepare(signature = pc , Foldstatistic = "logFC",species)
    ClusterProfilerOnthologies(stats,species,Prefix = paste0("PC",i),contraste =  paste0("PC",i)) 
  }
}

