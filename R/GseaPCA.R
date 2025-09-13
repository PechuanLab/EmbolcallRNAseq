
##' Run GSEA on principal components from PCA results
##'
##' For each of the top principal components, this function prepares a ranked gene list from the loadings, runs GSEA, and generates enrichment plots and tables.
##'
##' @param PCAToolsObject List or S3 object. Must contain a `loadings` matrix (genes × PCs) from PCA.
##' @param npcs Integer. Number of principal components to analyze. Default is 3.
##' @param species Character. Organism name (e.g., "Mus musculus", "Homo sapiens").
##'
##' @return Invisibly returns NULL. Side effects: generates GSEA plots and tables for each PC.
##' @export
##'
##' @examples
##' \dontrun{
##'   GseaPCA(PCAToolsObject, npcs = 3, species = "Mus musculus")
##' }

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

