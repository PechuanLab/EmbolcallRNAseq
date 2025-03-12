#' Remove Clustering Columns and Find Clusters
#'
#' This function removes columns from the metadata of a Seurat object that start with a specified prefix,
#' and then finds clusters at specified resolutions.
#'
#' @param seu A Seurat object.
#' @param prefix A string encoding the prefix of the columns to be removed from the metadata.
#' @param resolutions A numeric vector of resolutions at which to find clusters.
#' @return A Seurat object with updated clustering information.
#' @examples
#' @export
#' seu <- SeuratClustering(seu, "res_", c(0.1, 0.2, 0.3, 0.4, 0.5))
SeuratClustering <- function(seu, prefix, resolutions) {
  # Remove columns with the specified prefix from the metadata
  seu@meta.data <- seu@meta.data[, !grepl(paste0("^", prefix), colnames(seu@meta.data))]
  
  # Find clusters at specified resolutions
  for (res in resolutions) {
    seu <- FindClusters(seu, verbose = FALSE, resolution = res)
  }
  
  return(seu)
}
