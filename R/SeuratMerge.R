#' Merge Multiple Seurat Objects
#'
#' This function reads multiple Seurat objects from a list of file paths, merges them into a single Seurat object, and joins the layers.
#'
#' @param seurats A character vector of file paths to the Seurat objects to be merged.
#' @return A merged Seurat object with joined layers.
#' @export
#' @examples SeuratMerge(seurats)

SeuratMerge <- function(seurats) {
  # Read the first Seurat object
  seu <- readRDS(seurats[1])
  seu <- CreateSeuratObject(counts = seu[['RNA']]$counts, 
                            meta.data = seu@meta.data)
  
  # Loop through the remaining Seurat objects and merge them
  for (i in 2:length(seurats)) {
    seu1 <- readRDS(seurats[i])
    seu1 <- CreateSeuratObject(counts = seu1[['RNA']]$counts, 
                               meta.data = seu1@meta.data)
    seu <- merge(seu, seu1)
    rm(seu1)
  }
  
  # Join layers
  seu <- JoinLayers(seu)
  
  # Return the merged Seurat object
  return(seu)
}
