#' Process Seurat Object and Find Top Markers
#'
#' This function processes a Seurat object by normalizing and scaling the data, 
#' and then finds the top markers for each cluster.
#'
#' @seu seu A Seurat object to be processed.
#'
#' @return A data frame of top markers for each cluster.
#' @export
#'
#' @examples
#' \dontrun{
#'   # Assuming 'seu' is your Seurat object
#'   top_markers <- SeuratFindTopMarkers(seu)
#' }
SeuratFindTopMarkers <- function(seu) {
  # Set default assay to RNA
  DefaultAssay(seu) <- "RNA"
  
  # Normalize data
  seu <- NormalizeData(seu, scale.factor = median(seu@meta.data$nCount_RNA))
  
  # Scale data
  seu <- ScaleData(seu)
  
  # Find median of nCount_RNA
  median_nCount_RNA <- median(seu@meta.data$nCount_RNA)
  
  # Find all markers
  Markers <- FindAllMarkers(seu, only.pos = TRUE, logfc.threshold = 0.7, min.pct = 0.7, assay = "RNA")
  
  # Filter markers by adjusted p-value
  Markers <- Markers %>% filter(p_val_adj < 0.05)
  
  # Find top markers for each cluster
  topMarkers <- Markers %>% group_by(cluster) %>% top_n(avg_log2FC, n = 10)
  
  return(topMarkers)
}


