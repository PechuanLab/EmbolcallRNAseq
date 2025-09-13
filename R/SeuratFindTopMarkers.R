##' Find top markers for each cluster in a Seurat object
##'
##' Processes a Seurat object by normalizing and scaling the data, then identifies the top marker genes for each cluster using FindAllMarkers. Returns a data frame of the top markers per cluster, filtered by adjusted p-value and sorted by log fold change.
##'
##' @param seu Seurat object. The Seurat object to process and analyze.
##'
##' @return Data frame. Top marker genes for each cluster, including statistics such as log fold change and adjusted p-value.
##' @export
##'
##' @examples
##' \dontrun{
##'   # Assuming 'seu' is your Seurat object
##'   top_markers <- SeuratFindTopMarkers(seu)
##'   head(top_markers)
##' }
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


