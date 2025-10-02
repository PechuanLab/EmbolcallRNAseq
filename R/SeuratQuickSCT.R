#' Process Seurat Object with SCTransform, PCA, and UMAP
#'
#' This function processes a Seurat object by normalizing the data, scoring cell cycles, 
#' applying SCTransform, running PCA, finding neighbors, running UMAP, and finding clusters 
#' at multiple resolutions.
#'
#' @param seu A Seurat object to be processed.
#' @param ndims An integer specifying the number of dimensions to use for PCA and UMAP. Default is 30.
#' @param npcs Number of principal components to be calculated.
#' @param species Species for cell cycle scoring. One of "human" or "mouse". Default is "human".
#'
#' @return A processed Seurat object with normalized data, cell cycle scores, PCA, UMAP, and clusters.
#' @export
#'
#' @examples
#' \dontrun{
#'   # Assuming 'seurat_object' is your Seurat object
#'   seurat_object <- SeuratQuickSCT(seurat_object, ndims = 30, npcs = 50)
#' }
SeuratQuickSCT <- function(seu, ndims = 30, npcs = 50, species = c("human", "mouse")) {
  species <- match.arg(species)
  # Normalize data
  seu <- NormalizeData(seu, scale.factor = median(seu@meta.data$nCount_RNA))
  
  # Cell cycle scoring
  if (species == "mouse") {
    # Seurat's recommended approach for mouse is to title-case human CC genes
    # https://satijalab.org/seurat/articles/cell_cycle_vignette.html
    s.genes <- cc.genes$s.genes %>% tolower() %>% firstup()
    g2m.genes <- cc.genes$g2m.genes %>% tolower() %>% firstup()
  } else {
    # For human, use the original human gene symbols
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
  }
  seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
  
  # SCTransform
  seu <- SCTransform(seu, verbose = TRUE, conserve.memory = TRUE, 
                     vars.to.regress = c("S.Score", "G2M.Score",
                                         "percent_ribo", "percent_mito"))
  
  # Run PCA
  seu <- RunPCA(seu, npcs = npcs)
  
  # Find neighbors and run UMAP
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:ndims)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:ndims)
  
  # Clean up metadata
  seu@meta.data <- seu@meta.data[, !grepl(paste0("^", "SCT_"), colnames(seu@meta.data))]
  
  # Find clusters at multiple resolutions
  resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2)
  for (res in resolutions) {
    seu <- FindClusters(seu, verbose = FALSE, resolution = res)
  }
  
  return(seu)
}

