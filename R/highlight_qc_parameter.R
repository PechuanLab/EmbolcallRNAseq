#' Highlights the cells to be removed/preserved after hard filtering on a particular quality control parameter
#' 
#' @param seurat_object. Seurat object with precomputed UMAP
#' @param qc_parameter. QC parameter on the metadata to color by.
#' @param threshold. Filter value.
#'
#' @return UMAP with hightlighted cells to be abova and below a QC threshold.
#' @export
#'
#' @examples highlight_qc_parameter()
#' 
highlight_qc_parameter <- function(seurat_object, qc_parameter, threshold) {
  # Check if the QC parameter exists in the metadata
  if (!qc_parameter %in% colnames(seurat_object@meta.data)) {
    stop(paste("QC parameter", qc_parameter, "not found in metadata"))
  }
  
  # Add a metadata column indicating whether the QC parameter is below the threshold
  seurat_object$highlight <- ifelse(seurat_object[[qc_parameter]] < threshold, "Below Threshold", "Above Threshold")
  
  # Use scCustomize to highlight these cells in a UMAP plot
  DimPlot_scCustom(seurat_object, group.by = "highlight") + ggtitle(qc_parameter)+
    scale_color_manual(values = c("Below Threshold" = "red", "Above Threshold" = "grey"))
}