#' STAMP Atlas Score Heatmap. 
#' Provides a set of signature scores derived from the STAMP single cell RNA sequencing experiment to quantify relative abundance of the cell types composing the tumor microenvironment.
#' @param ExprMat. Normalized expression data matrix.
#' @param pure_markers. Marker set dataframe with two columns, one corresponding to cell types to score and another with the marker signature genes.
#' @param Csplit. Plotting parameter to passs to ComplexHeatmap. Split columns.
#' @param Rsplit. Plotting parameter to passs to ComplexHeatmap. Split rows.
#' @param colnamesS. Plotting parameter to passs to ComplexHeatmap. Show columns names.
#'
#' @return A clustered heatmap with PC score for each cell type and sample.
#' @export
#'
#' @examples AtlasScoreHeatmap(ExprMat,pure_markers,Csplit,Rsplit,colnamesS=F)
AtlasScoreHeatmap <- function(ExprMat,pure_markers,Csplit,Rsplit,colnamesS=F) {
  # Score the signatures
  pc_matrix = PCscoreMat(ExprMat,pure_markers)
  # Generate the plot
  h_pc=ComplexHeatmap::Heatmap(pc_matrix,
               show_column_names = colnamesS,
               cluster_columns = TRUE,
               show_column_dend = T,
               column_names_gp = gpar(fontsize = 10),
               column_split = Csplit,
               border_gp = gpar(col = "gray", lty = 2),
               column_title = paste0("PC score ","_matrix"),
               column_title_side = "top",
               show_row_names = TRUE,
               row_names_gp = gpar(fontsize = 7),
               show_row_dend = T,
               row_dend_side = "left",
               cluster_column_slices = T,
               row_title_rot = 0,
               row_split = Rsplit,
               row_names_side="right",
               name="Scaled PC",
               top_annotation = top_bar
  )
  return(h_pc)
}
