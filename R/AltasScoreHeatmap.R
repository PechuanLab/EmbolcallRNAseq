

#' Title
#'
#' @param ExprMat
#' @param pure_markers
#' @param Csplit
#' @param Rsplit
#' @param colnamesS
#'
#' @return A heatmap with PC score
#' @export
#'
#' @examples AtlasScoreHeatmap()
AtlasScoreHeatmap <- function(ExprMat,pure_markers,Csplit,Rsplit,colnamesS=F) {


  pc_matrix = PCscoreMat(ExprMat,pure_markers)
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
