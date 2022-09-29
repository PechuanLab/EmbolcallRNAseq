#' Title
#'
#' @param annotations
#' @param list_palettes
#' @return Heatmap annotations
#' @export
#'
#' @examples PrepareTopBar()

PrepareTopBar <- function(annotations,list_palettes) {
       # Name the palettes
       palette_names = colnames(annotations)
       names(list_palettes) = palette_names
       # unlist
       for (i in 1:length(list_palettes)) {
         names(list_palettes[[i]]) = levels(factor(annotations[,colnames(annotations)[i]]))

       }
       top_bar=ComplexHeatmap::HeatmapAnnotation(df=annotations,col=list_palettes)
       return(top_bar)
}
