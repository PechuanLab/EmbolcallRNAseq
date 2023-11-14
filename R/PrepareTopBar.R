#' Prepares the annotations for plotting with CompleHeatmap
#'
#' @param DiscreteAnnotations
#' @param DiscretePalettes
#' @param ContinuousPalettes
#' @param ContinuousAnnotations
#' @return Complex Heatmap annotations
#' @export
#'
#' @examples PrepareTopBar()

PrepareTopBar <- function(DiscreteAnnotations,DiscretePalettes,ContinuousAnnotations=NULL,ContinuousPalettes=NULL) {
       # Name the palettes
       palette_names = colnames(DiscreteAnnotations)
       names(DiscretePalettes) = palette_names
       # unlist
       for (i in 1:length(DiscretePalettes)) {
         names(DiscretePalettes[[i]]) = levels(factor(DiscreteAnnotations[,colnames(DiscreteAnnotations)[i]]))

       }
       # For the continuous
       names(ContinuousPalettes) = colnames(ContinuousAnnotations)
       # Bind Palettes
       list_palettes = c(DiscretePalettes,ContinuousPalettes)
       # Bind annotations
       if (is.null(ContinuousAnnotations)) {
              annotations = DiscreteAnnotations}
       else {
       annotations = cbind(DiscreteAnnotations,ContinuousAnnotations)}
       # Generate Top Bar
       top_bar=ComplexHeatmap::HeatmapAnnotation(df=annotations,col=list_palettes)
      
}
