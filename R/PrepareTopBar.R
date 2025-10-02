#' Prepare top bar annotation for ComplexHeatmap
#'
#' This function creates a ComplexHeatmap top annotation object from discrete and continuous sample annotations and their color palettes.
#'
#' @param DiscreteAnnotations Data frame. Discrete (categorical) sample annotations (columns = annotation types, rows = samples).
#' @param DiscretePalettes List. Named list of color vectors for each discrete annotation.
#' @param ContinuousAnnotations Data frame or NULL. Continuous (numeric) sample annotations. Default is NULL.
#' @param ContinuousPalettes List or NULL. Named list of color functions for each continuous annotation. Default is NULL.
#'
#' @return A ComplexHeatmap::HeatmapAnnotation object suitable for use as the `top_annotation` argument in ComplexHeatmap.
#' @export
#'
#' @examples
#' \dontrun{
#'   PrepareTopBar(DiscreteAnnotations, DiscretePalettes, ContinuousAnnotations = NULL, ContinuousPalettes = NULL)
#' }

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
