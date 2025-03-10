#' Plot PCA Standard Deviations from a Seurat Object
#'
#' This function takes a Seurat object and plots the standard deviations of the PCA reduction to aid in selecting the optimal number of principal components.
#'
#' @param seu A Seurat object containing PCA reduction.
#' @param ndims An integer specifying the number of dimensions to plot. Default is 30.
#'
#' @return A plot of the PCA standard deviations.
#' @export
#'
#' @examples
#' \dontrun{
#'   SeuratTalusPlot(seurat_object, 30)
#' }
SeuratTalusPlot <- function(seu, ndims = 30) {
  # Extract standard deviations from PCA reduction
  primer <- data.frame(sdev = seu@reductions$pca@stdev)
  # Plot using TalusPlot
  TalusPlot(primer, ndims)
}

