
#' Title
#'
#' @param mat
#' @param nvar
#'
#' @return Matrix
#' @export
#'
#' @examples NthMostVariableFeatures( )

NthMostVariableFeatures <- function(mat,nvar) {
  # Gets the top n most variable genes
  var_genes = apply(mat, 1, var)
  var_gens = var_genes[order(var_genes,decreasing = T)]
  topgenes = var_gens[1:nvar]
  topdat = mat[names(topgenes),]
  return(topdat)
}
