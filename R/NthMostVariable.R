
##' Select the n most variable features (genes)
##'
##' Returns a matrix of the top n most variable rows (genes) from an expression matrix.
##'
##' @param mat Numeric matrix. Rows are genes/features, columns are samples.
##' @param nvar Integer. Number of most variable features to select.
##'
##' @return A matrix containing only the n most variable rows (genes).
##' @export
##'
##' @examples
##' \dontrun{
##'   NthMostVariableFeatures(expr_matrix, nvar = 500)
##' }

NthMostVariableFeatures <- function(mat,nvar) {
  # Gets the top n most variable genes
  var_genes = apply(mat, 1, var)
  var_gens = var_genes[order(var_genes,decreasing = T)]
  topgenes = var_gens[1:nvar]
  topdat = mat[names(topgenes),]
  return(topdat)
}
