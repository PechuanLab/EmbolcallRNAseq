#' Score samples using a gene set
#'
#' This function computes a score for each sample in an expression matrix based on a provided gene set.
#' The score can be summarized using either the first principal component (PC) or the mean expression.
#'
#' @param ExprMat A numeric matrix of gene expression values with genes as rows and samples as columns.
#' @param GeneList A character vector of gene symbols to be used for scoring.
#' @param summarizationFunction A string indicating the summarization method: "PC" (default) for principal component analysis, or "mean" for average expression.
#'
#' @return A numeric vector of scores for each sample.
#' @export
#'
#' @examples
#' # Example usage:
#' # gsScore(expression_matrix, c("GeneA", "GeneB", "GeneC"))
#' # gsScore(expression_matrix, gene_set, summarizationFunction = "mean")

gsScore <- function(ExprMat, GeneList, summarizationFunction = "PC") {
  # Restrict to genes present in the expression matrix
  GeneList <- GeneList[GeneList %in% rownames(ExprMat)]
  gm <- ExprMat[GeneList, ]
  
  if (nrow(gm) == 1) {
    gss <- gm %>% c()
  } else {
    if (summarizationFunction == "PC") {
      pc <- stats::prcomp(t(gm), retx = TRUE)
      gss <- pc$x[, 1] * sign(cor(pc$x[, 1], colMeans(gm)))
    } else {
      gss <- colMeans(gm)
    }
  }
  
  return(gss)
}
