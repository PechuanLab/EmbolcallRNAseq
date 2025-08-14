#' Compute sample scores from a gene expression matrix using legacy summarization
#'
#' This function calculates a score for each sample based on the expression of a set of genes.
#' If multiple genes are provided, the score is computed using either the first principal component (PC)
#' or the mean expression across genes. If only one gene is provided, its expression values are returned directly.
#'
#' @param gm A numeric matrix of gene expression values with genes as rows and samples as columns.
#'           This matrix should contain only the genes of interest.
#' @param summarizationFunction A string specifying the summarization method: `"PC"` (default) for principal component analysis,
#'                               or `"mean"` for average expression.
#'
#' @return A numeric vector of scores for each sample.
#' @export
#'
#' @examples
#' # Using principal component summarization
#' gsScoreLegacy(gm = expression_matrix_subset, summarizationFunction = "PC")
#'
#' # Using mean summarization
#' gsScoreLegacy(gm = expression_matrix_subset, summarizationFunction = "mean")
gsScoreLegacy <- function(gm, summarizationFunction = "PC") {
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
