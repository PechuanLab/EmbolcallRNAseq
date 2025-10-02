#' Sweep theta values for projection score analysis
#'
#' Calculates projection scores for a range of theta values, where theta determines the fraction of most variable features to use. For each theta, the function selects the top variable features, computes the projection score, and returns a data frame of results. Useful for evaluating the effect of feature selection on projection stability.
#'
#' @param DataMatrix Matrix. Numeric data matrix (features x samples) to analyze.
#' @param NPCs Integer. Number of principal components to use. Default is 3.
#' @param nboot Integer. Number of bootstrap iterations. Default is 100.
#' @param thetas Numeric vector. Fraction(s) of most variable features to use (between 0 and 1). Default is c(0.25, 0.5, 0.75, 1).
#'
#' @return Data frame. Columns: Thetas (fraction of features used), ProjectionScores (score for each theta).
#' @export
#'
#' @examples
#' # Example with random data
#' set.seed(123)
#' mat <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' ps.df <- ThetaProjectionScore(mat, NPCs = 3, nboot = 50, thetas = c(0.1, 0.25, 0.5, 1))
#' 
ThetaProjectionScore <- function(DataMatrix, NPCs = 3, nboot = 100, thetas = c(0.25, 0.5, 0.75, 1)) {
  # Calculate the projections score given the thetas
  projection_scores = 1:length(thetas)
  for (i in 1:length(thetas)) {
    nvar = floor(thetas[i] * nrow(DataMatrix))
    topdat = NthMostVariableFeatures(DataMatrix, nvar)
    projection_scores[i] = ProjectionScore(DataMatrix = topdat, NPCs, nboot)
  }
  # Add zero
  thetas = c(0, thetas)
  projection_scores = c(0, projection_scores)
  ps.df = data.frame(Thetas = thetas, ProjectionScores = projection_scores)
  return(ps.df)
}

