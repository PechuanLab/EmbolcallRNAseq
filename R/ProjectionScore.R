
##' Compute projection score with bootstrapping
##'
##' Calculates the projection score as the difference between the explained variance (AlphaPE) and its bootstrapped average (BootAlpha) for a given data matrix. Useful for evaluating the stability of principal component projections.
##'
##' @param DataMatrix Matrix. Numeric data matrix (samples x features) to compute the projection score on.
##' @param NPCs Integer. Number of principal components to use. Default is 3.
##' @param nboot Integer. Number of bootstrap iterations. Default is 100.
##'
##' @return Numeric. Fraction of variance for projection score bootstrapped average.
##' @export
##'
##' @examples
##' # Example with random data
##' set.seed(123)
##' mat <- matrix(rnorm(1000), nrow = 100, ncol = 10)
##' ProjectionScore(mat, NPCs = 3, nboot = 50)
ProjectionScore <- function(DataMatrix, NPCs = 3, nboot = 100) {
  # Tau
  tau = AlphaPE(DataMatrix, NPCs) - BootAlpha(DataMatrix, NPCs, nboot)
  return(tau)
}

