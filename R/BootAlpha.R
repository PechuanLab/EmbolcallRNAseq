#' Generates bootstrap samples of alpha for projection score calculation.
#' 
#' This function takes a data matrix and resamples it to obtain a distribution for the fraction of variance explained by a given number of principal components.
#' 
#' @param DataMatrix A numeric matrix where rows represent samples and columns represent variables.
#' @param NPCs An integer specifying the number of principal components to consider. Default is 3.
#' @param nboot An integer specifying the number of bootstrap samples to generate. Default is 100.
#' 
#' @return A numeric value representing the bootstrapped average of the fraction of variance for the projection score.
#' @export
#' 
#' @examples
#' \dontrun{
#' # Example usage
#' DataMatrix <- matrix(rnorm(100), nrow = 10)
#' result <- BootAlpha(DataMatrix, NPCs = 3, nboot = 100)
#' print(result)
#' }
BootAlpha <- function(DataMatrix, NPCs = 3, nboot = 100) {
  # Dimension
  nr = nrow(DataMatrix)
  # Calculate alphas
  alphas = numeric(nboot)
  for (i in 1:nboot) {
    shuf = matrix(sample(DataMatrix, replace = TRUE), nr)
    alphas[i] = AlphaPE(shuf, NPCs)
  }
  alpha_bootstrap = mean(sqrt(alphas))
  return(alpha_bootstrap)
}


