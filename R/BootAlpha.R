#' Generates bootstrap samples of alpha for projection score calculation.
#' Takes a data matrix and resamples it to obtain a distribution for the fraction of variance explained by a given number of principal components.
#' @param DataMatrix
#' @param NPCs
#' @param nboot
#'
#' @return Fraction of variance for projection score bootstrapped average
#' @export
#'
#' @examples AlphaPE()
#' 
BootAlpha <- function(DataMatrix,NPCs =  3, nboot = 100){
  # Dimension
  nr = nrow(DataMatrix)
  # Calculate alphas
  alphas = 1:nboot
  for (i in 1:nboot) {
  shuf = matrix(sample(DataMatrix,replace = T),nr)
  alphas[i] = AlphaPE(shuf,NPCs)
    
  }
  alpha_bootstrap = mean((alphas)^(1/2))
  return(alpha_bootstrap)
}

