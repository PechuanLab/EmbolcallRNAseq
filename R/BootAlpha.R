#' Title
#'
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
  
  # Calculate alphas
  alphas = 1:nboot
  for (i in 1:nboot) {
  shuf = matrix(resample(DataMatrix),nr)
  alphas[i] = AlphaPE(shuf,NPCs)
    
  }
  alpha_bootstrap = mean((alphas)^(1/2))
  return(alpha_bootstrap)
}


