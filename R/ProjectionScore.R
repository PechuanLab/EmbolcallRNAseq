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
ProjectionScore <- function(DataMatrix,NPCs =  3, nboot = 100){
  # Tau
  tau = AlphaPE(DataMatrix,NPCs) - BootAlpha(DataMatrix,NPCs, nboot)
  return(tau)
}

