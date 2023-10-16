#' Alpha from the projection score calculation. 
#' Calculates the fraction of the variance explained by the chosen principal components and moderates it by taking the square root.
#' @param DataMatrix. Dataset in matrix form.
#' @param NPCs. Number of principal components that are considered informative for the score.
#'
#' @return alphaFinal. Fraction of variance explained by the principal components for the projection score calculation.
#' @export
#'
#' @examples AlphaPE(DataMatrix, NPCs = 3)
#' 
AlphaPE <- function(DataMatrix,NPCs =  3){
  # Calculate the svd
  svdRes = svd(DataMatrix)
  # Take Eigenvalues to estimate fraction of variance step.
  lambdas= (svdRes$d)^2
  alphak = sum(lambdas[1:NPCs])
  alphaFinal = (alphak/sum(lambdas))^(1/2)
  return(alphaFinal)
}
