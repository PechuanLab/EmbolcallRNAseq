#' Title
#'
#' @param DataMatrix
#' @param NPCs
#'
#' @return Fraction of variance for projection score
#' @export
#'
#' @examples AlphaPE()
#' 
AlphaPE <- function(DataMatrix,NPCs =  3){
  # svd
  svdRes = svd(DataMatrix)
  # Rasmuss Hennigson Talus Plot
  lambdas= (svdRes$d)^2
  alphak = sum(lambdas[1:NPCs])
  alphaFinal = (alphak/sum(lambdas))^(1/2)
  return(alphaFinal)
}
