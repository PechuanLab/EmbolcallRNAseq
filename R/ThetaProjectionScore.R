#' Title
#'
#' @param DataMatrix
#' @param NPCs
#' @param nboot
#' @param thetas
#'
#' @return Fraction of variance for projection score bootstrapped average
#' @export
#'
#' @examples ps.df = ThetaProjectionScore(DataMatrix,NPCs =  3, nboot = 100, thetas =c(0.005,0.01,0.02,0.05,0.08,0.1,0.12,0.15,0.17,0.18,0.2,0.23,0.25,0.3,0.5,0.75,1))
#' 
ThetaProjectionScore <- function(DataMatrix,NPCs =  3, nboot = 100, thetas =c(0.25,0.5,0.75,1)){
  
  # Calculate the projections score given the thetas
  projection_scores = 1:length(thetas)
  for (i in 1:length(thetas)) {
    nvar = floor(thetas[i]*nrow(DataMatrix))
    topdat = NthMostVariableFeatures(ExprMat,nvar) 
    projection_scores[i] = ProjectionScore(DataMatrix = topdat ,NPCs, nboot)
  }
  # Add zero
  thetas =c(0,thetas)
  projection_scores = c(0,projection_scores)
  ps.df = data.frame(Thetas = thetas, ProjectionScores = projection_scores)
  return(ps.df)
}

