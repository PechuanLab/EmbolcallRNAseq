#' Plot Projection Score Curve and Optimal Theta
#'
#' Plots the projection score as a function of theta (fraction of features used), overlays a smoothed spline, and marks the theta value that maximizes the projection score. Useful for visualizing variable selection in projection analysis.
#'
#' @param ps.df Data frame. Must contain columns \code{Thetas} (fraction of features used) and \code{ProjectionScores} (projection score for each theta).
#'
#' @return Invisibly returns the plot. Side effect: displays a ggplot2/ggpubr plot with the projection score curve and optimal theta.
#' @export
#'
#' @examples
#' \dontrun{
#'   ProjectionScorePlot(ps.df)
#' }
ProjectionScorePlot <- function(ps.df){
  # Get an approximate derivative
  spl = smooth.spline(ps.df$Thetas, ps.df$ProjectionScores, spar = 0.2)
  X = data.frame(t=seq(0,1.0,length=100) ) # make an ordered sequence
  names(X) = "Thetas"
  Y =  predict(spl,newdata=X) # calculate predictions for that sequence
  TagetTheta = Y$x[which.max(Y$y)]
  
  # Using ggscatter
  ggpubr::ggscatter(ps.df, x = "Thetas", y = "ProjectionScores",color = "royalblue") +
    geom_vline(xintercept = TagetTheta) +
    geom_line(aes(Y$x, Y$y), color = "royalblue") +
    xlab(TeX("$\\theta$")) +
    ylab("Projection Score") +
    annotate("text", x = TagetTheta + 0.025, y = 0.2, label = TagetTheta) +
    theme_linedraw() +
    labs_pubr()
}



