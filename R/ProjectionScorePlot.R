#' Title
#'
#' @param ps.df
#'
#' @return Projection Score Plot for quick variable selection
#' @export
#'
#' @examples ProjectionScorePlot(ps.df)
#' 
ProjectionScorePlot <- function(ps.df){
  # Get an approximate derivative
  spl = smooth.spline(ps.df$Thetas, ps.df$ProjectionScores, spar = 0.2)
  X = data.frame(t=seq(0,1.0,length=100) ) # make an ordered sequence
  names(X) = "Thetas"
  Y =  predict(spl,newdata=X) # calculate predictions for that sequence
  TagetTheta = Y$x[which.max(Y$y)]
  
  # Using ggscatter
  ggpubr::ggscatter(ps.df, x = "Thetas", y = "ProjectionScores",color = "royalblue")+geom_vline(xintercept = TagetTheta)+
    geom_line(aes(Y$x,Y$y),color = "royalblue")+xlab(TeX("$\\theta$"))+ylab("Projection Score")+
    annotate("text", x = TagetTheta+0.025, y = 0.2, label = TagetTheta)+ theme_linedraw()+labs_pubr()
    
}



