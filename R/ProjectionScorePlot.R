#' Title
#'
#' @param ps.df
#'
#' @return Projection Score Plot
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
  
  
  
  lines(X$Thetas,Y,type="l",main="Original fit") # check if predictions are correct
  
  dY <- diff(Y)/diff(X$Thetas)  # calculate the dY values for our function
  dX <- rowMeans(embed(X$Thetas,2)) # centers the X values for plotting
  plot(dX,dY,type="l",main="Derivative") # check if derivative is correct
  deriv=predict(spl, deriv=1,newdata = c(0))
  # Using ggscatter
  ggscatter(ps.df, x = "Thetas", y = "ProjectionScores")+geom_vline(xintercept = TagetTheta)+
    geom_line(aes(Y$x,Y$y))+xlab(TeX("$\\theta$"))+ylab("Projection Score")+
    annotate("text", x = TagetTheta+0.025, y = 0.2, label = TagetTheta)+ theme_linedraw()+labs_pubr()
    
}



