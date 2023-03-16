#' Title
#'
#' @param PCAToolsObject
#' @param cutoff
#'
#' @return The talus plot of a PCAtol pca object
#' @export
#'
#' @examples TalusPlot()
#' 
TalusPlot <- function(PCAToolsObject,cutoff){
  # Rasmuss Hennigson Talus Plot
  sigma= PCAToolsObject$sdev
  n = length(PCAToolsObject$sdev)
  # Get the shifted singular values
  sigmashift =  c(0,sigma)
  sigma = c(sigma,0)
  Logk = log(sigmashift)-log(sigma)
  # Remove non-sensical numbers
  Logk =  Logk[c(-1,-length(Logk))]
  pcs =c(1:(n-1))
  df = data.frame(PC = pcs, Logk = Logk)
  ggline(df, x = "PC",y="Logk",color = "royalblue",xlab = "Principal Component")+
         ylab(TeX("log $\\frac{\\sigma_k}{\\sigma_{k+1}}$")) + ylim(0,Logk[1])+
    geom_vline(xintercept = cutoff,linetype="dotted", size =1)
  
}
