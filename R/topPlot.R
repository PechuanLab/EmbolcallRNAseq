
#' Title
#'
#' @param signature1
#' @param ExprMat
#' @param Treatment
#' @param contraste
#'
#' @return Jitter plot
#' @export
#'
#' @examples topPlot( )

topPlot <- function(signature1,ExprMat,Treatment,contraste){
  topUp = signature1 %>% dplyr::top_n(piFC,n=5)
  topDown = signature1 %>% dplyr::top_n(piFC,n=-5)
  top10 = rbind(topUp,topDown)
  for (j in 1:10) {
    gensmbl=top10$symbol[j]
    nice.col <-RColorBrewer::brewer.pal(6,name="Spectral")
    pdf(paste(Prefix,contraste,gensmbl,".pdf",sep=""))
    print(graphics::stripchart(ExprMat[gensmbl,]~Treatment,vertical=TRUE,
                     las=2,cex.axis=0.35,pch=16,cex=1.3,
                     col=nice.col,method="jitter",ylab="Normalised LogExpression",
                     main=paste(contraste,gensmbl,".pdf",sep="_"))
    )
    dev.off()
  }
}
