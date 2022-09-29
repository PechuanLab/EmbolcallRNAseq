

#' Title
#'
#' @param signature
#' @param Prefix
#' @param contraste
#'
#' @return Volcano plot
#' @export
#'
#' @examples VolcanoWrap( )
VolcanoWrap <- function(signature,Prefix,contraste) {
  # VolcanoPlot
  xli = max(abs(signature$logFC))+0.05
  yli = (-log10(min(abs(signature$adj.P.Val))))+0.05
  png(paste(Prefix,contraste,"volcano.png",sep=""), units="in", width=5, height=5, res=300)
  print(EnhancedVolcano::EnhancedVolcano(signature,
                        lab = signature$symbol,
                        x = 'logFC',
                        y = 'adj.P.Val',
                        title = "",
                        subtitle = "",
                        pCutoff = 0.05,
                        FCcutoff = 1,
                        pointSize = 0.3,
                        labSize = 3.0,
                        xlim =c(-xli,xli),
                        ylim = c(0,yli),
                        colAlpha = 0.3,
                        legendPosition = "none")
  )
  dev.off()

}
