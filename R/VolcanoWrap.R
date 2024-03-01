#' Wrapps EnhancedVolcanoPlot
#'
#' @param signature DE result.
#' @param Prefix  String to append to title.
#' @param contraste Title of the plot.
#' @param listofgenes List of genes to highlight.
#' @return Volcano plot
#' @export
#'
#' @examples VolcanoWrap(signature,Prefix,contraste,listofgenes = NULL)
#'
VolcanoWrap <- function(signature,Prefix,contraste,listofgenes = NULL) {
  # VolcanoPlot
  xli = max(abs(signature$logFC))+0.05
  yli = (-log10(min(abs(signature$adj.P.Val))))+0.05
  png(paste(Prefix,contraste,"volcano.png",sep=""), units="in", width=5, height=5, res=300)
  print(EnhancedVolcano::EnhancedVolcano(signature,
                        lab = signature$symbol,
                        x = 'logFC',
                        y = 'adj.P.Val',
                        title = contraste,
                        subtitle = "",
                        pCutoff = 0.05,
                        FCcutoff = 2,
                        pointSize = 0.3,
                        labSize = 3.0,
                        xlim =c(-xli,xli),
                        ylim = c(0,yli),
                        colAlpha = 0.3,
                        legendPosition = "none",
                        selectLab = listofgenes)
  )
  dev.off()

}
