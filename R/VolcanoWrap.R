#' Wrapper for EnhancedVolcano Plot
#'
#' This function creates and saves a volcano plot using the EnhancedVolcano package, based on differential expression results.
#'
#' @param signature Data frame. Differential expression results, must contain columns \code{logFC}, \code{adj.P.Val}, and \code{symbol}.
#' @param Prefix Character. String to prepend to the output file name.
#' @param contraste Character. Title of the plot and part of the output file name.
#' @param listofgenes Character vector or NULL. List of gene symbols to highlight on the plot (optional).
#'
#' @return Invisibly returns NULL. Side effect: saves a PNG file with the volcano plot.
#' @export
#'
#' @examples
#' \dontrun{
#'   VolcanoWrap(signature, Prefix = "Exp1_", contraste = "Tumor_vs_Normal", listofgenes = c("GeneA", "GeneB"))
#' }
VolcanoWrap <- function(signature, Prefix, contraste, listofgenes = NULL) {
  # VolcanoPlot
  xli = max(abs(signature$logFC)) + 0.05
  yli = (-log10(min(abs(signature$adj.P.Val)))) + 0.05
  png(paste(Prefix, contraste, "volcano.png", sep = ""), units = "in", width = 5, height = 5, res = 300)
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
                        xlim = c(-xli, xli),
                        ylim = c(0, yli),
                        colAlpha = 0.3,
                        legendPosition = "none",
                        selectLab = listofgenes)
  )
  dev.off()
  invisible(NULL)
}
