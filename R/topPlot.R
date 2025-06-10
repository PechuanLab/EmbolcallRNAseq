#' Title
#'
#' @param signature1 DataFrame containing gene signatures with columns including 'symbol' and 'piFC'.
#' @param ExprMat Matrix of normalized expression values with genes as rows and samples as columns.
#' @param Treatment Vector indicating the treatment group for each sample in ExprMat.
#' @param contraste String representing the contrast condition for the plot title and file name.
#'
#' @return Jitter plot saved as a PDF file for the top 10 genes (5 upregulated and 5 downregulated).
#' @export
#'
#' @examples topPlot(signature1, ExprMat, Treatment, contraste)

topPlot <- function(signature1, ExprMat, Treatment, contraste) {
  # Select top 5 upregulated and top 5 downregulated genes based on piFC
  topUp = signature1 %>% dplyr::top_n(piFC, n = 5)
  topDown = signature1 %>% dplyr::top_n(piFC, n = -5)
  top10 = rbind(topUp, topDown)
  
  # Loop through the top 10 genes and create jitter plots
  for (j in 1:10) {
    gensmbl = top10$symbol[j]
    nice.col <- RColorBrewer::brewer.pal(6, name = "Spectral")
    
    # Save the plot as a PDF file
    pdf(paste(Prefix, contraste, gensmbl, ".pdf", sep = ""))
    print(graphics::stripchart(ExprMat[gensmbl, ] ~ Treatment, vertical = TRUE,
                               las = 2, cex.axis = 0.35, pch = 16, cex = 1.3,
                               col = nice.col, method = "jitter", ylab = "Normalised LogExpression",
                               main = paste(contraste, gensmbl, ".pdf", sep = "_")))
    dev.off()
  }
}
