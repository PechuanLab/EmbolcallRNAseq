#' Scores samples with a genesets
#'
#' @param ExprMat
#' @param summarizationFunction
#' @param GeneList
#'
#' @return PC score
#' @export
#'
#' @examples gsScore( )

gsScore <- function(ExprMat,GeneList,summarizationFunction="PC") {
  # Restrict to rowspace of expression matrix
  GeneList = GeneList[GeneList %in% rownames(ExprMat)]
  gm  = ExprMat[GeneList,]
  if (nrow(gm) == 1){
    gss = gm %>% c()
  } else {
    if (summarizationFunction == "PC") {
      pc <- stats::prcomp(t(gm),
                   retx=TRUE)
      gss <- pc$x[,1] * sign(cor(pc$x[,1], colMeans(gm)))
    } else {
      gss <- colMeans(gm)
    }
  }
  return(gss)
}
