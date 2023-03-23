#' Title
#'
#' @param gm
#' @param summarizationFunction
#'
#' @return PC score
#' @export
#'
#' @examples gsScore( )

gsScore <- function(gm, summarizationFunction="PC") {

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
