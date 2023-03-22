#' Title
#'
#' @param y.all
#' @param nsamples
#'
#' @return Projection Score Plot
#' @export
#'
#' @examples FilterLowlyExpressed(y.all,nsamples = 3)
#' 
FilterLowlyExpressed <- function(y.all,nsamples = 3){
  # Filter lowly expressed genes
  libsize = mean(y.all$samples$lib.size)
  cpmthresh = 15/libsize*10^6
  keep.exprs = rowSums(edgeR::cpm(y.all$counts) >= cpmthresh) > nsamples
  yf =  y.all[keep.exprs,, keep.lib.sizes=FALSE]
  return(yf)
}



