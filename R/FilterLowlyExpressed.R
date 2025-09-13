##' Filter lowly expressed genes from DGEList
##'
##' Removes genes with low counts per million (CPM) across samples from an edgeR DGEList object.
##'
##' @param y.all edgeR DGEList object containing counts and sample info.
##' @param nsamples Integer. Minimum number of samples in which a gene must be expressed above threshold. Default is 3.
##' @param nreads Integer. Minimum read count threshold for CPM calculation. Default is 15.
##'
##' @return A filtered DGEList object with lowly expressed genes removed.
##' @export
##'
##' @examples
##' \dontrun{
##'   FilterLowlyExpressed(y.all, nsamples = 3, nreads = 15)
##' }
FilterLowlyExpressed <- function(y.all,nsamples = 3,nreads = 15){
  # Filter lowly expressed genes
  libsize = mean(y.all$samples$lib.size)
  cpmthresh = nreads/libsize*10^6
  keep.exprs = rowSums(edgeR::cpm(y.all$counts) >= cpmthresh) > nsamples
  yf =  y.all[keep.exprs,, keep.lib.sizes=FALSE]
  return(yf)
}



