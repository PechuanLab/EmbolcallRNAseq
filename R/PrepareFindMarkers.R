
##' Format Seurat FindMarkers output for Embolcall
##'
##' Renames columns and adds a gene symbol column to a Seurat FindMarkers result for downstream compatibility.
##'
##' @param FindMarkersResult Data frame. Output from Seurat's FindMarkers or FindAllMarkers function.
##'
##' @return Data frame with columns renamed to `logFC`, `adj.P.Val`, and a new `symbol` column with gene names.
##' @export
##'
##' @examples
##' \dontrun{
##'   formatted <- PrepareFindMarkers(FindMarkersResult)
##' }

PrepareFindMarkers <- function(FindMarkersResult) {
       # Format accordingly
       colnames(FindMarkersResult)[c(2,5)] = c("logFC","adj.P.Val")
       FindMarkersResult$symbol = rownames(FindMarkersResult)
       return(FindMarkersResult)
}
