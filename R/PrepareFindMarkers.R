#' Title
#'
#' @param FindMarkersResult
#' @return Seurat to Embolcall
#' @export
#'
#' @examples PrepareFindMarkers()

PrepareFindMarkers <- function(FindMarkersResult) {
       # Format accordingly
       colnames(FindMarkersResult)[c(2,5)] = c("logFC","adj.P.Val")
       FindMarkersResult$symbol = rownames(FindMarkersResult)
       return(FindMarkersResult)
}
