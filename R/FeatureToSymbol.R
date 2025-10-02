#' Update feature names in ExpressionSet to gene symbols
#'
#' Replaces row names in an ExpressionSet object with gene symbols from its featureData.
#'
#' @param es.all An ExpressionSet or SummarizedExperiment object with featureData containing a "symbol" column.
#'
#' @return The input object with updated row names (gene symbols).
#' @export
#'
#' @examples
#' \dontrun{
#'   FeaturesToSymbol(eset)
#' }
FeaturesToSymbol <- function(es.all) {
  # Change to symbol on an ESET
  annotations = es.all@featureData@data
  rownames(es.all) = scuttle::uniquifyFeatureNames(ID= rownames(es.all),names=annotations[rownames(es.all), "symbol"])
  return(es.all)
}
