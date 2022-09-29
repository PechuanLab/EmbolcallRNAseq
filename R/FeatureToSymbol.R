
#' Title
#'
#' @param es.all
#'
#' @return Gene names
#' @export
#'
#' @examples FeaturesToSymbol()
FeaturesToSymbol <- function(es.all) {
  # Change to symbol on an ESET
  annotations = es.all@featureData@data
  rownames(es.all) = scuttle::uniquifyFeatureNames(ID= rownames(es.all),names=annotations[rownames(es.all), "symbol"])
  return(es.all)
}
