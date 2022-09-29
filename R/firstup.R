

#' Title
#'
#' @param x
#'
#' @return firstup
#' @export
#'
#' @examples firstup( )
firstup <- function(x) {
  # First letter to capital
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


