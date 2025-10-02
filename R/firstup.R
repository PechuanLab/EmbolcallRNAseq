#' Capitalize the First Letter of Each String
#'
#' Converts the first character of each element in a character vector to uppercase, leaving the rest unchanged.
#'
#' @param x Character vector. The strings to capitalize.
#'
#' @return Character vector with the first letter of each element capitalized.
#' @export
#'
#' @examples
#' firstup(c("gene", "example", "test"))
firstup <- function(x) {
  # First letter to capital
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


