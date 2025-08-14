#' Assign multiple values to multiple variables
#'
#' This function acts as a flexible wrapper to assign values from a list or function to multiple variables.
#' It supports assigning a single value to one variable, or multiple values to multiple variables.
#' If fewer values are provided than variables, the remaining variables are assigned `NULL`.
#'
#' @param lhs A list of variable names (unquoted) to assign values to.
#' @param rhs A list of values, a function, or a formula to assign to the variables in `lhs`.
#'
#' @return Invisibly returns `NULL`. The primary effect is assignment in the parent frame.
#' @export
#'
#' @examples
#' # Assign multiple values
#' ListSupport(a, b, c) <- list(1, 2, 3)
#'
#' # Assign a single value
#' ListSupport(x) <- 42
#'
#' # Assign fewer values than variables
#' ListSupport(p, q, r) <- list("apple", "banana")
ListSupport <- function(lhs, rhs) {
  frame <- parent.frame()
  lhs <- as.list(substitute(lhs))
  
  if (length(lhs) > 1)
    lhs <- lhs[-1]
  
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir = frame)
    return(invisible(NULL)) 
  }
  
  if (is.function(rhs) || is(rhs, 'formula'))
    rhs <- list(rhs)
  
  if (length(lhs) > length(rhs))
    rhs <- c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  
  for (i in seq_along(lhs))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir = frame)
  
  return(invisible(NULL)) 
}
