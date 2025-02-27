#' Generate All Possible Contrast Combinations
#'
#' This function generates all possible contrast combinations for a given design matrix.
#' It takes a vector of column names from the design matrix and returns a matrix of contrasts.
#'
#' @param column_names A character vector containing the names of the columns in the design matrix.
#' @return A matrix of contrasts suitable for use with the limma package.
#' @examples
#' # Example usage
#' design_matrix_columns <- c("A", "B", "C")
#' contrast_matrix <- generate_contrasts(design_matrix_columns)
#' print(contrast_matrix)
#' @export
generate_contrasts <- function(column_names) {
  # Create all possible pairs of contrasts
  contrast_pairs <- combn(column_names, 2, simplify = FALSE)
  
  # Generate contrast expressions
  contrast_expressions <- sapply(contrast_pairs, function(pair) {
    paste(pair[1], "-", pair[2])
  })
  
  # Create the contrast matrix using makeContrasts
  contrast_matrix <- makeContrasts(contrasts = contrast_expressions, levels = column_names)
  
  return(contrast_matrix)
}
