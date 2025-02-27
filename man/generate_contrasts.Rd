% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_contrasts.R
\name{generate_contrasts}
\alias{generate_contrasts}
\title{Generate All Possible Contrast Combinations}
\usage{
generate_contrasts(column_names)
}
\arguments{
\item{column_names}{A character vector containing the names of the columns in the design matrix.}
}
\value{
A matrix of contrasts suitable for use with the limma package.
}
\description{
This function generates all possible contrast combinations for a given design matrix.
It takes a vector of column names from the design matrix and returns a matrix of contrasts.
}
\examples{
# Example usage
design_matrix_columns <- c("A", "B", "C")
contrast_matrix <- generate_contrasts(design_matrix_columns)
print(contrast_matrix)
}
