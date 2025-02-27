% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BulkPlot.R
\name{BulkPlot}
\alias{BulkPlot}
\title{Master Wrapper for Multiple Functions}
\usage{
BulkPlot(signature, Prefix, contraste, species, Foldstatistic = "logFC")
}
\arguments{
\item{signature}{A data frame containing gene expression data, including columns for log fold change (logFC) and adjusted p-values (adj.P.Val).}

\item{Prefix}{A character string to be used as a prefix for output files.}

\item{contraste}{A character string specifying the contrast condition.}

\item{species}{A character string indicating the species (e.g., "human", "mouse").}

\item{Foldstatistic}{A character string specifying the fold change statistic to use (default is "logFC").}
}
\value{
This function does not return a value but generates several output files and plots.
}
\description{
This function serves as a master wrapper for multiple functions related to differential gene expression analysis.
}
\examples{
# Example usage:
BulkPlot(signature, Prefix = "Sample", contraste = "ConditionA_vs_ConditionB", species = "human")
}
