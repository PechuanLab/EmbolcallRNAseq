% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyGSEA.R
\name{tidyGSEA}
\alias{tidyGSEA}
\title{Tidy GSEA Results by Selecting Top Pathways}
\usage{
tidyGSEA(fgseaRes, npaths = 10)
}
\arguments{
\item{fgseaRes}{A data frame or tibble containing the results from `fgsea::fgsea()`,
including columns `NES` and `p.adjust`.}

\item{npaths}{Integer. Number of top upregulated and downregulated pathways to return. Default is 10.}
}
\value{
A tibble containing the top `npaths` upregulated and downregulated pathways,
sorted by descending NES.
}
\description{
This function filters and selects the top upregulated and downregulated pathways
from a `fgsea` result table based on normalized enrichment score (NES) and adjusted p-value.
It returns a tidy tibble containing the most significant pathways.
}
\examples{
tidyGSEA(fgseaRes, npaths = 10)
}
