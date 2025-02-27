#' Loops BulkPlot through a set of limma outputs
#'
#' This function loops through a set of contrasts from limma outputs and generates various plots and files for differential gene expression analysis.
#'
#' @param ExprMat Normalized expression matrix.
#' @param Prefix String to append to the contrasts.
#' @param Contrastes A vector of contrast names.
#' @param annotations Annotations for the genes (not used in the function but can be included for completeness).
#' @param Treatment Treatment conditions (not used in the function but can be included for completeness).
#' @param species A character string indicating the species (e.g., "human", "mouse").
#' @param logfold Log fold change threshold for differential expression.
#'
#' @return This function does not return a value but generates several output files and plots.
#' @export
#'
#' @examples
#' # Example usage:
#' BulkPlots(ExprMat, Prefix = "Sample", Contrastes = c("ConditionA_vs_ConditionB", "ConditionC_vs_ConditionD"), annotations, Treatment, species = "human", logfold = 1)
BulkPlots <- function(ExprMat, Prefix, Contrastes, annotations, Treatment, species, logfold) {

  # Put everything in a folder and be tidy
  dir = getwd()
  onestan = paste(getwd(), paste(Prefix, "DEPlots", sep = "_"), sep = "/")
  dir.create(onestan)
  setwd(onestan)

  # Loop for all contrasts
  for (i in 1:length(Contrastes)) {
    # Get the contrast and create the directory
    contraste = Contrastes[i]
    onestan2 = paste(onestan, contraste, sep = "/")
    dir.create(onestan2)
    setwd(onestan2)
  
    # Prepare the fitted limma object
    fit3 = limma::treat(fit2, lfc = 0)
    signature = limma::topTreat(fit3, coef = contraste, n = Inf, adjust.method = "fdr", sort.by = "P") 
    # Save the whole thing in case you would like to run something else
    write.csv(signature, paste0(Prefix, contraste, "allGenes.csv"))
    # Bulk Plot
    try(BulkPlot(signature, Prefix, contraste, species), silent = TRUE)
  }
}

