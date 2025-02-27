#' Master Wrapper for Multiple Functions
#'
#' This function serves as a master wrapper for multiple functions related to differential gene expression analysis.
#'
#' @param signature A data frame containing gene expression data, including columns for log fold change (logFC) and adjusted p-values (adj.P.Val).
#' @param Prefix A character string to be used as a prefix for output files.
#' @param contraste A character string specifying the contrast condition.
#' @param species A character string indicating the species (e.g., "human", "mouse").
#' @param Foldstatistic A character string specifying the fold change statistic to use (default is "logFC").
#' @return This function does not return a value but generates several output files and plots.
#' @export
#'
#' @examples
#' # Example usage:
#' BulkPlot(signature, Prefix = "Sample", contraste = "ConditionA_vs_ConditionB", species = "human")
BulkPlot <- function(signature, Prefix, contraste, species, Foldstatistic = "logFC") {
	
	# Some formatting
	signature$symbol = rownames(signature)
	signature = signature %>% dplyr::mutate(piFC =  logFC * (-log10(adj.P.Val)))
	
	# Save significant genes
  	signature1 = signature %>% dplyr::filter(adj.P.Val < 0.05) %>% 
                     dplyr::arrange(desc(abs(piFC))) %>% drop_na()
	write.csv(signature1, paste(Prefix, contraste, "DESignificant.csv", sep = "_"))

	# Differential expression heatmap
	try(DEHeatMap(Contrastes = contraste, lfc = 2, ExprMat), silent = TRUE)
	
	# Volcano plot
    VolcanoWrap(signature, Prefix, contraste)
	
	# MA plot
    ma_wrap(signature, contraste, Prefix)
	
	# Gene set overrepresentation analysis
    try(GostWrap(signature = signature), silent = TRUE)

	# Prepare for GSEA
    stats = GSEAPrepare(signature, Foldstatistic, species)
	
	# Run GSEA
	ClusterProfilerOnthologies(stats, species, Prefix, contraste)
}
