#' Title
#'
#' @param signature
#' @param Prefix
#' @param contraste
#' @param species
#' @param Foldstatistic
#'
#' @return Master wrapper for multiple functions
#' @export
#'
#' @examples BulkPlot()

BulkPlot <- function(signature,Prefix,contraste,species,Foldstatistic = "logFC") {
	
	# Some formatting
	signature$symbol = rownames(signature)
	signature = signature %>% dplyr::mutate(piFC =  logFC *(-log10(adj.P.Val)))
	
	# Save significant genes
    signature1 = signature %>% dplyr::filter(adj.P.Val<0.05) %>% 
                     dplyr::arrange(desc(abs(piFC))) %>%  drop_na()
	write.csv(signature1,paste(Prefix,contraste,"DESignificant.csv",sep="_"))

	# Differential expression hmp
	DEHeatMap(Contrastes = contraste ,lfc = 2)
	# VolcanoPlot
    VolcanoWrap(signature,Prefix,contraste)
	# MA plot
    ma_wrap(signature,contraste,Prefix)
    # Gost
    GostWrap(signature=signature)

	# Prepare for GSEA
    stats = GSEAPrepare(signature, Foldstatistic ,species) 
	# Run GSEA
    ClusterProfilerOnthologies(stats,species,Prefix,contraste) 
}
