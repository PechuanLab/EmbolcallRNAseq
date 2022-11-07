#' Wrapper for gostplot
#'
#' @param signature
#' @param species
#'
#' @return Wrapper for gostplot
#' @export
#'
#' @examples GostWrap()

GostWrap <- function(signature,species =  "mmusculus") {
	# Get significant genes
	significant =  signature %>% dplyr::filter(adj.P.Val < 0.05 & abs(logFC)>1)
	# Calculate gost result
	gostres = gprofiler2::gost(query = rownames(significant), organism = species)
	# Plot the top pathways for each category
	topPath = gostres$result %>% group_by(source) %>% slice_min(p_value,n = 1)
	p = gprofiler2::gostplot(gostres, capped = F, interactive = F)
	pp = gprofiler2::publish_gostplot(p, highlight_terms = topPath$term_id, 
                       width = NA, height = NA, filename = NULL )

    pdf("GOstPlot.pdf", width = 12, height =12)    
    print(pp)
    dev.off()              
	# Save significant pathways
	saveRDS(gostres,"GOstPlot_results.rds")

}
