
##' Run and plot GO/Pathway enrichment with gprofiler2::gost
##'
##' This function filters a differential expression signature for significant genes, runs gprofiler2::gost for pathway enrichment, plots the top pathways, and saves both the plot and results.
##'
##' @param signature Data frame. Differential expression results with rownames as gene symbols and columns including `logFC` and `adj.P.Val`.
##' @param species Character. Organism code for gprofiler2 (e.g., "mmusculus" for mouse, "hsapiens" for human). Default is "mmusculus".
##'
##' @return Invisibly returns NULL. Side effects: saves a PDF plot ("GOstPlot.pdf") and an RDS file ("GOstPlot_results.rds").
##' @export
##'
##' @examples
##' \dontrun{
##'   GostWrap(signature, species = "mmusculus")
##' }

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
