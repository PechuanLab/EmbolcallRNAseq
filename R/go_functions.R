#' Save GO enrichment results to CSV
#'
#' Combines a list of GO enrichment results into a single data frame and writes it to a CSV file.
#'
#' @param enrich_list List of data frames, each containing GO enrichment results for a cluster.
#' @param filename Character string. Prefix for the output CSV file.
#'
#' @return A data frame of combined enrichment results (invisible). Also writes a CSV file to disk.
#' @export
#'
#' @examples
#' \dontrun{
#'   go_enrich_save(list_of_results, filename = "my_analysis")
#' }
go_enrich_save <- function(enrich_list, filename) {
  #  e.g  go_enrich_save(enrich_list = eout$go_result)
  ## this function is to output csv file of GO enrichment result ####
  ### input is list object ###
  num_len = length(enrich_list)
  go_path = paste0(paste(filename, "GO_out", "k", num_len, sep = "_"), ".csv")
  go_csv = dplyr::bind_rows(enrich_list, .id = "cluster_label")
  write.csv(go_csv, file = go_path)
  return(go_csv)
}

#' Annotate clusters with top GO term
#'
#' For each cluster, assigns the top GO term description to the marker label vector.
#'
#' @param marker_lab Integer or character vector. Cluster assignments for each gene.
#' @param enrich_list List of data frames, each containing GO enrichment results for a cluster.
#' @param top_go_name Integer. Which GO term to use (row index in enrichment result). Default is 1.
#'
#' @return A vector of GO term descriptions, same length as marker_lab.
#' @export
#'
#' @examples
#' \dontrun{
#'   add_go_heat(marker_lab, enrich_list, top_go_name = 1)
#' }
add_go_heat <- function(marker_lab, enrich_list, top_go_name = 1) {
  ## this function is used to produce stuff that will be used in Heatmap function ##
  ## for row split and row name. ####
  ## top_go_name is number of GO name ranging 1 to 5.
  ## If fdr p value are quite similar and small (e-8 or e-7), you may be interested in other GO name ##
  ## output is for option "row_split  ###
  c = length(enrich_list)
  tn = top_go_name
  g_in = marker_lab
  go_result = enrich_list
  for (i in 1:c) {
    g_in[g_in == i] = go_result[[i]]$Description[tn]
  }
  return(g_in)
}

#' Convert GO term list to tidy data frame
#'
#' Takes a list of GO term assignments and returns a tidy data frame with gene and GO term columns.
#'
#' @param go_title List or vector of GO term assignments per gene.
#'
#' @return A data frame with columns Gene and GOTerm.
#' @export
#'
#' @examples
#' \dontrun{
#'   SaveGo(go_title)
#' }
SaveGo <- function(go_title) {
  # Saves Go terms
  tosaveGo = do.call(cbind.data.frame, as.list(go_title))
  tosaveGoDF = tidyr::gather(tosaveGo, Gene, GOTerm, factor_key = TRUE) %>%
    dplyr::arrange(factor(GOTerm))
  return(tosaveGoDF)
}

#' Highlight top genes for annotation
#'
#' Selects top genes from GO enrichment results for annotation in heatmaps, optionally converting IDs.
#'
#' @param ngenes Integer. Number of top genes to select per cluster. Default is 5.
#' @param g_csv Data frame. Combined GO enrichment results (from go_enrich_save).
#' @param db OrgDb object. Organism database for ID conversion.
#' @param Gomet Character. Enrichment method ("Wiki", "Kegg", "Custom").
#'
#' @return Integer vector of row indices for selected genes in the scaled matrix.
#' @export
#'
#' @examples
#' \dontrun{
#'   HLGenes(ngenes = 5, g_csv, db, Gomet = "Wiki")
#' }
HLGenes <- function(ngenes = 5, g_csv, db, Gomet) {
  # Highlights genes randomly
  geneses = paste0("gene", 1:ngenes)
  genes = g_csv %>%
    dplyr::group_by(cluster_label) %>%
    dplyr::filter(FDR_Pvalue == min(FDR_Pvalue)) %>%
    dplyr::ungroup() %>%
    tidyr::separate(geneID, geneses, "/")
  geneVec = genes %>%
    dplyr::select(-cluster_label, -Description, -FDR_Pvalue) %>%
    unique() %>%
    as.matrix() %>%
    c() %>%
    na.omit() %>%
    as.character()
  if (Gomet %in% c("Wiki", "Kegg", "Custom")) {
    geneVec = clusterProfiler::bitr(gene = geneVec, fromType = "ENTREZID", toType = c("SYMBOL"),
                                   OrgDb = db, drop = TRUE)
    geneVec = geneVec$SYMBOL
  }
  geneVecIdx = which(rownames(scaled_mat) %in% geneVec)
  return(geneVecIdx)
}

#' Clean GSEA result names for readability
#'
#' Removes prefixes and underscores from GSEA result names, IDs, and descriptions for easier plotting.
#'
#' @param gsea_result GSEA result object (from clusterProfiler or similar).
#'
#' @return The modified GSEA result object with cleaned names and descriptions.
#' @export
#'
#' @examples
#' \dontrun{
#'   CleanName(gsea_result)
#' }
CleanName <- function(gsea_result) {
  # Clean GSEA name
  names(gsea_result@geneSets) = stringr::str_replace_all(stringr::str_remove(names(gsea_result@geneSets), "HALLMARK_"), "_", " ")
  gsea_result@result$ID = stringr::str_replace_all(stringr::str_remove(gsea_result@result$ID, "HALLMARK_"), "_", " ")
  gsea_result@result$Description = stringr::str_replace_all(stringr::str_remove(gsea_result@result$Description, "HALLMARK_"), "_", " ")
  return(gsea_result)
}
