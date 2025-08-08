#' Perform Gene Set Enrichment Analysis (GSEA) and Generate Plots
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) using WikiPathways and MSigDB gene sets, and generates plots to visualize the results. It also writes the results to CSV files for MSigDB gene sets.
#'
#' @param stats A named numeric vector representing the gene list for GSEA. The names should correspond to gene identifiers, and the values should represent the ranking metric (e.g., log fold change or statistical significance).
#' @param species A string specifying the species being analyzed (e.g., "human", "mouse"). This parameter is not directly used in the function but may be relevant for upstream data preparation.
#' @param Prefix A string used as a prefix for output file names and plot titles.
#' @param contraste A string specifying the contrast or condition being analyzed. This is used in output file names and plot titles.
#'
#' @return This function generates plots and writes CSV files as output. It does not return any R objects directly.
#' @export
#'
#' @examples
#' # Example usage:
#' # Assuming `stats` is a named numeric vector of gene rankings,
#' # `species` is "human", `Prefix` is "Experiment1", and `contraste` is "ConditionA":
#' ClusterProfilerOnthologies(stats, "human", "Experiment1", "ConditionA")
#'
ClusterProfilerOnthologies <- function(stats, species, Prefix, contraste) {
  # Stats
  geneList = stats

  # WikiPathways
  wpid2name = wp2gene %>% dplyr::select(wpid, name) # TERM2NAME
  wpid2geneonly = wp2gene %>% dplyr::select(wpid, gene) # TERM2GENE
  ewp2 = try(clusterProfiler::GSEA(
    geneList, 
    TERM2GENE = wpid2geneonly,
    TERM2NAME = wpid2name, 
    verbose = FALSE, 
    pAdjustMethod = "fdr",
    by = "fgsea",
    pvalueCutoff = 0.2
  ), silent = TRUE)

  # Change to symbol
  ewp2 = DOSE::setReadable(ewp2, dbGSEA, keyType = "ENTREZID")
  ClusterProfilerPlots(ewp2, Prefix, contraste, DataBase = "Wiki", geneList)

  # MsigDB
  testfor = msigDB$gs_cat %>% unique()

  for (k in 1:length(testfor)) {
    # Test for all the msigDB gene sets
    subdbname = testfor[k]
    m_t2g = msigDB %>% dplyr::filter(gs_cat == subdbname) %>%
      dplyr::select(gs_name, entrez_gene)

    # GSEA
    gseares = try(clusterProfiler::GSEA(
      geneList,
      TERM2GENE = m_t2g,
      by = "fgsea",
      pvalueCutoff = 0.2,
      pAdjustMethod = "fdr"
    ), silent = TRUE)

    gseares = DOSE::setReadable(gseares, dbGSEA, keyType = "ENTREZID")
    write.csv(gseares, paste(Prefix, contraste, subdbname, ".csv", sep = "_"))

    # Plots
    try(ClusterProfilerPlots(gseares, Prefix, contraste, DataBase = subdbname, geneList), silent = TRUE)
  }
}