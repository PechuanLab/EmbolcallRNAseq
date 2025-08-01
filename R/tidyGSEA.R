#' Tidy GSEA Results by Selecting Top Pathways
#'
#' This function filters and selects the top upregulated and downregulated pathways
#' from a `fgsea` result table based on normalized enrichment score (NES) and adjusted p-value.
#' It returns a tidy tibble containing the most significant pathways.
#'
#' @param fgseaRes A data frame or tibble containing the results from `fgsea::fgsea()`,
#' including columns `NES` and `p.adjust`.
#' @param npaths Integer. Number of top upregulated and downregulated pathways to return. Default is 10.
#'
#' @return A tibble containing the top `npaths` upregulated and downregulated pathways,
#' sorted by descending NES.
#'
#' @export
#'
#' @examples tidyGSEA(fgseaRes, npaths = 10)

tidyGSEA <- function(fgseaRes, npaths = 10) {
  topUP <- fgseaRes %>%
    dplyr::filter(NES > 0) %>%
    dplyr::top_n(n = npaths, wt = p.adjust / abs(NES)) %>%
    tibble::as_tibble()

  topDown <- fgseaRes %>%
    dplyr::filter(NES < 0) %>%
    dplyr::top_n(n = npaths, wt = p.adjust / abs(NES)) %>%
    tibble::as_tibble()

  fgseaResTidy <- rbind(topUP, topDown) %>%
    dplyr::arrange(desc(NES))

  return(fgseaResTidy)
}
