#' Loads gene sets for gene onthology analysis and gene set enrichment anlysis.
#'
#' @param species. Dataset reference species either "Mus musculus" or "Homo sapiens".
#' @param date. To be downloaded rWikipPathWays database date.
#' @return A list of gene set databases.
#' @examples 
#' LoadGSEADB(species = "Mus musculus" ,date="20241110")
#' @export
LoadGSEADB <- function(species, date = "20251010") {
  # Warn user to ensure the provided date is up-to-date
  warning(sprintf("Please ensure 'date' is up-to-date (YYYYMMDD). Provided date: %s", as.character(date)), call. = FALSE)

  #Wiki Pathways
  wpgmt = rWikiPathways::downloadPathwayArchive(organism = species, format = "gmt", date = date)
  wp2gene = clusterProfiler::read.gmt(wpgmt) %>% 
    tidyr::separate(term, c("name","version","wpid","org"), "%")
  # MiSig Databases
  msigDB = msigdbr::msigdbr(species = species)
  # Keggg
  if (species == "Mus musculus") {
    dbGSEA = org.Mm.eg.db::org.Mm.eg.db
    keggdb = "mmu"
  } else if (species == "Homo sapiens") {
    dbGSEA = org.Hs.eg.db::org.Hs.eg.db
    keggdb = "hsa"
  } else {
    stop("Unsupported species. Use 'Mus musculus' or 'Homo sapiens'.")
  }
  return(list(wp2gene, msigDB, dbGSEA, keggdb))
}
