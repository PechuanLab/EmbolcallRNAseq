#' Loads gene sets for gene onthology analysis and gene set enrichment anlysis.
#'
#' @param species. Dataset reference species either "Mus musculus" or "Homo sapiens".
#' @param date. To be downloaded rWikipPathWays database date.
#' @return A list of gene set databases.
#' @examples 
#' LoadGSEADB(species = "Mus musculus" ,date="20231010")
#' @export
LoadGSEADB <- function(species,date="current"){
  #Wiki Pathways
  wpgmt = rWikiPathways::downloadPathwayArchive(organism=species, format = "gmt", date= date)
  wp2gene = clusterProfiler::read.gmt(wpgmt) %>% 
  tidyr::separate(term, c("name","version","wpid","org"), "%")
  # MiSig Databases
  msigDB = msigdbr::msigdbr(species = species)
  # Keggg
  if(species == "Mus musculus"){
    dbGSEA = org.Mm.eg.db::org.Mm.eg.db
    keggdb = "mmu"
  } else if (species =="Homo sapiens"){
    dbGSEA = org.Hs.eg.db::org.Hs.eg.db
    keggdb = "hsa"
  }
  return(list(wp2gene,msigDB,dbGSEA,keggdb))
}
