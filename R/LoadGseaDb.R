#' Title
#'
#' @param species
#'
#' @return GSEA database
#' @export
#'
#' @examples  LoadGSEADB( )
LoadGSEADB <- function(species){

  # Load all databases
  #Wiki Pathways
  wpgmt = rWikiPathways::downloadPathwayArchive(organism=species, format = "gmt")
  wp2gene = clusterProfiler::read.gmt(wpgmt)
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
