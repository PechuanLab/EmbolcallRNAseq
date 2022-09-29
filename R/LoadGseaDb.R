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
  m_h =msigdbr::msigdbr(species = species,category = "H")
  m_kegg = msigdbr::msigdbr(species = species,subcategory = "CP:KEGG")
  m_c2 = msigdbr::msigdbr(species = species,category = "C2")
  m_c5 = msigdbr::msigdbr(species = species,category = "C5")
  m_c7 = msigdbr::msigdbr(species = species,category = "C7")
  if(species == "Mus musculus"){
    dbGSEA = org.Mm.eg.db::org.Mm.eg.db
    cellMarkers = 'http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Mouse_cell_markers.txt'
    keggdb = "mmu"
  } else if (species =="Homo sapiens"){
    dbGSEA = org.Hs.eg.db::org.Hs.eg.db
    cellMarkers = 'http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt'
    keggdb = "hsa"
  }
  return(list(wp2gene,m_h,m_kegg,m_c2,m_c5,m_c7,dbGSEA,cellMarkers,keggdb))
}
