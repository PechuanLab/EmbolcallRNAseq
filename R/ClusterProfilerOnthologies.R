#' Title
#'
#' @param stats
#' @param species
#' @param Prefix
#' @param contraste
#'
#' @return Plots
#' @export
#'
#' @examples ClusterProfilerOnthologies()

ClusterProfilerOnthologies <- function(stats,species,Prefix,contraste) {

  #GSEA
  geneList = stats
  wp2gene = wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
  wpid2gene = wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
  wpid2name = wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
  ewp2 =  try(clusterProfiler::GSEA(geneList, TERM2GENE = wpid2gene,
                   TERM2NAME = wpid2name, verbose=FALSE, pvalueCutoff = 0.2), silent =T)
  ewp2 = DOSE::setReadable(ewp2, dbGSEA, keyType = "ENTREZID")
  ClusterProfilerPlots(ewp2,Prefix,contraste,DataBase="Wiki",geneList)
  # Cell Markers
  cell_markers = vroom::vroom(cellMarkers) %>%
    tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>%
    dplyr::select(cellMarker, geneID) %>%
    dplyr::mutate(geneID = strsplit(geneID, ', '))
  y = try(clusterProfiler::GSEA(geneList, TERM2GENE=cell_markers, pvalueCutoff = 0.2),silent =T)
  y = DOSE::setReadable(y, dbGSEA, keyType = "ENTREZID")
  ClusterProfilerPlots(y,Prefix,contraste,DataBase ="cellMarkers",geneList)
  # misgiDB
  m_df <- msigdbr::msigdbr(species = species)
  testfor = c("H","C2","C3","C5","C7")
  for (k in 1:length(testfor)) {
    subdb = testfor[k]
    m_t2g = msigdbr::msigdbr(species = species, category = subdb) %>%
      dplyr::select(gs_name, entrez_gene)
    lol = try(clusterProfiler::GSEA(geneList,TERM2GENE=m_t2g,pvalueCutoff = 0.2),silent =T)
    lol = DOSE::setReadable(lol, dbGSEA, keyType = "ENTREZID")
    write.csv(lol,paste(Prefix,contraste,subdb,".csv",sep="_"))
    try(ClusterProfilerPlots(lol,Prefix,contraste,DataBase = subdb,geneList),silent =T)
  }
  # KEGG
  wut = try(clusterProfiler::gseKEGG(geneList = geneList,organism = keggdb,pvalueCutoff = 0.2),silent =T)
  wut = DOSE::setReadable(wut, dbGSEA, keyType = "ENTREZID")
  ClusterProfilerPlots(wut,Prefix,contraste,DataBase = "KEGG",geneList)

}
