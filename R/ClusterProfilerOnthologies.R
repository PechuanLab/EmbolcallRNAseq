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
  
  # Stats
  geneList = stats

  # WikiPathways
  wpid2name = wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
  wpid2geneonly = wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
  ewp2 =  try(clusterProfiler::GSEA(geneList, 
                   TERM2GENE = wpid2geneonly,
                   TERM2NAME = wpid2name, 
                   verbose=FALSE, 
                   pAdjustMethod = "fdr",
                   by = "fgsea",
                   pvalueCutoff = 0.2), silent =T)
  ewp2 = DOSE::setReadable(ewp2, dbGSEA, keyType = "ENTREZID")
  ClusterProfilerPlots(ewp2,Prefix,contraste,DataBase="Wiki",geneList)

  # MsigDB
  testfor = msigDB$gs_cat %>% unique()
  for (k in 1:length(testfor)) {
    # Test for all the msigDB gene sets
    subdbname = testfor[k]
    m_t2g = msigDB %>% filter(gs_cat == subdbname) %>%
    dplyr::select(gs_name, entrez_gene)
    gseares = try(clusterProfiler::GSEA(geneList,
      TERM2GENE=m_t2g,
      by = "fgsea",
      pvalueCutoff = 0.2,
      pAdjustMethod = "fdr"),silent =T)
    gseares = DOSE::setReadable(gseares, dbGSEA, keyType = "ENTREZID")
    write.csv(gseares,paste(Prefix,contraste,subdb,".csv",sep="_"))
    try(ClusterProfilerPlots(lol,Prefix,contraste,DataBase = subdb,geneList),silent =T)
  }

  # KEGG
  wut = try(clusterProfiler::gseKEGG(geneList = geneList,organism = keggdb,pvalueCutoff = 0.2),silent =T)
  wut = DOSE::setReadable(wut, dbGSEA, keyType = "ENTREZID")
  ClusterProfilerPlots(wut,Prefix,contraste,DataBase = "KEGG",geneList)

}
