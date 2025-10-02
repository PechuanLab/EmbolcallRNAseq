#' Prepare ranked gene list for GSEA
#'
#' Converts a differential expression signature table into a named numeric vector for GSEA, mapping gene symbols to Entrez IDs.
#'
#' @param signature Data frame. Must contain columns for gene symbols, log fold change, and adjusted p-value.
#' @param Foldstatistic Character. Name of the column to use as ranking statistic (default: "logFC").
#' @param species Character. Organism name ("Mus musculus" or "Homo sapiens").
#'
#' @return Named numeric vector of ranking statistics, names are Entrez gene IDs.
#' @export
#'
#' @examples
#' \dontrun{
#'   stats <- GSEAPrepare(signature, Foldstatistic = "logFC", species = "Mus musculus")
#' }

GSEAPrepare <- function(signature, Foldstatistic = "logFC",species ) {
  
  # Prepare for GSEA
  fold = signature[,Foldstatistic]
  gene = signature$symbol
  p_val_adj = signature$adj.P.Val
  dat = data.frame(LogFC = fold, Gene = gene, padj = p_val_adj)

  if(species == "Mus musculus"){
    dbGSEA = org.Mm.eg.db::org.Mm.eg.db
  } else if (species =="Homo sapiens"){
    dbGSEA = org.Hs.eg.db::org.Hs.eg.db
  }
  dat$GenENTRZ = AnnotationDbi::mapIds(dbGSEA, as.character(signature$symbol), 'ENTREZID', 'SYMBOL')
  dat = dat %>% tidyr::drop_na()
  #Prepare
  stats = dat$LogFC
  names(stats) = as.character(dat$GenENTRZ)
  stats = stats[order(stats, decreasing = T)]
  return(stats)
}
