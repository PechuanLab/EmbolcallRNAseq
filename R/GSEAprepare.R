

#' Title
#'
#' @param signature
#' @param Foldstatistic
#' @param species
#'
#' @return Gene vector
#' @export
#'
#' @examples GSEAPrepare( )

GSEAPrepare <- function(signature, Foldstatistic = "piFC",species ) {

  # Prepare for GSEA
  fold = signature[,Foldstatistic]
  gene = signature$symbol
  p_val_adj = signature$adj.P.Val
  dat = data.frame(LogFC =fold, Gene = gene, padj = p_val_adj)

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
