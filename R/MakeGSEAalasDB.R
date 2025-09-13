
##' Prepare marker annotation data frame for GSEA
##'
##' Filters marker data frame by adjusted p-value and log2 fold change, adds Entrez IDs, and returns annotated data frame for GSEA.
##'
##' @param marker_df Data frame. Must contain columns for gene, cluster, p_val_adj, and avg_log2FC.
##' @param p_val_adj Numeric. Adjusted p-value threshold for filtering. Default is 0.05.
##' @param avg_log2FC Numeric. Absolute log2 fold change threshold for filtering. Default is 1.
##'
##' @return Data frame with marker info and Entrez IDs, suitable for GSEA input.
##' @export
##'
##' @examples
##' \dontrun{
##'   MakeGSEAtlasDB(marker_df, p_val_adj = 0.05, avg_log2FC = 1)
##' }
MakeGSEAtlasDB <-function(marker_df,p_val_adj = 0.05,avg_log2FC=1){

  marker_df = tibble::as_tibble(marker_df)
  topMarkers =  marker_df %>% dplyr::group_by(cluster) %>%dplyr::filter(p_val_adj<0.05)
  topMarkers =topMarkers %>% dplyr::ungroup() %>% dplyr:filter(abs(avg_log2FC)>1)
  # Change to ENTREZID
  mouse_db=org.Mm.eg.db::org.Mm.eg.db
  symbol=topMarkers$gene
  gene_obj=AnnotationDbi :: select(mouse_db,keys=symbol,keytype = "SYMBOL",
                                   columns = c("ENTREZID","SYMBOL"))
  add_gid=cbind(topMarkers,gene_obj)
  add_gid=add_gid %>% tidyr::drop_na(ENTREZID)
  return(add_gid)
}
