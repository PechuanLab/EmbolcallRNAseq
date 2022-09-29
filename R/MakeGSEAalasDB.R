
#' Title
#'
#' @param marker_df
#' @param p_val_adj
#' @param avg_log2FC
#'
#' @return Output annotation
#' @export
#'
#' @examples MakeGSEAAtlasDB( )
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
