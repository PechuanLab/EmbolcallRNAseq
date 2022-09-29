#' Title
#'
#' @param scaled_mat
#' @param Gomet
#' @param nk_Go
#' @param organismRef
#' @param nHL
#' @param Csplit
#' @param gf
#' @param title
#' @param clusCol
#' @param top_annotation
#' @param cust_db
#' @param hl
#' @param showrownames
#'
#' @return GO Heatmap
#' @export
#'
#' @examples GOHeatMap( )

GOHeatMap <- function(scaled_mat,Gomet,nk_Go,organismRef,
                      nHL=5,Csplit,gf = 5,title = "  ",
                      clusCol=T,top_annotation,cust_db=NULL,hl = T,
                     showrownames = F){
  # Wraps up the former functions for compact coding, returns a ggplot2 complexheatmap
  # Annotate by Go Terms
  #  Get the corrrect species DB
  if(organismRef=="Mus musculus"){
    db=org.Mm.eg.db::org.Mm.eg.db
    orgKegg = "mmu"
  }else if(organismRef=="Homo sapiens")
  {
    db=org.Hs.eg.db::org.Hs.eg.db
    orgKegg = "hsa"}
  # Run GO term analysis
  eout = enrich_out(expmat=scaled_mat,nk=nk_Go,scale = FALSE,organism =organismRef,
                  top_go = 3,Method = Gomet,cust_db=cust_db)   ##
  # Prepare the heatmap annotation
  go_title= add_go_heat(marker_lab = eout$marker_label,enrich_list =eout$go_result,top_go_name = 1)
  # Save the outcome for details
  tosaveGo = SaveGo(go_title)
  write.csv(tosaveGo,paste0(Gomet,"KmeanGo.csv"))
  g_csv = go_enrich_save(enrich_list = eout$go_result,filename="Exploratory")
  # Highlight some n genes as an annotation
  geneVecIdx = HLGenes(n=nHL,g_csv,db,Gomet= Gomet)

  if (hl==T) {


    genesHL =ComplexHeatmap::rowAnnotation(foo =ComplexHeatmap::anno_mark(at = geneVecIdx,
                                            labels = rownames(scaled_mat)[geneVecIdx],
                                            labels_gp = gpar(fontsize = gf)))
  } else{

    genesHL = NULL
  }
  # Heatmap
  GoHM =ComplexHeatmap::Heatmap(scaled_mat,
                  show_column_names = F,
                  cluster_columns = clusCol,
                  show_column_dend = T,
                  column_names_gp = gpar(fontsize = 2),
                  column_split = Csplit,
                  #column_order = sam_order,
                  # column_order = col_order_id,
                  border_gp = gpar(col = "gray", lty = 2),
                  column_title = title,
                  column_title_side = "top",
                  show_row_names = showrownames,
                  right_annotation = genesHL,
                  row_names_gp = gpar(fontsize = 2.5),
                  show_row_dend = TRUE,
                  row_dend_side = "left",
                  row_split = go_title,
                  #cluster_rows =cluster_within_group(t(scaled_mat), group$cluster),
                  cluster_row_slices = TRUE,
                  cluster_column_slices = T,
                  row_title_rot = 0,
                  #row_title = cluster_name,
                  name="Scaled TMM",
                  top_annotation = top_annotation

  )
  return(GoHM)
}
