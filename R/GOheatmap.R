#' Generate a GO term-annotated heatmap for gene expression clusters
#'
#' This function creates a heatmap of a scaled gene expression matrix, clusters genes, annotates clusters with top GO terms, and highlights selected genes. It supports multiple enrichment methods and custom annotations.
#'
#' @param scaled_mat Numeric matrix. Scaled gene expression matrix (genes × samples).
#' @param Gomet Character. Enrichment method: "GO", "Wiki", "Kegg", or "Custom".
#' @param nk_Go Integer. Number of gene clusters for enrichment analysis.
#' @param organismRef Character. Organism name: "Mus musculus" or "Homo sapiens".
#' @param nHL Integer. Number of genes to highlight per cluster. Default is 5.
#' @param Csplit Integer or factor. Column split for ComplexHeatmap.
#' @param gf Numeric. Font size for gene highlight labels. Default is 5.
#' @param title Character. Title for the heatmap. Default is blank.
#' @param clusCol Logical. Whether to cluster columns. Default is TRUE.
#' @param top_annotation ComplexHeatmap annotation object. Additional top annotation for the heatmap.
#' @param cust_db Data frame. Custom TERM2GENE mapping for enrichment if Gomet = "Custom". Default is NULL.
#' @param hl Logical. Whether to highlight genes. Default is TRUE.
#' @param showrownames Logical. Whether to show row (gene) names. Default is FALSE.
#' @param showcolnames Logical. Whether to show column (sample) names. Default is FALSE.
#'
#' @return A ComplexHeatmap::Heatmap object with GO term row splits and gene highlights.
#' @export
#'
#' @examples
#' \dontrun{
#'   GOHeatMap(scaled_mat, Gomet = "GO", nk_Go = 8, organismRef = "Mus musculus", nHL = 5, Csplit = 2, gf = 5, title = "Top Genes", clusCol = TRUE, top_annotation = NULL)
#' }

GOHeatMap <- function(scaled_mat,Gomet,nk_Go,organismRef,
                      nHL=5,Csplit,gf = 5,title = "  ",
                      clusCol=T,top_annotation,cust_db=NULL,hl = T,
                     showrownames = F,showcolnames = F){
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
                  show_column_names = showcolnames,
                  cluster_columns = clusCol,
                  show_column_dend = T,
                  column_names_gp = gpar(fontsize = 5),
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
