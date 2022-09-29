#' Title
#'
#' @param enrich_list
#' @param filename
#'
#' @return Dowload files
#' @export
#'
#' @examples go_enrich_save( )
go_enrich_save<-function(enrich_list, filename){
  #  e.g  go_enrich_save(enrich_list = eout$go_result)
  ## this function is to output csv file of GO enrichment result ####
  ### input is list object ###
  num_len=length(enrich_list)
  go_path=paste0(paste(filename,"GO_out","k",num_len,sep = "_"),".csv")
  go_csv=dplyr::bind_rows(enrich_list, .id = "cluster_label")
  write.csv(go_csv,file=go_path)
  return(go_csv)
}


#' Title
#'
#' @param marker_lab
#' @param enrich_list
#' @param top_go_name
#'
#' @return GO sescription
#' @export
#'
#' @examples add_go_heat( )
add_go_heat<-function(marker_lab,enrich_list,top_go_name=1){
  ## this function is used to produce stuff that will be used in Heatmap function ##
  ## for row split and row name. ####
  ## top_go_name is number of GO name ranging 1 to 5.
  ## If fdr p value are quite similar and small (e-8 or e-7), you may be interested in other GO name ##
  ## output is for option "row_split  ###
  c=length(enrich_list)
  tn=top_go_name
  g_in=marker_lab
  go_result=enrich_list
  for (i in 1:c){
    g_in[g_in==i]=go_result[[i]]$Description[tn]
  }
  return(g_in)
}

#' Title
#'
#' @param go_title
#'
#' @return GO term
#' @export
#'
#' @examples SaveGo( )
SaveGo <- function(go_title){
  # Saves Go terms
  tosaveGo = do.call(cbind.data.frame, as.list(go_title))
  tosaveGoDF =tidyr::gather(tosaveGo, Gene, GOTerm, factor_key=TRUE) %>%dplyr::arrange(factor(GOTerm))
  return(tosaveGoDF)
}


#' Title
#'
#' @param ngenes
#' @param g_csv
#' @param db
#' @param Gomet
#'
#' @return Gene list
#' @export
#'
#' @examples HLGenes( )
HLGenes <-function(ngenes=5,g_csv,db,Gomet){
  # Highlights genes randomly
  geneses = paste0("gene",1:ngenes)
  genes = g_csv %>% dplyr::group_by(cluster_label) %>% dplyr::filter(FDR_Pvalue == min(FDR_Pvalue)) %>%
    dplyr::ungroup() %>%tidyr::separate(geneID,geneses,"/")
  geneVec= genes %>% dplyr::select(-cluster_label,-Description,-FDR_Pvalue) %>% unique() %>%
    as.matrix() %>% c() %>% na.omit() %>% as.character()
  if (Gomet %in% c("Wiki","Kegg","Custom")) {
    geneVec =clusterProfiler::bitr(gene=geneVec,fromType="ENTREZID",toType = c("SYMBOL"),
                   OrgDb =db, drop = T)
    geneVec = geneVec$SYMBOL
  }
  geneVecIdx = which(rownames(scaled_mat) %in% geneVec)
  return(geneVecIdx)
}



#' Title
#'
#' @param gsea_result
#'
#' @return GSEA Description
#' @export
#'
#' @examples CleanName( )
CleanName <- function(gsea_result) {
  # Clean GSEA name
  names(gsea_result@geneSets) = str_replace_all(str_remove(names(gsea_result@geneSets),"HALLMARK_"),"_"," ")
  gsea_result@result$ID = str_replace_all(str_remove(gsea_result@result$ID,"HALLMARK_"),"_"," ")
  gsea_result@result$Description = str_replace_all(str_remove(gsea_result@result$Description,"HALLMARK_"),"_"," ")
  return(gsea_result)
}
