
#' Title
#'
#' @param ExprMat
#' @param pure_markers
#'
#' @return PC score matrix
#' @export
#'
#' @examples PCscoreMat( )

PCscoreMat <- function(ExprMat,pure_markers){

  sam_col=ncol(ExprMat)
  cluster_row=pure_markers$cluster %>% dplyr::n_distinct()
  pc_df = matrix(nrow=cluster_row,ncol = sam_col)
  cluster_name=pure_markers%>% dplyr::group_by(cluster) %>% dplyr::summarise()
  cluster_name=cluster_name$cluster
  rownames(pc_df)=cluster_name
  colnames(pc_df)=colnames(ExprMat)
  cluster_no=c()

  # More preparation
  for (i in 1:cluster_row){
    # i=1
    marker=pure_markers %>% dplyr::filter(cluster==cluster_name[i]) %>% dplyr::select(gene)
    marker=marker$gene

    gene_in=intersect(rownames(ExprMat),marker)

    if(length(gene_in)!=0){
      if(length(gene_in) == 1) {
        marker_mat=t(as.matrix(ExprMat[gene_in,]))
      }else {
        marker_mat=ExprMat[gene_in,]
      }
      marker_score=gsScore(gm=marker_mat)
      marker_score=scale(marker_score)
      marker_score=as.vector(marker_score)
      pc_df[cluster_name[i],]=marker_score
    }else{
      cluster_no=c(cluster_no,cluster_name[i])
    }

  }

  row_keep=setdiff(cluster_name,cluster_no)
  pc_matrix=pc_df[row_keep,]
  return(pc_matrix)
}

