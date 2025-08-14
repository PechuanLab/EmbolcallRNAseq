
#' Compute a principal component score matrix for clusters
#'
#' This function calculates a score matrix for each cluster based on expression of marker genes.
#' For each cluster, it extracts the relevant genes from the expression matrix and computes a score
#' using the `gsScoreLegacy()` function. The scores are scaled and returned as a matrix.
#'
#' @param ExprMat A numeric matrix of gene expression values with genes as rows and samples as columns.
#' @param pure_markers A data frame with columns `gene` and `cluster`, specifying marker genes for each cluster.
#'
#' @return A matrix of scaled scores with clusters as rows and samples as columns.
#'         Clusters with no matching genes in the expression matrix are excluded.
#' @export
#'
#' @examples
#' # PCscoreMat(expression_matrix, marker_df)

PCscoreMat <- function(ExprMat, pure_markers) {
  sam_col <- ncol(ExprMat)
  cluster_row <- pure_markers$cluster %>% dplyr::n_distinct()
  
  pc_df <- matrix(nrow = cluster_row, ncol = sam_col)
  cluster_name <- pure_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise() %>%
    dplyr::pull(cluster)
  
  rownames(pc_df) <- cluster_name
  colnames(pc_df) <- colnames(ExprMat)
  
  cluster_no <- c()
  
  for (i in seq_len(cluster_row)) {
    marker <- pure_markers %>%
      dplyr::filter(cluster == cluster_name[i]) %>%
      dplyr::pull(gene)
    
    gene_in <- intersect(rownames(ExprMat), marker)
    
    if (length(gene_in) != 0) {
      if (length(gene_in) == 1) {
        marker_mat <- t(as.matrix(ExprMat[gene_in, ]))
      } else {
        marker_mat <- ExprMat[gene_in, ]
      }
      
      marker_score <- gsScoreLegacy(gm = marker_mat)
      marker_score <- scale(marker_score)
      pc_df[cluster_name[i], ] <- as.vector(marker_score)
    } else {
      cluster_no <- c(cluster_no, cluster_name[i])
    }
  }
  
  row_keep <- setdiff(cluster_name, cluster_no)
  pc_matrix <- pc_df[row_keep, ]
  
  return(pc_matrix)
}

