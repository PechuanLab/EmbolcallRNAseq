#' Title
#'
#' @param expmat
#' @param nk
#' @param scale
#' @param organism
#' @param top_go
#' @param Method
#' @param cust_db
#'
#' @return Objects to GOHeatmap
#' @export
#'
#' @examples enrich_out()
enrich_out<-function(expmat,nk,scale=TRUE,organism,
                     top_go=3,
                     Method = c("GO","Wiki","Kegg","Custom"),
                     cust_db=NULL){
  ##   this function is used to annotate row with GO name   ###
  ## Output of this function is a list including 2 object   ###
  ## parameter instruction: top_go : you can choose how many top least p value as well as other infor. ###
  go_result=vector(mode = "list",length=nk)
  if(isTRUE(scale)){
    mat = t(scale(t(expmat)))
  }else{mat=expmat}

  if(organism=="Mus musculus"){
    db=org.Mm.eg.db::org.Mm.eg.db
    orgKegg = "mmu"
  }else if(organism=="Homo sapiens")
  {
    db=org.Hs.eg.db::org.Hs.eg.db
    orgKegg = "hsa"}

  group = cluster::pam(mat, k = nk)
  g_in=group$clustering

  for (i in 1:nk){
    gene_vec=names(g_in)[g_in ==i]
    gene_cluster<-clusterProfiler::bitr(gene=gene_vec,fromType="SYMBOL",toType = c("ENTREZID","ENSEMBL"),
                       OrgDb =db, drop = T)

    if(Method=="GO"){
      over_rep=clusterProfiler::enrichGO(gene=gene_cluster$ENTREZID,OrgDb = db,
                        ont = "ALL",pAdjustMethod = "fdr",
                        pvalueCutoff = 0.01,qvalueCutoff = 0.05,
                        readable = TRUE)
      enrich_out=data.frame(Description=over_rep@result$Description[1:top_go],
                            FDR_Pvalue=over_rep@result$p.adjust[1:top_go],
                            geneID=over_rep@result$geneID[1:top_go])
      go_result[[i]]=enrich_out
    }else if(Method=="Wiki"){
      over_rep=clusterProfiler::enrichWP(gene=gene_cluster$ENTREZID,
                        organism = organism,
                        pAdjustMethod = "fdr",
                        pvalueCutoff = 0.01,qvalueCutoff = 0.05
      )
      enrich_out=data.frame(Description=over_rep@result$Description[1:top_go],
                            FDR_Pvalue=over_rep@result$p.adjust[1:top_go],
                            geneID=over_rep@result$geneID[1:top_go])
      go_result[[i]]=enrich_out
    }else if(Method=="Kegg"){

      over_rep=clusterProfiler::enrichKEGG(gene=gene_cluster$ENTREZID,
                          organism = orgKegg,
                          pAdjustMethod = "fdr",
                          pvalueCutoff = 0.01,qvalueCutoff = 0.05
      )
      enrich_out=data.frame(Description=over_rep@result$Description[1:top_go],
                            FDR_Pvalue=over_rep@result$p.adjust[1:top_go],
                            geneID=over_rep@result$geneID[1:top_go])
      go_result[[i]]=enrich_out
    } else if(Method=="Custom"){

      over_rep=clusterProfiler::enricher(gene_cluster$ENTREZID,pvalueCutoff = 0.1,minGSSize =1,TERM2GENE=cust_db)

      enrich_out=data.frame(Description=over_rep@result$Description[1],
                            FDR_Pvalue=over_rep@result$p.adjust[1],
                            geneID=over_rep@result$geneID[1])
      go_result[[i]]=enrich_out


    }
  }
  out_list=list("go_result"=go_result,"marker_label"=group$clustering)

  return(out_list)
}
