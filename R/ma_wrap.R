
#' Generate and save an MA plot for differential expression results
#'
#' This function creates an MA plot (log2 fold change vs. mean expression) from a differential expression signature table and saves it as a PNG file. It uses ggpubr::ggmaplot for visualization, highlights significant genes, and handles infinite/NA values.
#'
#' @param signature Data frame containing differential expression results. Must include columns: `AveExpr` (average expression), `logFC` (log2 fold change), `adj.P.Val` (adjusted p-value), and `symbol` (gene names).
#' @param contrast Character string. Name of the contrast/condition, used in the output filename.
#' @param Prefix Character string. Prefix for the output filename.
#'
#' @return Invisibly returns NULL. The function is called for its side effect of saving a PNG MA plot to disk.
#' @export
#'
#' @examples
#' \dontrun{
#'   ma_wrap(signature, contrast = "Treatment_vs_Control", Prefix = "Experiment1_")
#' }

ma_wrap<-function(signature,contrast,Prefix){
  basemean = signature$AveExpr
  fold = signature$logFC
  p_val_adj = signature$adj.P.Val

  dat = data.frame(name=signature$symbol,baseMean =basemean, log2FoldChange = fold, padj = p_val_adj)
  rownames(dat)=signature$symbol

  dat$baseMean=2**(dat$baseMean)

  dat=do.call(data.frame,lapply(dat, function(x) replace(x, is.infinite(x),NA)))

  dat=dat %>% tibble::as_tibble()
  dat = dat %>% tidyr::drop_na()

  dat=as.data.frame(dat)
  png(paste(Prefix,contrast,"MA_plot.png",sep=""), units="in", width=5, height=5, res=300)
  print(
    ggpubr::ggmaplot(dat,
             fdr = 0.05, fc = 2, size = 0.4,
             palette = c("#B31B21", "#1465AC", "darkgray"),
             genenames = as.vector(dat$name),
             select.top.method = c( "padj","fc"),
             legend = "top", top = 25,
             font.label = c("bold", 8), label.rectangle = TRUE,
             font.legend = "bold",
             font.main = "bold",
             ggtheme = ggplot2::theme_classic())

  )
  dev.off()
}
