#' Title
#'
#' @param Signature
#' @param Pheno
#' @param my_comparisons
#' @param Factor
#' @param Palette
#' @return PC score plot
#' @export
#'
#' @examples GeneScorePlot(Signature = pcSig, Pheno = yf$samples,  my_comparisons = list(c("TumorDay1","SkinDay8") ,c("TumorDay3","SkinDay8"),c("TumorDay8","SkinDay8"),c("PoreDay1","SkinDay8"),c("PoreDay3","SkinDay8"),c("PoreDay8","SkinDay8")))

GeneScorePlot <- function(Signature,Pheno,my_comparisons,Factor, Palette) {
  # Add to the phenotype information
  Pheno$Score = Signature
  # Plot
  ggboxplot(Pheno, x = Factor, y = "Score",
            color = Factor, palette = Palette ,add="jitter")+
    stat_compare_means(comparisons = my_comparisons,method ="t.test",
                       method.args = list(alternative = "two.sided"),
                       p.adjust.method = "fdr",var.equal=F,
                       label="p.signif")+
    ylab("Gene Score")+
    xlab("")+
    theme(axis.text.x = element_text(angle = 90,hjust = 1),legend.position = "none")+
    font("x.text", size = 11)+font("y.text", size = 11)

}
