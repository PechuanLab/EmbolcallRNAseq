#' Plot gene signature scores by group with statistical comparisons
#'
#' This function creates a boxplot (with jitter) of gene signature scores (e.g., PC scores) grouped by a specified phenotype variable, and adds statistical comparisons between specified groups.
#'
#' @param Signature Numeric vector. Gene signature scores for each sample (length must match number of rows in \code{Pheno}).
#' @param Pheno Data frame. Phenotype/sample metadata, one row per sample.
#' @param my_comparisons List of character vectors. Each vector specifies two group names (as in \code{Factor}) to compare.
#' @param Factor Character. Column name in \code{Pheno} to use for grouping and coloring.
#' @param Palette Character vector. Colors to use for each group in the plot.
#'
#' @return A ggplot object (boxplot with jitter and statistical comparisons).
#' @export
#'
#' @examples
#' \dontrun{
#'   GeneScorePlot(
#'     Signature = pcSig,
#'     Pheno = yf$samples,
#'     my_comparisons = list(c("TumorDay1", "SkinDay8")),
#'     Factor = "Treatment",
#'     Palette = c("red", "blue")
#'   )
#' }

GeneScorePlot <- function(Signature, Pheno, my_comparisons, Factor, Palette) {
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
