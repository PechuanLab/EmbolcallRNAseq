#' Generate a Heatmap for Differentially Expressed Genes
#'
#' This function generates a heatmap for differentially expressed (DE) genes based on a specified log fold change (lfc) threshold. It uses the `treat` function to identify DE genes and scales the expression matrix for visualization. The heatmap can also include annotations based on GO terms or KEGG pathways.
#'
#' @param Contrastes A string specifying the contrast or condition being analyzed. This is used to extract DE genes and name the output files.
#' @param lfc A numeric value specifying the log fold change threshold for identifying DE genes. Default is 2.
#' @param ExprMat A numeric matrix of gene expression values, where rows correspond to genes and columns correspond to samples.
#'
#' @return A heatmap plot for DE genes is generated and saved as a PDF file. The function does not return an R object directly.
#' @export
#'
#' @examples
#' # Example usage:
#' # Assuming `Contrastes` is "ConditionA_vs_Control", `lfc` is 2, and `ExprMat` is a matrix of expression data:
#' DEHeatMap("ConditionA_vs_Control", lfc = 2, ExprMat)
#'
DEHeatMap <- function(Contrastes, lfc = 2, ExprMat) {
  # Identify DE genes using treat with the specified log fold change threshold
  FitTreat = treat(fit2, lfc = lfc)
  signature = topTreat(
    FitTreat, 
    coef = Contrastes, 
    n = Inf, 
    adjust.method = "fdr", 
    sort.by = "P", 
    p.value = 0.05
  )

  # Subset the expression matrix to include only DE genes
  topdat = ExprMat[rownames(signature),]

  # Scale the expression matrix for heatmap visualization
  scaled_mat = t(scale(t(topdat)))

  # Generate a heatmap with GO term or KEGG annotations
  g = GOHeatMap(
    scaled_mat = scaled_mat, 
    Gomet = "Wiki", 
    nk = 8,
    organismRef = "Mus musculus", 
    title = paste0("Top ", lfc, Contrastes, " Genes"),
    nHL = 5, 
    Csplit = 2, 
    top_annotation = top_bar
  )

  # Save the heatmap as a PDF
  pdf(paste0(Contrastes, "TopDE_genes.pdf"), width = 15)
  print(g)
  dev.off()
}