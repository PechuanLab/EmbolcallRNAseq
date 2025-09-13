
##' Plot STAMP deconvolution signature scores by compartment
##'
##' Generates a boxplot of STAMP atlas signature scores for different cell types within a specified compartment, using principal component scores and phenotype information. Useful for visualizing deconvolution results across experimental groups.
##'
##' @param ExprMat Matrix. Expression matrix (genes x samples) for deconvolution.
##' @param Pheno Data frame. Phenotype/sample metadata to merge with scores.
##' @param my_comparisons List. List of group comparisons for plotting (not directly used in function, but included for compatibility).
##' @param Factor Character. Column name in `Pheno` to use for coloring groups.
##' @param Palette Character vector. Colors to use for group levels in the plot.
##' @param Compartment Character. Compartment to subset cell types for plotting (e.g., "Lymphoid", "Stroma", "Myeloid", "Other").
##'
##' @return A ggplot2 object: boxplot of STAMP atlas signature scores by cell type and group.
##' @export
##'
##' @examples
##' # Example usage (requires proper input objects):
##' # STAMPDeconvolPlot(ExprMat, Pheno, my_comparisons, Factor = "Group", Palette = c("red", "blue"), Compartment = "Lymphoid")

STAMPDeconvolPlot <- function(ExprMat,Pheno,my_comparisons,Factor, Palette, Compartment) {
  # Add to the phenotype information
  pc_matrix = PCscoreMat(ExprMat,pure_markers)
  df_atlas = t(pc_matrix) %>% as.data.frame()
  df_atlas = cbind(df_atlas,Pheno) 
  df_atlas = pivot_longer(df_atlas,cols = Cd4:moDC, names_to = "CellType",values_to = "Score")
  # subset to compartment
  if (Compartment == "Lymphoid") {
    df_atlas = df_atlas %>% filter(CellType %in% c("Cd4","Cd8","Tregs","Nk:Ncr1","gdTCell","ILC","NkTcells"))
  } else  if (Compartment == "Stroma") {
    df_atlas = df_atlas %>% filter(CellType %in% c("KPP","Keratinocytes","Chondrocyte","Schwann","Platelet","Melanocytes"))
  }else  if (Compartment == "Myeloid") {
    df_atlas = df_atlas %>% filter(CellType %in% c("Fibroblasts","VEC","LEC","Pericytes"))
  }else if (Compartment == "Other") {
    df_atlas = df_atlas %>% filter(CellType %in% c("cDC1","moDC","Monocyte","Monocyte",
                                                   "Langerhans","DC2","PlasmocytoidDC","Neutrophils","Macrophage"))
  }
  # Plot
  ggboxplot(df_atlas, x = "CellType", y = "Score",
            color = Factor, palette = Palette ,add="jitter")+
    ylab("STAMP Atlas Signature Score")+
    xlab(Compartment)+
    theme(axis.text.x = element_text(angle = 90,hjust = 1),legend.position = "top")+
    font("x.text", size = 11)+font("y.text", size = 11)

}

