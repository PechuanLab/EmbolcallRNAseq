#' Title
#'
#' @param Signature
#' @param Pheno
#' @param my_comparisons
#' @param Factor
#' @param Palette
#' @param Compartment
#' @return PC score plot
#' @export
#'
#' @examples GeneScorePlot(Signature = pcSig, Pheno = yf$samples,  my_comparisons = list(c("TumorDay1","SkinDay8") ,c("TumorDay3","SkinDay8"),c("TumorDay8","SkinDay8"),c("PoreDay1","SkinDay8"),c("PoreDay3","SkinDay8"),c("PoreDay8","SkinDay8")))

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

