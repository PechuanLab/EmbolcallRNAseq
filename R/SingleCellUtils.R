#' Subsets the neutropgils based on known markers (only mouse)
#' Takes a data matrix and resamples it to obtain a distribution for the fraction of variance explained by a given number of principal components.
#' @param DataMatrix
#' @param NPCs
#' @param nboot
#'
#' @return Fraction of variance for projection score bootstrapped average
#' @export
#'
#' @examples AlphaPE()
#' 
QCFilter <- function(DataMatrix,NPCs =  3, nboot = 100){
  # Extreact S100a8 and 
  nFeathLow = 500
nCountHigh = 200000
riboT = 40
mitoT = 15

Idents(seu) = seu@meta.data$SampleID
# Inspect the thresholds 
p1 = VlnPlot_scCustom(seu,"nFeature_RNA", pt.size = 0, log = T,
                      plot_boxplot =T) +
  geom_hline(yintercept = nFeathLow)+ggtitle(paste0("nFeature_RNA=",nFeathLow))
p1
p2 = VlnPlot_scCustom(seu,"nFeature_RNA", pt.size = 0, log = F,
                      plot_boxplot =T) +
  geom_hline(yintercept = nFeathLow)+ggtitle(paste0("nFeature_RNA=",nFeathLow))
p2
p3 = VlnPlot_scCustom(seu,"nCount_RNA", pt.size = 0, log = T,
                      plot_boxplot =T) +
  geom_hline(yintercept = nCountHigh)+ggtitle(paste0("nCount_RNA=",nCountHigh))
p3
p4 = VlnPlot_scCustom(seu,"nCount_RNA", pt.size = 0, log = F,
                      plot_boxplot =T) +
  geom_hline(yintercept = nCountHigh)+ggtitle(paste0("nCount_RNA=",nCountHigh))
p4

p5 =VlnPlot_scCustom(seu,"percent_ribo", pt.size = 0, log = F,  
                     plot_boxplot =T) +
  geom_hline(yintercept = riboT)+ggtitle(paste0("percent_ribo=",riboT))
p5
p6 = VlnPlot_scCustom(seu ,"percent_mito", pt.size = 0, log =F,
                      plot_boxplot =T) +
  geom_hline(yintercept = mitoT)+ggtitle(paste0("percent_mito=",mitoT))
p6

# save the plots
pdf(paste0(SampleName,"QCRest.pdf"), width = 15, height = 12)
plot_grid(p1,p2,p3,p4,p5,p6)
dev.off()

}
