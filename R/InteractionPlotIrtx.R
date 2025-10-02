#' Plot gene expression interaction over time and group
#'
#' This function creates an interaction plot for a single gene, showing its expression over time and across groups (clones), with error bars and group-specific normalization. The plot is saved as a PDF.
#'
#' @param gensmbl Character. Gene symbol to plot.
#' @param x.factor Numeric or factor. Time or x-axis variable for the plot.
#' @param trace.factor Factor. Grouping variable (e.g., clone or condition).
#' @param ExprMat Numeric matrix. Expression matrix with genes as rows and samples as columns.
#' @param cpalette Character vector. Colors for each group/clone.
#'
#' @return Invisibly returns NULL. Side effect: saves a PDF plot named "<gene>_Interaction.pdf".
#' @export
#'
#' @examples
#' \dontrun{
#'   InteractionPlotIrtx("GeneA", x.factor = time, trace.factor = group, ExprMat = exprs, cpalette = c("red", "blue"))
#' }


InteractionPlotIrtx <-function(gensmbl,x.factor,trace.factor,ExprMat,cpalette) {
  # Interaction Plot for Iratxe's style
  # Filename
  filename = paste(gensmbl,"_","Interaction.pdf",sep="")
  # Get the Gene
  x = ExprMat[gensmbl,]
  # Prepare for the plot
  df=data.frame(cpm=x,time = x.factor,clone=trace.factor)
  aggregate = df %>% dplyr::group_by(time,clone) %>% dplyr::summarise(Mean = mean(cpm),SE = sd(cpm)/sqrt(n())) %>% dplyr:: ungroup()
  normaggregate = aggregate %>% dplyr:: group_by(clone) %>%
    dplyr::mutate(Normal = Mean-Mean[time==0][1L]) %>% dplyr::ungroup()
  #  Plot
  pdf(filename)
  print(
    ggpubr::ggline(normaggregate,"time","Normal",
           color = "clone", palette = gplots::col2hex(ClonePalette),numeric.x.axis=T)+
      xlim(c(0,72))+ylab("Centered Expression")+xlab("Time (Hours)")+
      geom_errorbar(aes(x = time, ymin = Normal-SE, ymax = Normal+SE), width=.1,
                    position=position_dodge(0.05))+ggtitle(gensmbl)
  )
  dev.off()
}
