
#' Title
#'
#' @param gensmbl
#' @param x.factor
#' @param trace.factor
#' @param ExprMat
#' @param cpalette
#'
#' @return save ggplot
#' @export
#'
#' @examples InteractionPlotIrtx( )


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
