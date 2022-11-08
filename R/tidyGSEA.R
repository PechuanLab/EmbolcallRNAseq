#' Title
#'
#' @param fgseaRes
#' @param npaths
#'
#' @return NES table
#' @export
#'
#' @examples tidyGSEA( )

tidyGSEA <- function(fgseaRes,npaths = 10 ) {
  topUP =  fgseaRes %>%
    dplyr::filter(NES>0) %>% dplyr::top_n(p.adjust/abs(NES),n = npaths) %>%
    tibble::as_tibble()
  topDown =  fgseaRes %>%
    dplyr::filter(NES<0) %>% dplyr::top_n(p.adjust/abs(NES),n = npaths) %>%
    tibble::as_tibble()
  fgseaResTidy = rbind(topUP,topDown)
  fgseaResTidy = fgseaResTidy %>%  dplyr::arrange(desc(NES))
  return(fgseaResTidy)
}
