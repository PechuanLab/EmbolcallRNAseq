#' Get a Color Palette
#'
#' This function returns a specified color palette.
#'
#' @param palette A character string specifying the palette to use. Options are "palette1", "palette2", and "palette3".
#' @return A vector of color hex codes.
#' @examples
#' ColorPalettes("palette1")
#' ColorPalettes("palette2")
ColorPalettes <- function(palette = "palette1") {
  palettes <- list(
    "palette1" = c("#FF5733", "#33FF57", "#3357FF"),
    "palette2" = c("#F0E442", "#0072B2", "#D55E00"),
    "palette3" = c("#CC79A7", "#56B4E9", "#009E73")
  )
  
  if (!palette %in% names(palettes)) {
    stop("Invalid palette name. Choose from 'palette1', 'palette2', or 'palette3'.")
  }
  
  return(palettes[[palette]])
}
