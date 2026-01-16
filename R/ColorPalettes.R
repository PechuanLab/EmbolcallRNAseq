#' Get a Color Palette
#'
#' Returns one of the predefined color palettes used across the package.
#'
#' @param palette Character. Name of the palette to return. Options:
#'   "TreatmentPalette", "ClusterPaletteCoarse", "CellcylcePal".
#'   Default is "TreatmentPalette".
#' @return A character vector of color hex codes.
#' @export
#' @examples
#' ColorPalettes("TreatmentPalette")
#' ColorPalettes("ClusterPaletteCoarse")
#' ColorPalettes("CellcylcePal")
ColorPalettes <- function(palette = "TreatmentPalette") {
  palettes <- list(
    TreatmentPalette = c("gray", "#8AC926", "forestgreen"),
    ClusterPaletteCoarse = c("sienna4", "firebrick", "gray47", "royalblue", "skyblue",
                             "#FF7F00", "gold", "#6B8E23", "salmon", "tan1"),
    CellcylcePal = c("#ff595e", "#1982C4", "#8AC926", "#FFCA3A")
  )

  if (!palette %in% names(palettes)) {
    stop(sprintf("Invalid palette name. Choose from: %s", paste(names(palettes), collapse = ", ")))
  }

  return(palettes[[palette]])
}
