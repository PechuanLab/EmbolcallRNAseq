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
    CellcylcePal = c("#ff595e", "#1982C4", "#8AC926", "#FFCA3A"),
    Big30 = c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
  "#3182bd", "#e6550d", "#31a354", "#756bb1", "#636363",
  "#9ecae1", "#fdae6b", "#a1d99b", "#bcbddc", "#969696",
  "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff9896"
)
  )

  if (!palette %in% names(palettes)) {
    stop(sprintf("Invalid palette name. Choose from: %s", paste(names(palettes), collapse = ", ")))
  }

  return(palettes[[palette]])
}
