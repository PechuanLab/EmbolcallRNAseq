#' Compute and plot PROGENy pathway activities for Seurat clusters
#'
#' Runs PROGENy on a Seurat object, adds the pathway activity assay, scales it,
#' summarizes activities by cluster (or provided identity column), and saves a
#' heatmap of averaged pathway activities per cluster.
#'
#' @param seu A Seurat object.
#' @param id_col Optional character. Column name in Seurat metadata to use as cell
#'   identities (if provided, \code{Idents(seu)} will be set to this column).
#' @param organism Character. "Human" or "Mouse" (passed to PROGENy). Default "Human".
#' @param top Integer. Number of top target genes per pathway for PROGENy. Default 100.
#' @param perm Integer. Number of permutations for PROGENy. Default 1.
#' @param return_assay Logical. If TRUE, PROGENy results are returned as a new assay.
#'   Default TRUE.
#' @param assay_name Character. Name to use for the PROGENy assay. Default "progeny".
#' @param scale_assay Logical. If TRUE, scales the PROGENy assay using Seurat::ScaleData. Default TRUE.
#' @param outfile Character. Path to output PDF heatmap file. Default "Progeny_pathway_by_detailedID.pdf".
#' @param width Numeric. Width (in inches) of the output PDF. Default 12.
#' @param top_bar Optional. Top annotation object for ComplexHeatmap::Heatmap (passed through). Default NULL.
#'
#' @return A list with elements \code{seu} (the Seurat object with added assay) and
#'   \code{summary} (data.frame with average and sd pathway activities per cluster).
#' @export
#'
#' @examples
#' \dontrun{
#'   out <- ProgenySeurat(seu, id_col = "ClusterFine", organism = "Human",
#'                       top = 100, outfile = "Progeny.pdf")
#'   head(out$summary)
#' }
ProgenySeurat <- function(seu,
                          id_col = NULL,
                          organism = c("Human", "Mouse"),
                          top = 100,
                          perm = 1,
                          return_assay = TRUE,
                          assay_name = "progeny",
                          scale_assay = TRUE,
                          outfile = "Progeny_pathway_by_detailedID.pdf",
                          width = 12,
                          top_bar = NULL) {
  organism <- match.arg(organism)

  # Optionally set identities from metadata column
  if (!is.null(id_col)) {
    if (!id_col %in% colnames(seu@meta.data)) {
      stop("id_col not found in Seurat object's meta.data")
    }
    Seurat::Idents(seu) <- seu@meta.data[[id_col]]
  }

  # Compute PROGENy pathway activities and add as assay
  seu <- progeny::progeny(seu,
                         scale = FALSE,
                         organism = organism,
                         top = top,
                         perm = perm,
                         return_assay = return_assay,
                         assay_name = assay_name)

  # Scale PROGENy assay if requested
  if (scale_assay) {
    seu <- Seurat::ScaleData(seu, assay = assay_name)
  }

  # Prepare cluster mapping
  CellsClusters <- data.frame(
    Cell = names(Seurat::Idents(seu)),
    CellType = as.character(Seurat::Idents(seu)),
    stringsAsFactors = FALSE
  )

  # Extract scaled PROGENy scores (cells x pathways)
  progeny_mat <- as.data.frame(t(Seurat::GetAssayData(seu, slot = "scale.data", assay = assay_name)))
  progeny_scores_df <- progeny_mat %>%
    tibble::rownames_to_column("Cell") %>%
    tidyr::pivot_longer(-Cell, names_to = "Pathway", values_to = "Activity")

  # Join with cluster assignments and summarize
  progeny_scores_df <- dplyr::inner_join(progeny_scores_df, CellsClusters, by = "Cell")
  summarized_progeny_scores <- progeny_scores_df %>%
    dplyr::group_by(Pathway, CellType) %>%
    dplyr::summarise(avg = mean(Activity, na.rm = TRUE),
                     std = sd(Activity, na.rm = TRUE),
                     .groups = "drop")

  # Prepare matrix for heatmap: pathways as rows, clusters as columns
  summarized_progeny_scores_df <- summarized_progeny_scores %>%
    dplyr::select(-std) %>%
    tidyr::pivot_wider(names_from = CellType, values_from = avg) %>%
    tibble::column_to_rownames("Pathway") %>%
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

  progeny_mat_plot <- as.matrix(summarized_progeny_scores_df) %>% t()

  # Color mapping
  paletteLength <- 100
  myColor <- grDevices::colorRampPalette(c("steelblue1", "white", "tomato"))(paletteLength)

  progenyBreaks <- c(
    seq(min(progeny_mat_plot, na.rm = TRUE), 0, length.out = ceiling(paletteLength / 2) + 1),
    seq(max(progeny_mat_plot, na.rm = TRUE) / paletteLength,
        max(progeny_mat_plot, na.rm = TRUE),
        length.out = floor(paletteLength / 2))
  )

  # Draw heatmap to PDF
  grDevices::pdf(outfile, width = width)
  ComplexHeatmap::Heatmap(progeny_mat_plot,
                          cluster_columns = TRUE,
                          cluster_rows = TRUE,
                          column_names_gp = grid::gpar(fontsize = 10),
                          row_names_gp = grid::gpar(fontsize = 10),
                          column_split = min(3, ncol(progeny_mat_plot)),
                          row_title = "cluster_%s",
                          row_title_gp = grid::gpar(fontsize = 0),
                          row_split = min(7, nrow(progeny_mat_plot)),
                          top_annotation = top_bar,
                          border_gp = grid::gpar(col = "gray", lty = 2),
                          row_dend_gp = grid::gpar(col = "gray"),
                          name = "Pathway Activities",
                          col = myColor,
                          heatmap_legend_param = list(at = c(min(progeny_mat_plot, na.rm = TRUE), 0, max(progeny_mat_plot, na.rm = TRUE))))
  grDevices::dev.off()

  return(list(seu = seu, summary = summarized_progeny_scores))
}
