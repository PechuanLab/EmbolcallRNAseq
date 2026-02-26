#' Annotate Seurat object with SingleR using a reference dataset
#'
#' Converts Seurat -> SingleCellExperiment, loads/accepts a reference (SCE/Seurat/path),
#' runs SingleR and adds predicted labels to Seurat metadata.
#'
#' @param seu Seurat object (test).
#' @param reference Either a SingleCellExperiment, Seurat object, or path to an rds file containing the reference Seurat/SCE.
#' @param reference_label_col Character. Column name in reference colData to use as labels. Default "CellType".
#' @param assay Character. Assay in the Seurat test object to convert to SCE. Default "RNA".
#' @param assay.type.ref Character. assay.type.ref passed to SingleR (e.g. "logcounts"). Default "logcounts".
#' @param output_col Character. Metadata column name to add to Seurat with predicted labels. Default "SingleR_labels".
#' @param return_singleR Logical. If TRUE return SingleR result in output. Default TRUE.
#'
#' @return If return_singleR TRUE, list(seu = Seurat object with added metadata, singleR = SingleR result).
#'         Otherwise returns the modified Seurat object.
#' @export
AnnotateWithSingleR <- function(seu,
                                reference,
                                reference_label_col = "CellType",
                                assay = "RNA",
                                assay.type.ref = "logcounts",
                                output_col = "SingleR_labels",
                                return_singleR = TRUE) {
  if (!requireNamespace("SingleR", quietly = TRUE)) {
    stop("Package 'SingleR' is required but not installed.")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required but not installed.")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed.")
  }

  # convert test Seurat to SCE
  sce_test <- tryCatch(
    as.SingleCellExperiment(seu, assay = assay),
    error = function(e) stop("Failed to convert 'seu' to SingleCellExperiment: ", e$message)
  )

  # load/convert reference
  ref_obj <- reference
  if (is.character(reference) && length(reference) == 1) {
    ref_obj <- readRDS(reference)
  }

  # convert to SCE if Seurat provided
  if (inherits(ref_obj, "Seurat")) {
    ref_sce <- tryCatch(
      as.SingleCellExperiment(ref_obj, assay = assay),
      error = function(e) stop("Failed to convert reference Seurat to SingleCellExperiment: ", e$message)
    )
  } else if (inherits(ref_obj, "SingleCellExperiment")) {
    ref_sce <- ref_obj
  } else {
    stop("Reference must be a path, Seurat object, or SingleCellExperiment.")
  }

  # check label column exists
  if (!reference_label_col %in% colnames(SingleCellExperiment::colData(ref_sce))) {
    stop("reference_label_col not found in reference colData.")
  }

  labels_vec <- SingleCellExperiment::colData(ref_sce)[, reference_label_col]

  # run SingleR
  singleR_res <- SingleR::SingleR(test = sce_test,
                                  ref = ref_sce,
                                  labels = labels_vec,
                                  assay.type.ref = assay.type.ref)

  # predicted labels
  predicted_labels <- as.character(singleR_res$labels)
  names(predicted_labels) <- rownames(singleR_res) # should match test cell names

  # add to Seurat metadata
  seu <- Seurat::AddMetaData(seu, metadata = predicted_labels, col.name = output_col)

  if (isTRUE(return_singleR)) {
    return(list(seu = seu, singleR = singleR_res))
  } else {
    return(seu)
  }
}