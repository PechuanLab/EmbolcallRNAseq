#' Calculate projection score from a Seurat object
#'
#' Extracts an assay matrix from a Seurat object and passes it to
#' \code{ProjectionScore}. By default, the function uses the active assay and the
#' normalized \code{data} layer/slot. If variable features are present, they are
#' used by default to keep the matrix focused on informative genes.
#'
#' @param seu Seurat object.
#' @param NPCs Integer. Number of principal components to use. Default is 3.
#' @param nboot Integer. Number of bootstrap iterations. Default is 100.
#' @param assay Character. Assay to use. Defaults to the object's active assay.
#' @param layer Character. Seurat v5 layer to extract. If \code{NULL}, the
#'   function tries \code{slot} as a layer first and then as a Seurat v4 slot.
#' @param slot Character. Seurat v4 slot fallback. Default is \code{"data"}.
#' @param features Character vector of features to use. If supplied, this takes
#'   precedence over \code{use.variable.features}.
#' @param use.variable.features Logical. Use \code{VariableFeatures(seu)} when
#'   \code{features} is not supplied and variable features are available.
#'   Default is \code{TRUE}.
#' @param strict.features Logical. If \code{TRUE}, error when any requested
#'   feature is absent from the assay matrix. If \code{FALSE}, absent features
#'   are dropped with a warning. Default is \code{FALSE}.
#'
#' @return Numeric projection score.
#' @export
#'
#' @examples
#' \dontrun{
#' score <- SeuratProjectionScore(seu, assay = "RNA", layer = "data")
#' score_all_features <- SeuratProjectionScore(seu, use.variable.features = FALSE)
#' }
SeuratProjectionScore <- function(seu,
                                  NPCs = 3,
                                  nboot = 100,
                                  assay = NULL,
                                  layer = NULL,
                                  slot = "data",
                                  features = NULL,
                                  use.variable.features = TRUE,
                                  strict.features = FALSE) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("The SeuratObject package is required to use SeuratProjectionScore.", call. = FALSE)
  }

  if (!inherits(seu, "Seurat")) {
    stop("`seu` must be a Seurat object.", call. = FALSE)
  }

  if (!is.numeric(NPCs) || length(NPCs) != 1 || NPCs < 1) {
    stop("`NPCs` must be a positive integer.", call. = FALSE)
  }

  if (!is.numeric(nboot) || length(nboot) != 1 || nboot < 1) {
    stop("`nboot` must be a positive integer.", call. = FALSE)
  }

  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(seu)
  }

  assay_matrix <- .SeuratProjectionScoreGetAssayData(
    seu = seu,
    assay = assay,
    layer = layer,
    slot = slot
  )

  if (is.null(features) && isTRUE(use.variable.features)) {
    features <- tryCatch(
      SeuratObject::VariableFeatures(seu, assay = assay),
      error = function(e) character()
    )
    if (length(features) == 0) {
      features <- NULL
    }
  }

  if (!is.null(features)) {
    features <- unique(as.character(features))
    available_features <- rownames(assay_matrix)
    missing_features <- setdiff(features, available_features)

    if (length(missing_features) > 0) {
      msg <- paste0(
        length(missing_features),
        " requested feature(s) were not found in assay `",
        assay,
        "`: ",
        paste(utils::head(missing_features, 10), collapse = ", "),
        if (length(missing_features) > 10) ", ..." else ""
      )
      if (isTRUE(strict.features)) {
        stop(msg, call. = FALSE)
      }
      warning(msg, call. = FALSE)
    }

    features <- intersect(features, available_features)
    if (length(features) == 0) {
      stop("No requested features were found in the selected assay matrix.", call. = FALSE)
    }

    assay_matrix <- assay_matrix[features, , drop = FALSE]
  }

  data_matrix <- as.matrix(assay_matrix)
  storage.mode(data_matrix) <- "double"

  if (nrow(data_matrix) == 0 || ncol(data_matrix) == 0) {
    stop("The selected assay matrix must have at least one feature and one cell.", call. = FALSE)
  }

  if (!all(is.finite(data_matrix))) {
    stop("The selected assay matrix contains non-finite values.", call. = FALSE)
  }

  ProjectionScore(DataMatrix = data_matrix, NPCs = as.integer(NPCs), nboot = as.integer(nboot))
}

.SeuratProjectionScoreGetAssayData <- function(seu, assay, layer, slot) {
  if (!is.null(layer)) {
    return(tryCatch(
      SeuratObject::GetAssayData(object = seu, assay = assay, layer = layer),
      error = function(layer_error) {
        tryCatch(
          SeuratObject::GetAssayData(object = seu, assay = assay, slot = layer),
          error = function(slot_error) {
            stop(
              "Could not extract layer/slot `", layer, "` from assay `", assay, "`. ",
              "Layer error: ", conditionMessage(layer_error), "; ",
              "slot error: ", conditionMessage(slot_error),
              call. = FALSE
            )
          }
        )
      }
    ))
  }

  tryCatch(
    SeuratObject::GetAssayData(object = seu, assay = assay, layer = slot),
    error = function(layer_error) {
      tryCatch(
        SeuratObject::GetAssayData(object = seu, assay = assay, slot = slot),
        error = function(slot_error) {
          stop(
            "Could not extract layer/slot `", slot, "` from assay `", assay, "`. ",
            "Layer error: ", conditionMessage(layer_error), "; ",
            "slot error: ", conditionMessage(slot_error),
            call. = FALSE
          )
        }
      )
    }
  )
}
