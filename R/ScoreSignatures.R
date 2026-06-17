#' Score curated gene signatures
#'
#' Scores directional or non-directional gene signatures on an expression matrix
#' or a Seurat object. For Seurat objects, scores can be added directly to cell
#' metadata. Directional signatures are represented as \code{up - down}.
#'
#' @param object Expression matrix/data frame with genes as rows and samples or
#'   cells as columns, or a Seurat object.
#' @param signatures Named list of signatures. Each entry can be a character
#'   vector, or a list with \code{up} and optional \code{down} character vectors.
#'   If \code{NULL}, signatures are loaded from \code{\link{PaperSignatures}}
#'   using \code{paper} and \code{species}.
#' @param method Character. Scoring backend. \code{"seurat"} uses
#'   \code{Seurat::AddModuleScore}; \code{"ucell"} uses
#'   \code{UCell::ScoreSignatures_UCell}; \code{"zscore"} averages row-z-scored
#'   genes; \code{"mean"} averages expression; \code{"pc"} uses \code{gsScore}.
#'   If \code{NULL}, Seurat objects default to \code{"seurat"} and matrices
#'   default to \code{"zscore"}.
#' @param paper Character. Paper registry key used when \code{signatures = NULL}.
#'   Default \code{"moore2026_crc_states"}.
#' @param species Character. Species used when loading paper signatures.
#'   Supported values are \code{"mouse"} and \code{"human"}.
#' @param assay Character. Seurat assay to score. Defaults to active assay.
#' @param layer Optional Seurat v5 layer to use for filtering and non-Seurat
#'   methods. For \code{method = "seurat"}, this is passed as the
#'   \code{slot} argument to \code{Seurat::AddModuleScore}.
#' @param slot Character. Seurat v4 slot/layer fallback. Defaults to
#'   \code{"counts"} for UCell and \code{"data"} otherwise.
#' @param min_genes Integer. Minimum number of retained genes required across
#'   \code{up} and \code{down}. Default \code{3}.
#' @param add_to_metadata Logical. For Seurat objects, add scores to metadata
#'   and return the modified object. Defaults to \code{TRUE} for Seurat objects.
#' @param prefix Optional character prefix for score names.
#' @param suffix Character suffix for score names. Default \code{"_score"}.
#' @param nbin,ctrl,seed,search Arguments passed to
#'   \code{Seurat::AddModuleScore}.
#' @param maxRank,ncores,BPPARAM Arguments passed to
#'   \code{UCell::ScoreSignatures_UCell}.
#' @param verbose Logical. If \code{TRUE}, report retained signature counts.
#' @param ... Additional arguments passed to the selected backend.
#'
#' @return For matrices, a signature-by-sample score matrix. For Seurat objects,
#'   returns the Seurat object with metadata columns when
#'   \code{add_to_metadata = TRUE}; otherwise returns a signature-by-cell matrix.
#'   Signature metadata and gene-overlap summaries are stored as attributes on
#'   returned matrices.
#' @export
#'
#' @examples
#' \dontrun{
#' sigs <- PaperSignatures("moore2026")
#' score_mat <- ScoreSignatures(expr_mat, sigs, method = "zscore")
#' seu <- ScoreSignatures(seu, paper = "moore2026", species = "mouse", method = "seurat")
#' seu <- ScoreSignatures(seu, sigs, method = "ucell", assay = "RNA")
#' }
ScoreSignatures <- function(object,
                            signatures = NULL,
                            method = NULL,
                            paper = "moore2026_crc_states",
                            species = c("mouse", "human"),
                            assay = NULL,
                            layer = NULL,
                            slot = NULL,
                            min_genes = 3,
                            add_to_metadata = NULL,
                            prefix = NULL,
                            suffix = "_score",
                            nbin = 24,
                            ctrl = 100,
                            seed = 1,
                            search = FALSE,
                            maxRank = 1500,
                            ncores = 1,
                            BPPARAM = NULL,
                            verbose = TRUE,
                            ...) {
  is_seurat <- inherits(object, "Seurat")
  species <- match.arg(species)

  if (is.null(method)) {
    method <- if (is_seurat) "seurat" else "zscore"
  }
  method <- .ScoreSignaturesNormalizeMethod(method)

  if (is.null(slot)) {
    slot <- if (identical(method, "ucell")) "counts" else "data"
  }
  if (is.null(add_to_metadata)) {
    add_to_metadata <- is_seurat
  }
  if (!is.numeric(min_genes) || length(min_genes) != 1 || min_genes < 1) {
    stop("`min_genes` must be a positive integer.", call. = FALSE)
  }

  if (is.null(signatures)) {
    signatures <- PaperSignatures(paper = paper, species = species)
  }

  signatures <- .ScoreSignaturesValidate(signatures)

  if (is_seurat) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("The SeuratObject package is required to score Seurat objects.", call. = FALSE)
    }
    if (is.null(assay)) {
      assay <- SeuratObject::DefaultAssay(object)
    }
    assay_matrix <- .ScoreSignaturesGetAssayData(
      seu = object,
      assay = assay,
      layer = layer,
      slot = slot
    )
  } else {
    if (!is.matrix(object) && !is.data.frame(object)) {
      stop("`object` must be an expression matrix/data frame or a Seurat object.", call. = FALSE)
    }
    assay_matrix <- as.matrix(object)
    storage.mode(assay_matrix) <- "double"
  }

  if (is.null(rownames(assay_matrix)) || anyNA(rownames(assay_matrix))) {
    stop("The expression matrix must have gene names as row names.", call. = FALSE)
  }
  if (is.null(colnames(assay_matrix))) {
    colnames(assay_matrix) <- paste0("Sample", seq_len(ncol(assay_matrix)))
  }
  if (nrow(assay_matrix) == 0 || ncol(assay_matrix) == 0) {
    stop("The selected expression matrix must have at least one gene and one sample/cell.", call. = FALSE)
  }

  filtered <- .ScoreSignaturesFilterSignatures(
    signatures = signatures,
    available_features = rownames(assay_matrix),
    min_genes = as.integer(min_genes)
  )
  signatures <- filtered$signatures

  if (isTRUE(verbose)) {
    message(
      "Scoring ", length(signatures), " signature(s) with method `", method, "`."
    )
  }

  score_matrix <- switch(
    method,
    seurat = {
      if (!is_seurat) {
        stop("`method = \"seurat\"` requires a Seurat object.", call. = FALSE)
      }
      .ScoreSignaturesScoreSeurat(
        seu = object,
        signatures = signatures,
        assay = assay,
        slot = if (!is.null(layer)) layer else slot,
        nbin = nbin,
        ctrl = ctrl,
        seed = seed,
        search = search,
        ...
      )
    },
    ucell = .ScoreSignaturesScoreUCell(
      assay_matrix = assay_matrix,
      signatures = signatures,
      maxRank = maxRank,
      ncores = ncores,
      BPPARAM = BPPARAM,
      ...
    ),
    zscore = .ScoreSignaturesScoreMatrix(
      assay_matrix = assay_matrix,
      signatures = signatures,
      method = "zscore"
    ),
    mean = .ScoreSignaturesScoreMatrix(
      assay_matrix = assay_matrix,
      signatures = signatures,
      method = "mean"
    ),
    pc = .ScoreSignaturesScoreMatrix(
      assay_matrix = assay_matrix,
      signatures = signatures,
      method = "pc"
    )
  )

  score_names <- .ScoreSignaturesScoreNames(rownames(score_matrix), prefix = prefix, suffix = suffix)
  filtered$metadata$score_name <- score_names[filtered$metadata$signature_name]
  rownames(score_matrix) <- unname(score_names[rownames(score_matrix)])

  attr(score_matrix, "signature_metadata") <- filtered$metadata
  attr(score_matrix, "signature_overlap") <- filtered$overlap

  if (is_seurat && isTRUE(add_to_metadata)) {
    score_df <- as.data.frame(t(score_matrix), check.names = FALSE)
    object <- SeuratObject::AddMetaData(object = object, metadata = score_df)
    return(object)
  }

  score_matrix
}

.ScoreSignaturesNormalizeMethod <- function(method) {
  method <- tolower(gsub("_", "", as.character(method[1])))
  if (method %in% c("seurat", "module", "addmodulescore", "addmodule")) {
    return("seurat")
  }
  if (method %in% c("ucell", "u", "umodule")) {
    return("ucell")
  }
  if (method %in% c("zscore", "z", "scaledmean")) {
    return("zscore")
  }
  if (method %in% c("mean", "average", "avg")) {
    return("mean")
  }
  if (method %in% c("pc", "pca", "pc1")) {
    return("pc")
  }
  stop("`method` must be one of `seurat`, `ucell`, `zscore`, `mean`, or `pc`.", call. = FALSE)
}

.ScoreSignaturesValidate <- function(signatures) {
  if (is.character(signatures)) {
    signatures <- list(Signature = signatures)
  }
  if (!is.list(signatures) || length(signatures) == 0) {
    stop("`signatures` must be a non-empty character vector or named list.", call. = FALSE)
  }
  if (is.null(names(signatures))) {
    names(signatures) <- paste0("Signature", seq_along(signatures))
  }
  empty_names <- names(signatures) == "" | is.na(names(signatures))
  names(signatures)[empty_names] <- paste0("Signature", which(empty_names))
  names(signatures) <- make.unique(names(signatures))

  validated <- lapply(seq_along(signatures), function(i) {
    signature <- signatures[[i]]
    signature_name <- names(signatures)[i]

    if (is.character(signature)) {
      up <- signature
      down <- character()
      metadata <- list()
    } else if (is.list(signature)) {
      up <- .ScoreSignaturesFirstNonNull(
        signature$up,
        signature$positive,
        signature$genes,
        signature$features
      )
      down <- .ScoreSignaturesFirstNonNull(signature$down, signature$negative, character())
      metadata <- signature$metadata
      if (is.null(metadata)) {
        metadata <- list()
      }
    } else {
      stop("Signature `", signature_name, "` must be a character vector or list.", call. = FALSE)
    }

    up <- unique(as.character(up))
    down <- unique(as.character(down))
    up <- up[!is.na(up) & nzchar(up)]
    down <- down[!is.na(down) & nzchar(down)]

    if (length(up) == 0 && length(down) == 0) {
      stop("Signature `", signature_name, "` has no genes.", call. = FALSE)
    }

    list(up = up, down = down, metadata = metadata)
  })
  stats::setNames(validated, names(signatures))
}

.ScoreSignaturesFirstNonNull <- function(...) {
  values <- list(...)
  for (value in values) {
    if (!is.null(value)) {
      return(value)
    }
  }
  NULL
}

.ScoreSignaturesFilterSignatures <- function(signatures, available_features, min_genes) {
  filtered <- list()
  overlap_rows <- list()
  metadata_rows <- list()

  for (signature_name in names(signatures)) {
    signature <- signatures[[signature_name]]
    up_present <- intersect(signature$up, available_features)
    down_present <- intersect(signature$down, available_features)
    retained_n <- length(up_present) + length(down_present)

    overlap_rows[[signature_name]] <- data.frame(
      signature_name = signature_name,
      original_up_n = length(signature$up),
      retained_up_n = length(up_present),
      missing_up_n = length(setdiff(signature$up, available_features)),
      original_down_n = length(signature$down),
      retained_down_n = length(down_present),
      missing_down_n = length(setdiff(signature$down, available_features)),
      stringsAsFactors = FALSE
    )

    if (retained_n < min_genes) {
      next
    }

    filtered[[signature_name]] <- list(
      up = up_present,
      down = down_present,
      metadata = signature$metadata
    )
    metadata_rows[[signature_name]] <- .ScoreSignaturesMetadataRow(
      signature_name = signature_name,
      metadata = signature$metadata
    )
  }

  if (length(filtered) == 0) {
    stop(
      "No signatures retained at least ", min_genes,
      " gene(s) present in the expression matrix.",
      call. = FALSE
    )
  }

  overlap <- do.call(rbind, overlap_rows)
  rownames(overlap) <- overlap$signature_name
  metadata <- do.call(rbind, metadata_rows)
  rownames(metadata) <- metadata$signature_name

  dropped <- length(signatures) - length(filtered)
  if (dropped > 0) {
    warning(dropped, " signature(s) dropped because too few genes were present.", call. = FALSE)
  }

  list(signatures = filtered, metadata = metadata, overlap = overlap)
}

.ScoreSignaturesMetadataRow <- function(signature_name, metadata) {
  get_chr <- function(name) {
    value <- metadata[[name]]
    if (is.null(value) || length(value) == 0) {
      return(NA_character_)
    }
    paste(as.character(value), collapse = "; ")
  }
  data.frame(
    signature_name = signature_name,
    source = get_chr("source"),
    doi = get_chr("doi"),
    state = get_chr("state"),
    species = get_chr("species"),
    context = get_chr("context"),
    interpretation = get_chr("interpretation"),
    stringsAsFactors = FALSE
  )
}

.ScoreSignaturesScoreNames <- function(signature_names, prefix = NULL, suffix = "_score") {
  score_names <- signature_names
  if (!is.null(prefix) && nzchar(prefix)) {
    score_names <- paste0(prefix, "_", score_names)
  }
  if (!is.null(suffix) && nzchar(suffix)) {
    score_names <- paste0(score_names, suffix)
  }
  stats::setNames(make.unique(score_names), signature_names)
}

.ScoreSignaturesScoreMatrix <- function(assay_matrix, signatures, method) {
  score_matrix <- vapply(
    signatures,
    function(signature) {
      up_score <- .ScoreSignaturesScoreGeneSet(assay_matrix, signature$up, method)
      down_score <- .ScoreSignaturesScoreGeneSet(assay_matrix, signature$down, method)
      up_score - down_score
    },
    FUN.VALUE = numeric(ncol(assay_matrix))
  )
  score_matrix <- t(score_matrix)
  colnames(score_matrix) <- colnames(assay_matrix)
  score_matrix
}

.ScoreSignaturesScoreGeneSet <- function(assay_matrix, genes, method) {
  if (length(genes) == 0) {
    return(rep(0, ncol(assay_matrix)))
  }
  gene_matrix <- as.matrix(assay_matrix[genes, , drop = FALSE])

  if (identical(method, "mean")) {
    return(colMeans(gene_matrix))
  }

  if (identical(method, "zscore")) {
    scaled <- t(scale(t(gene_matrix)))
    scaled[is.na(scaled)] <- 0
    return(colMeans(scaled))
  }

  if (identical(method, "pc")) {
    return(gsScore(ExprMat = gene_matrix, GeneList = rownames(gene_matrix), summarizationFunction = "PC"))
  }

  stop("Unsupported matrix scoring method.", call. = FALSE)
}

.ScoreSignaturesScoreSeurat <- function(seu,
                                         signatures,
                                         assay,
                                         slot,
                                         nbin,
                                         ctrl,
                                         seed,
                                         search,
                                         ...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The Seurat package is required for `method = \"seurat\"`.", call. = FALSE)
  }

  up_sets <- lapply(signatures, `[[`, "up")
  down_sets <- lapply(signatures, `[[`, "down")
  up_scores <- .ScoreSignaturesAddModuleScore(
    seu = seu,
    gene_sets = up_sets,
    assay = assay,
    slot = slot,
    nbin = nbin,
    ctrl = ctrl,
    seed = seed,
    search = search,
    ...
  )

  has_down <- lengths(down_sets) > 0
  down_scores <- matrix(
    0,
    nrow = length(signatures),
    ncol = ncol(up_scores),
    dimnames = dimnames(up_scores)
  )
  if (any(has_down)) {
    scored_down <- .ScoreSignaturesAddModuleScore(
      seu = seu,
      gene_sets = down_sets[has_down],
      assay = assay,
      slot = slot,
      nbin = nbin,
      ctrl = ctrl,
      seed = seed,
      search = search,
      ...
    )
    down_scores[rownames(scored_down), ] <- scored_down
  }

  up_scores - down_scores
}

.ScoreSignaturesAddModuleScore <- function(seu, gene_sets, assay, slot, nbin, ctrl, seed, search, ...) {
  temp_name <- "__SignatureScore_"
  while (any(paste0(temp_name, seq_along(gene_sets)) %in% colnames(seu@meta.data))) {
    temp_name <- paste0(temp_name, sample.int(9999, 1), "_")
  }

  scored <- Seurat::AddModuleScore(
    object = seu,
    features = gene_sets,
    assay = assay,
    slot = slot,
    name = temp_name,
    nbin = nbin,
    ctrl = ctrl,
    seed = seed,
    search = search,
    ...
  )

  score_cols <- paste0(temp_name, seq_along(gene_sets))
  missing_cols <- setdiff(score_cols, colnames(scored@meta.data))
  if (length(missing_cols) > 0) {
    stop("Seurat::AddModuleScore did not return the expected score columns.", call. = FALSE)
  }

  score_matrix <- t(as.matrix(scored@meta.data[, score_cols, drop = FALSE]))
  rownames(score_matrix) <- names(gene_sets)
  colnames(score_matrix) <- colnames(seu)
  score_matrix
}

.ScoreSignaturesScoreUCell <- function(assay_matrix, signatures, maxRank, ncores, BPPARAM, ...) {
  if (!requireNamespace("UCell", quietly = TRUE)) {
    stop("The UCell package is required for `method = \"ucell\"`.", call. = FALSE)
  }
  if (is.null(BPPARAM)) {
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      stop("The BiocParallel package is required for default UCell serial execution.", call. = FALSE)
    }
    BPPARAM <- BiocParallel::SerialParam()
  }

  up_sets <- lapply(signatures, `[[`, "up")
  down_sets <- lapply(signatures, `[[`, "down")
  up_scores <- .ScoreSignaturesUCellMatrix(
    assay_matrix = assay_matrix,
    gene_sets = up_sets,
    maxRank = maxRank,
    ncores = ncores,
    BPPARAM = BPPARAM,
    ...
  )

  has_down <- lengths(down_sets) > 0
  down_scores <- matrix(
    0,
    nrow = length(signatures),
    ncol = ncol(up_scores),
    dimnames = dimnames(up_scores)
  )
  if (any(has_down)) {
    scored_down <- .ScoreSignaturesUCellMatrix(
      assay_matrix = assay_matrix,
      gene_sets = down_sets[has_down],
      maxRank = maxRank,
      ncores = ncores,
      BPPARAM = BPPARAM,
      ...
    )
    down_scores[rownames(scored_down), ] <- scored_down
  }

  up_scores - down_scores
}

.ScoreSignaturesUCellMatrix <- function(assay_matrix, gene_sets, maxRank, ncores, BPPARAM, ...) {
  scores <- UCell::ScoreSignatures_UCell(
    matrix = assay_matrix,
    features = gene_sets,
    maxRank = maxRank,
    name = "",
    ncores = ncores,
    BPPARAM = BPPARAM,
    ...
  )
  score_matrix <- t(as.matrix(scores))
  rownames(score_matrix) <- names(gene_sets)
  score_matrix
}

.ScoreSignaturesGetAssayData <- function(seu, assay, layer, slot) {
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
