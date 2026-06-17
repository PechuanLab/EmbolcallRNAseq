#' Compute MSigDB pathway activities for a Seurat object
#'
#' Scores MSigDB gene sets in a Seurat object and stores the resulting
#' pathway-by-cell matrix as a new assay. By default, Hallmark pathways are
#' loaded from MSigDB for the inferred organism.
#'
#' @param seu Seurat object.
#' @param species Character. One of \code{"auto"}, \code{"human"},
#'   \code{"mouse"}, \code{"Homo sapiens"}, or \code{"Mus musculus"}. When
#'   \code{"auto"}, the species is inferred from overlap between assay features
#'   and human/mouse MSigDB genes. Default \code{"auto"}.
#' @param collection Character vector. MSigDB collection(s) to use. Default
#'   \code{"H"} for Hallmark.
#' @param subcollection Optional character vector of MSigDB subcollection(s) to
#'   use.
#' @param gene_sets Optional named list of gene symbols. When supplied, these
#'   gene sets are used instead of fetching MSigDB.
#' @param method Character. Scoring backend: \code{"seurat"} for
#'   \code{Seurat::AddModuleScore} or \code{"ucell"} for UCell rank-based
#'   scores. Aliases \code{"module"}, \code{"addmodulescore"}, \code{"u"}, and
#'   \code{"u_module"} are accepted.
#' @param assay Character. Assay to score. Defaults to the active assay.
#' @param layer Optional Seurat v5 layer to use when extracting data for
#'   filtering and UCell scoring. For \code{method = "seurat"}, this value is
#'   passed to \code{Seurat::AddModuleScore} as its \code{slot} argument.
#' @param slot Character. Seurat v4 slot/layer fallback. Defaults to
#'   \code{"data"} for Seurat module scores and \code{"counts"} for UCell.
#' @param assay_name Character. Name of the assay to add. Default
#'   \code{"Pathway_Activities"}.
#' @param min_genes Integer. Minimum number of genes from a pathway that must be
#'   present in the selected assay. Default \code{5}.
#' @param overwrite Logical. If \code{TRUE}, replace an existing assay named
#'   \code{assay_name}. Default \code{TRUE}.
#' @param nbin Integer. Number of expression bins passed to
#'   \code{Seurat::AddModuleScore}. Default \code{24}.
#' @param ctrl Integer. Number of control genes passed to
#'   \code{Seurat::AddModuleScore}. Default \code{100}.
#' @param seed Integer. Random seed passed to \code{Seurat::AddModuleScore}.
#'   Default \code{1}.
#' @param search Logical. Passed to \code{Seurat::AddModuleScore}. Default
#'   \code{FALSE}.
#' @param maxRank Integer. Maximum rank passed to
#'   \code{UCell::ScoreSignatures_UCell}. Default \code{1500}.
#' @param ncores Integer. Number of cores passed to UCell. Default \code{1}.
#' @param BPPARAM Optional BiocParallel parameter object for UCell. If
#'   \code{NULL}, \code{BiocParallel::SerialParam()} is used.
#' @param verbose Logical. If \code{TRUE}, report selected species and retained
#'   pathway counts. Default \code{TRUE}.
#' @param ... Additional arguments passed to the selected scoring backend.
#'
#' @return The input Seurat object with a new assay containing pathway
#'   activities in the \code{data} layer/slot. Original MSigDB pathway names and
#'   retained gene counts are stored in the pathway assay feature metadata.
#' @export
#'
#' @examples
#' \dontrun{
#' seu <- MSigDBPathwayActivities(seu)
#' seu <- MSigDBPathwayActivities(seu, species = "mouse", method = "ucell")
#' seu <- MSigDBPathwayActivities(seu, collection = "C2", subcollection = "CP:REACTOME")
#' }
MSigDBPathwayActivities <- function(seu,
                                    species = "auto",
                                    collection = "H",
                                    subcollection = NULL,
                                    gene_sets = NULL,
                                    method = c("seurat", "ucell"),
                                    assay = NULL,
                                    layer = NULL,
                                    slot = NULL,
                                    assay_name = "Pathway_Activities",
                                    min_genes = 5,
                                    overwrite = TRUE,
                                    nbin = 24,
                                    ctrl = 100,
                                    seed = 1,
                                    search = FALSE,
                                    maxRank = 1500,
                                    ncores = 1,
                                    BPPARAM = NULL,
                                    verbose = TRUE,
                                    ...) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("The SeuratObject package is required.", call. = FALSE)
  }
  if (!inherits(seu, "Seurat")) {
    stop("`seu` must be a Seurat object.", call. = FALSE)
  }
  method <- .MSigDBPathwayActivitiesNormalizeMethod(method)
  if (is.null(slot)) {
    slot <- if (method == "ucell") "counts" else "data"
  }
  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(seu)
  }
  if (assay_name %in% SeuratObject::Assays(seu) && !isTRUE(overwrite)) {
    stop("Assay `", assay_name, "` already exists. Use `overwrite = TRUE` to replace it.", call. = FALSE)
  }
  if (!is.numeric(min_genes) || length(min_genes) != 1 || min_genes < 1) {
    stop("`min_genes` must be a positive integer.", call. = FALSE)
  }

  assay_matrix <- .MSigDBPathwayActivitiesGetAssayData(
    seu = seu,
    assay = assay,
    layer = layer,
    slot = slot
  )
  available_features <- rownames(assay_matrix)
  if (length(available_features) == 0) {
    stop("The selected assay contains no features.", call. = FALSE)
  }

  if (is.null(gene_sets)) {
    species <- .MSigDBPathwayActivitiesNormalizeSpecies(species)
    if (identical(species, "auto")) {
      species <- .MSigDBPathwayActivitiesInferSpecies(
        features = available_features,
        collection = collection,
        subcollection = subcollection
      )
    }
    gene_sets <- .MSigDBPathwayActivitiesFetchGeneSets(
      species = species,
      collection = collection,
      subcollection = subcollection
    )
  } else {
    gene_sets <- .MSigDBPathwayActivitiesValidateGeneSets(gene_sets)
    species <- .MSigDBPathwayActivitiesNormalizeSpecies(species, allow_auto = TRUE)
  }

  filtered <- .MSigDBPathwayActivitiesFilterGeneSets(
    gene_sets = gene_sets,
    available_features = available_features,
    min_genes = as.integer(min_genes)
  )
  gene_sets <- filtered$gene_sets
  pathway_metadata <- filtered$metadata

  if (isTRUE(verbose)) {
    message(
      "Scoring ", length(gene_sets), " pathway(s) with method `", method, "`",
      if (!identical(species, "auto")) paste0(" for ", species) else "",
      "."
    )
  }

  default_assay <- SeuratObject::DefaultAssay(seu)
  score_matrix <- switch(
    method,
    seurat = .MSigDBPathwayActivitiesScoreSeurat(
      seu = seu,
      gene_sets = gene_sets,
      assay = assay,
      slot = if (!is.null(layer)) layer else slot,
      nbin = nbin,
      ctrl = ctrl,
      seed = seed,
      search = search,
      ...
    ),
    ucell = .MSigDBPathwayActivitiesScoreUCell(
      assay_matrix = assay_matrix,
      gene_sets = gene_sets,
      maxRank = maxRank,
      ncores = ncores,
      BPPARAM = BPPARAM,
      ...
    )
  )
  SeuratObject::DefaultAssay(seu) <- default_assay

  score_matrix <- .MSigDBPathwayActivitiesPrepareScoreMatrix(score_matrix)
  pathway_metadata <- pathway_metadata[rownames(score_matrix), , drop = FALSE]
  rownames(score_matrix) <- .MSigDBPathwayActivitiesSeuratFeatureNames(rownames(score_matrix))
  rownames(pathway_metadata) <- rownames(score_matrix)
  pathway_metadata$method <- method
  pathway_metadata$species <- species
  pathway_metadata$collection <- paste(collection, collapse = ",")
  pathway_metadata$subcollection <- if (is.null(subcollection)) NA_character_ else paste(subcollection, collapse = ",")

  pathway_assay <- SeuratObject::CreateAssayObject(data = score_matrix)
  pathway_assay <- SeuratObject::AddMetaData(pathway_assay, metadata = pathway_metadata)
  seu[[assay_name]] <- pathway_assay
  SeuratObject::DefaultAssay(seu) <- default_assay
  seu
}

.MSigDBPathwayActivitiesNormalizeMethod <- function(method) {
  method <- tolower(as.character(method[1]))
  if (method %in% c("seurat", "module", "addmodulescore", "add_module_score")) {
    return("seurat")
  }
  if (method %in% c("ucell", "u", "u_module", "umodule", "u-module")) {
    return("ucell")
  }
  stop("`method` must be one of `seurat` or `ucell`.", call. = FALSE)
}

.MSigDBPathwayActivitiesNormalizeSpecies <- function(species, allow_auto = TRUE) {
  species <- tolower(gsub("_", " ", as.character(species[1])))
  if (allow_auto && species %in% c("auto", "automatic", "infer", "inferred")) {
    return("auto")
  }
  if (species %in% c("human", "homo sapiens", "hs", "hsa")) {
    return("Homo sapiens")
  }
  if (species %in% c("mouse", "mus musculus", "mm", "mmu")) {
    return("Mus musculus")
  }
  stop("`species` must be `auto`, `human`, `mouse`, `Homo sapiens`, or `Mus musculus`.", call. = FALSE)
}

.MSigDBPathwayActivitiesInferSpecies <- function(features, collection, subcollection) {
  human_sets <- tryCatch(
    .MSigDBPathwayActivitiesFetchGeneSets("Homo sapiens", collection, subcollection),
    error = function(e) NULL
  )
  mouse_sets <- tryCatch(
    .MSigDBPathwayActivitiesFetchGeneSets("Mus musculus", collection, subcollection),
    error = function(e) NULL
  )
  human_overlap <- if (is.null(human_sets)) 0 else sum(features %in% unique(unlist(human_sets, use.names = FALSE)))
  mouse_overlap <- if (is.null(mouse_sets)) 0 else sum(features %in% unique(unlist(mouse_sets, use.names = FALSE)))
  if (human_overlap > mouse_overlap) {
    return("Homo sapiens")
  }
  if (mouse_overlap > human_overlap) {
    return("Mus musculus")
  }

  feature_sample <- utils::head(features, 5000)
  human_case <- sum(grepl("[A-Z]", feature_sample) & !grepl("[a-z]", feature_sample))
  mouse_case <- sum(grepl("^[A-Z][a-z]", feature_sample))
  if (human_case > mouse_case) {
    return("Homo sapiens")
  }
  if (mouse_case > human_case) {
    return("Mus musculus")
  }

  warning("Could not confidently infer species from assay features; using human MSigDB.", call. = FALSE)
  "Homo sapiens"
}

.MSigDBPathwayActivitiesFetchGeneSets <- function(species, collection, subcollection) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("The msigdbr package is required to fetch MSigDB gene sets.", call. = FALSE)
  }
  msigdbr_args <- list(species = species)
  msigdbr_formals <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% msigdbr_formals) {
    msigdbr_args$collection <- collection
  } else if ("category" %in% msigdbr_formals) {
    msigdbr_args$category <- collection
  }
  if (!is.null(subcollection)) {
    if ("subcollection" %in% msigdbr_formals) {
      msigdbr_args$subcollection <- subcollection
    } else if ("subcategory" %in% msigdbr_formals) {
      msigdbr_args$subcategory <- subcollection
    }
  }

  msigdb <- do.call(msigdbr::msigdbr, msigdbr_args)
  collection_col <- .MSigDBPathwayActivitiesFirstColumn(msigdb, c("gs_collection", "gs_cat", "collection"))
  subcollection_col <- .MSigDBPathwayActivitiesFirstColumn(msigdb, c("gs_subcollection", "gs_subcat", "subcollection"))
  pathway_col <- .MSigDBPathwayActivitiesFirstColumn(msigdb, c("gs_name", "name", "pathway"))
  gene_col <- .MSigDBPathwayActivitiesFirstColumn(msigdb, c("gene_symbol", "db_gene_symbol", "human_gene_symbol"))

  if (!is.null(collection_col)) {
    msigdb <- msigdb[msigdb[[collection_col]] %in% collection, , drop = FALSE]
  }
  if (!is.null(subcollection) && !is.null(subcollection_col)) {
    msigdb <- msigdb[msigdb[[subcollection_col]] %in% subcollection, , drop = FALSE]
  }
  if (is.null(pathway_col) || is.null(gene_col)) {
    stop("Could not identify pathway and gene columns in msigdbr output.", call. = FALSE)
  }
  if (nrow(msigdb) == 0) {
    stop("No MSigDB gene sets were found for the requested collection/subcollection.", call. = FALSE)
  }

  gene_sets <- split(as.character(msigdb[[gene_col]]), as.character(msigdb[[pathway_col]]))
  .MSigDBPathwayActivitiesValidateGeneSets(gene_sets)
}

.MSigDBPathwayActivitiesValidateGeneSets <- function(gene_sets) {
  if (!is.list(gene_sets) || is.null(names(gene_sets)) || any(!nzchar(names(gene_sets)))) {
    stop("`gene_sets` must be a named list of character vectors.", call. = FALSE)
  }
  gene_sets <- lapply(gene_sets, function(genes) {
    genes <- unique(stats::na.omit(as.character(genes)))
    genes[nzchar(genes)]
  })
  gene_sets <- gene_sets[lengths(gene_sets) > 0]
  if (length(gene_sets) == 0) {
    stop("No non-empty gene sets were supplied.", call. = FALSE)
  }
  gene_sets
}

.MSigDBPathwayActivitiesFilterGeneSets <- function(gene_sets, available_features, min_genes) {
  original_counts <- lengths(gene_sets)
  gene_sets <- lapply(gene_sets, intersect, y = available_features)
  retained_counts <- lengths(gene_sets)
  keep <- retained_counts >= min_genes
  if (!any(keep)) {
    stop(
      "No pathways retained at least ", min_genes,
      " gene(s) present in the selected assay.",
      call. = FALSE
    )
  }
  dropped <- sum(!keep)
  if (dropped > 0) {
    warning(dropped, " pathway(s) dropped because too few genes were present in the assay.", call. = FALSE)
  }
  gene_sets <- gene_sets[keep]
  metadata <- data.frame(
    pathway_name = names(gene_sets),
    original_n_genes = as.integer(original_counts[keep]),
    retained_n_genes = as.integer(retained_counts[keep]),
    stringsAsFactors = FALSE,
    row.names = names(gene_sets)
  )
  list(gene_sets = gene_sets, metadata = metadata)
}

.MSigDBPathwayActivitiesScoreSeurat <- function(seu, gene_sets, assay, slot, nbin, ctrl, seed, search, ...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The Seurat package is required for `method = \"seurat\"`.", call. = FALSE)
  }
  temp_name <- "__PathwayActivities_"
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

.MSigDBPathwayActivitiesScoreUCell <- function(assay_matrix, gene_sets, maxRank, ncores, BPPARAM, ...) {
  if (!requireNamespace("UCell", quietly = TRUE)) {
    stop("The UCell package is required for `method = \"ucell\"`.", call. = FALSE)
  }
  if (is.null(BPPARAM)) {
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      stop("The BiocParallel package is required for default UCell serial execution.", call. = FALSE)
    }
    BPPARAM <- BiocParallel::SerialParam()
  }
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

.MSigDBPathwayActivitiesGetAssayData <- function(seu, assay, layer, slot) {
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

.MSigDBPathwayActivitiesPrepareScoreMatrix <- function(score_matrix) {
  score_matrix <- as.matrix(score_matrix)
  storage.mode(score_matrix) <- "double"
  if (any(!is.finite(score_matrix))) {
    stop("Pathway score matrix contains non-finite values.", call. = FALSE)
  }
  score_matrix
}

.MSigDBPathwayActivitiesSeuratFeatureNames <- function(pathway_names) {
  make.unique(gsub("_", "-", pathway_names, fixed = TRUE))
}

.MSigDBPathwayActivitiesFirstColumn <- function(data, columns) {
  found <- intersect(columns, colnames(data))
  if (length(found) == 0) {
    return(NULL)
  }
  found[1]
}
