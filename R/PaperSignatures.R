#' Curated paper-derived gene signatures
#'
#' Returns curated signatures collected from papers. Signatures are represented
#' as named lists with \code{up}, optional \code{down}, and \code{metadata}
#' fields so they can be scored by \code{\link{ScoreSignatures}}.
#'
#' The first registry entry captures the central WNT/MAPK cell-state model from
#' Moore et al. 2026. These are intentionally compact marker signatures derived
#' from the paper text and figure labels; they can be expanded later with the
#' full supplementary gene modules when those tables are available in the
#' package.
#'
#' @param paper Character. Paper key to return. Currently accepts
#'   \code{"moore2026_crc_states"}, \code{"moore2026"}, or \code{"all"}.
#' @param species Character. \code{"mouse"} returns mouse-style gene symbols;
#'   \code{"human"} returns uppercase human-style symbols.
#' @param include_metadata Logical. If \code{TRUE}, keep per-signature metadata.
#'
#' @return Named list of signatures. Each signature is a list with \code{up},
#'   \code{down}, and optionally \code{metadata}.
#' @export
#'
#' @examples
#' sigs <- PaperSignatures("moore2026")
#' names(sigs)
#' sigs$Moore2026_WNT_StemLike$up
PaperSignatures <- function(paper = c("moore2026_crc_states", "moore2026", "all"),
                            species = c("mouse", "human"),
                            include_metadata = TRUE) {
  paper <- match.arg(paper)
  species <- match.arg(species)

  signatures <- .PaperSignaturesRegistry()
  if (!identical(paper, "all")) {
    keep <- vapply(
      signatures,
      function(signature) paper %in% signature$metadata$aliases,
      FUN.VALUE = logical(1)
    )
    signatures <- signatures[keep]
  }

  signatures <- lapply(signatures, .PaperSignaturesConvertSpecies, species = species)

  if (!isTRUE(include_metadata)) {
    signatures <- lapply(signatures, function(signature) {
      signature$metadata <- NULL
      signature
    })
  }

  signatures
}

.PaperSignaturesRegistry <- function() {
  source_moore2026 <- list(
    paper_key = "moore2026_crc_states",
    aliases = c("moore2026_crc_states", "moore2026"),
    source = "Moore et al. Nature Genetics 2026",
    title = paste(
      "Dynamic transitioning between MAPK-driven and WNT-driven cell states",
      "drives intestinal cancer and shapes therapy response"
    ),
    doi = "10.1038/s41588-026-02611-0",
    species = "mouse",
    context = "intestinal colorectal cancer",
    curation = paste(
      "Compact marker signature curated from the article text and figure labels;",
      "expand with supplementary gene modules when available."
    )
  )

  list(
    Moore2026_WNT_StemLike = .PaperSignature(
      up = c("Lgr5", "Axin2", "Nkd1", "Smoc2", "Sox9", "Ascl2", "Apcdd1"),
      down = c("Dusp6", "Anxa10", "Fosl1", "Etv4", "Etv5", "Spry2", "Spry4"),
      metadata = c(
        source_moore2026,
        list(
          state = "WNT-driven stem-like",
          interpretation = paste(
            "High scores mark a WNT/stem-like state associated with Lgr5+",
            "identity, tumor initiation, and persistence after MAPK inhibition."
          ),
          evidence = paste(
            "Paper highlights WNT target/stem genes including Lgr5, Axin2,",
            "Nkd1 and Smoc2, with inverse behavior relative to MAPK/YAP genes."
          )
        )
      )
    ),
    Moore2026_MAPK_TA_Like = .PaperSignature(
      up = c(
        "Dusp6", "Anxa10", "Fosl1", "Fos", "Jun", "Etv4", "Etv5",
        "Spry2", "Spry4", "Ccnd3"
      ),
      down = c("Lgr5", "Axin2", "Nkd1", "Smoc2", "Sox9", "Ascl2", "Apcdd1"),
      metadata = c(
        source_moore2026,
        list(
          state = "MAPK-driven transit-amplifying-like",
          interpretation = paste(
            "High scores mark a MAPK/YAP, AP-1-associated proliferative",
            "state linked to tumor growth and regeneration."
          ),
          evidence = paste(
            "Paper marks the MAPK state with Dusp6, Anxa10/module 13,",
            "AP-1/FOSL1 activity and proliferation genes including Ccnd3."
          )
        )
      )
    ),
    Moore2026_TA_Proliferation = .PaperSignature(
      up = c("Ccnd3", "Mki67", "Top2a", "Pcna", "Ccna2", "Mcm5"),
      down = character(),
      metadata = c(
        source_moore2026,
        list(
          state = "transit-amplifying/proliferation",
          interpretation = paste(
            "High scores mark cycling or transit-amplifying-like tumor cells;",
            "use together with the MAPK_TA_Like score to separate pathway",
            "activation from generic proliferation."
          ),
          evidence = paste(
            "Paper identifies proliferation and transit-amplifying-like tumor",
            "programs, including Ccnd3, in MAPK-shifted cells."
          )
        )
      )
    )
  )
}

.PaperSignature <- function(up, down = character(), metadata = list()) {
  list(
    up = unique(as.character(up)),
    down = unique(as.character(down)),
    metadata = metadata
  )
}

.PaperSignaturesConvertSpecies <- function(signature, species) {
  convert <- switch(
    species,
    mouse = function(x) x,
    human = toupper
  )
  signature$up <- unique(convert(signature$up))
  signature$down <- unique(convert(signature$down))
  if (!is.null(signature$metadata)) {
    signature$metadata$species <- species
  }
  signature
}
