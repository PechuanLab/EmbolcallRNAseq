#' Cell-type signature lookup (mouse or human)
#'
#' Return cell-type signature lists or a specific subsignature. Choose species
#' to get mouse-style (Title-case) or human-style (UPPERCASE) gene symbols.
#'
#' @param signature_name Character. Group name (e.g. "TCells") or subsignature
#'   (e.g. "Effector"). Case-insensitive.
#' @param species Character. "mouse" or "human". Determines symbol casing.
#'   Default: "mouse".
#' @return If a group name is provided, a named list of subsignatures; if a
#'   subsignature name is provided, a character vector of gene symbols. NULL if
#'   not found.
#' @export
scSignatures <- function(signature_name, species = c("mouse", "human")) {
  species <- match.arg(species)

  to_title <- function(vec) {
    vapply(vec, function(g) {
      g2 <- tolower(g)
      sub("^(.)", "\\U\\1", g2, perl = TRUE)
    }, FUN.VALUE = character(1))
  }

  # base (mouse-style) groups
  TCells <- list(
    Exhaustion = c("Pdcd1","Tigit","Tox","Lag3","Entpd1","Tnfrsf9","Havcr2","Ctla4"),
    Effector = c("Gzma","Gzmb","Gzmk","Gzmm","Prf1","Nkg7","Ifng","Tnf","Klrd1"),
    Coestimulatory = c("Tnfrsf4","Icos","Cd28"),
    NaiveMemory = c("Il7r","Tcf7","Sell","Ccr7","Lef1"),
    Th1 = c("Tbx21","Il12rb1","Il12rb2","Stat4"),
    Th2 = c("Gata3","Stat6","Il4ra")
  )

  Fibroblasts <- list(
    PanFibroblasts = c("Pdgfra","Postn","Dcn","Col1a1","Lum","Fbln1"),
    Mesothelial = c("Lgals7","Upk3b","Lrrn4","Ptger3","Wt1","Krt19"),
    Adventitial = c("Pi16","Dpt","Cd34","Clec3b","Ly6c1"),
    Myofibroblasts = c("Acta2","Tagln","Myh11","Cnn1","Lrrc15"),
    Inflammatory = c("Cxcl5","Rarres2","Ccl2","Ccl7","Cxcl12"),
    Retinol = c("Rbp1","Rbp4","Stra6","Lrat","Rdh10","Aldh1a1","Aldh1a2","Crabp1","Crabp2","Rara","Rxra","Cyp26b1","Rarres2"),
    Others = c("Grem1","Spp1")
  )

  Macrophages <- list(
    M1_Inflammatory = c("Il1b","Tnf","Nos2","Il6","Il12b"),
    M2_Alternative = c("Mrc1","Cd163","Arg1","Il10","Chil3"),
    Resident = c("Adgre1","C1qa","C1qb","C1qc","Csf1r")
  )

  NSCLC_TumorCells <- list(
    PulmonaryTFs = c("Nkx2-1","Sox9","Sox2","Foxa2","Cebpa","Cebpb","Foxa1"),
    AT1 = c("Hopx","Pdpn","Nkx2-1","Aqp5","Ager"),
    AT2 = c("Sftpc","Sftpb","Sftpa1","Lamp3","Napsa"),
    BasalCells = c("Krt5","Krt6a","Trp63","Upk1b"),
    CiliatedCells = c("Foxj1","Tppp3","Tuba1a","Cfap43","Ccdc78","Ccdc113"),
    PulmonaryNeuroendocrine = c("Ascl1","Calca","Grp","Cdh18","Nrxn1"),
    GobletCells = c("Tff2","Muc5ac","Muc5b","Spdef","Agr2"),
    ClubCells = c("Scgb1a1","Scgb3a2","Cyp2f2"),
    Neuroendocrine = c("Chga","Chgb"),
    Mesenchymal = c("Vim","Cdh2","Fn1","Snai1","Snai2"),
    Epithelial = c("Cdh1","Epcam","Krt19","Krt8")
  )

  GastricCells <- list(
    Parietal = c("Atp4a","Atp4b","Clic6","Ckb"),
    Chief = c("Pga5","Pgc","Cblif"),
    PitCells = c("Dmbt1","Lypd8","Agr2","Ctse","Dpcr1"),
    SurfaceMucous = c("Muc5ac","Tff1","Gkn1"),
    MucousNeck = c("Muc6","Tff2","Gkn2"),
    Enteroendocrine = c("Chga","Chgb","Tph1","Neurod1"),
    StemProgenitor = c("Lgr5","H19"),
    Tuft = c("Pou2f3","Dclk1","Trpm5"),
    SPEM = c("Cd44","Tff2","Gkn3")
  )

  Esophagus <- list(
    Basal = c("Krt5","Krt14","Trp63","Krt15"),
    Suprabasal = c("Krt4","Krt13","Ivl","Flg"),
    Differentiated = c("Krt1","Krt10","Lor"),
    Goblet_Barrett = c("Muc2","Tff3","Cdx2"),
    SmoothMuscle = c("Acta2","Myh11","Myl9"),
    Endothelial = c("Pecam1","Vwf","Eng")
  )

  Colon <- list(
    Enterocyte = c("Alpi","Slc26a3","Fabp2","Tmem37"),
    Goblet = c("Muc2","Fcgbp","Tff3","Agr2","Spdef","Zg16","Clca1","Spink4"),
    Paneth_or_PanethLike = c("Lyz","Ca4","Ca7","Spib"),
    Enteroendocrine = c("Chga","Chgb","Cpe","Neurod1","Pyy"),
    Progenitor_Regenerative = c("Sox9","Cdk6","Muc4","Fabp5","Pla2g2a","Lcn2"),
    TA_Cycling = c("Mki67","Top2a","Pcna","Ccna2","Mcm5"),
    TA_Cycling_Alt = c("Nusap1","Ube2c","Cenpf","Aurkb","Cdk1"),
    Stem_WNT = c("Lgr5","Ascl2","Smoc2","Axin2","Apcdd1","Nkd1"),
    Tuft = c("Pou2f3","Trpm5","Gfi1b","Dclk1")
  )

  PancreaticCancer <- list(
    Tumor_Epithelial = c("Krt19","Epcam","Sox9","Pdx1","Kras"),
    Cancer_Associated_Fibroblasts = c("Acta2","Pdgfra","Postn","Fap","Col1a1"),
    Immune = c("Cd68","Adgre1","Cd3e","Cd4","Cd8a"),
    Acinar = c("Cpa1","Ptf1a","Cela3b","Amy2a3"),
    Endocrine = c("Ins1","Gcg","Sst","Ppy")
  )

  base_groups <- list(
    TCells = TCells,
    Fibroblasts = Fibroblasts,
    Macrophages = Macrophages,
    NSCLC_TumorCells = NSCLC_TumorCells,
    GastricCells = GastricCells,
    Esophagus = Esophagus,
    Colon = Colon,
    PancreaticCancer = PancreaticCancer
  )

  # convert according to species
  if (species == "human") {
    groups <- lapply(base_groups, function(sub) lapply(sub, toupper))
  } else {
    groups <- lapply(base_groups, function(sub) lapply(sub, to_title))
  }

  if (missing(signature_name) || length(signature_name) != 1) {
    stop("Please provide a single signature_name (group or subsignature).")
  }

  sig_lower <- tolower(signature_name)

  # match group
  grp_match <- names(groups)[tolower(names(groups)) == sig_lower]
  if (length(grp_match) == 1) return(groups[[grp_match]])

  # match subsignature
  for (g in names(groups)) {
    subs <- names(groups[[g]])
    idx <- which(tolower(subs) == sig_lower)
    if (length(idx) == 1) return(groups[[g]][[idx]])
  }

  return(NULL)
}