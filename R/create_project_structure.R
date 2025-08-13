#' Create a Standard Project Directory Structure
#'
#' This function creates a predefined directory structure for a computational biology project
#' under a user-specified top-level directory. The structure varies depending on the project type:
#' "scRNAseq" or "BulkRNAseq".
#'
#' @param top_dir Character. The name of the top-level directory where the structure will be created.
#' @param project_type Character. Type of project: "scRNAseq" or "BulkRNAseq".
#' @return Creates directories on disk. Returns a message upon successful creation.
#' @export
#' @examples
#' create_project_structure(top_dir = "MyProject", project_type = "scRNAseq")

create_project_structure <- function(top_dir = "NewProject", project_type = "scRNAseq") {
  if (!(project_type %in% c("scRNAseq", "BulkRNAseq"))) {
    stop("❌ Invalid project_type. Choose either 'scRNAseq' or 'BulkRNAseq'.")
  }

  if (project_type == "scRNAseq") {
    dirs <- c(
      "01_data/01_rawdata",
      "01_data/02_cellbender_outputs",
      "01_data/03_individualseurats",
      "01_data/04_clusteringseurats",
      "01_data/05_clusteredseurats",
      "01_data/06_finalatlas",
      "01_data/07_additionalobjects",
      "02_src/01_cellranger",
      "02_src/02_cellbender",
      "02_src/03_qc",
      "02_src/04_clustering",
      "02_src/05_analysis",
      "03_qc",
      "04_clustering",
      "05_analysis",
      "06_figures",
      "07_slides",
      "08_paper"
    )
  } else if (project_type == "BulkRNAseq") {
    dirs <- c(
      "01_data/01_raw_fastq",
      "01_data/02_trimmed_fastq",
      "01_data/03_counts",
      "01_data/04_deglist",
      "02_src/01_alignment",
      "02_src/02_quantification",
      "02_src/03_differential_expression",
      "03_qc",
      "04_analysis/02_EDA",
      "04_analysis/02_DE_analysis",
      "05_figures",
      "06_slides",
      "07_paper"
    )
  }

  for (d in dirs) {
    dir.create(file.path(top_dir, d), recursive = TRUE, showWarnings = FALSE)
  }

  message("✅ ", project_type, " directory structure successfully created under: ", top_dir)
}

