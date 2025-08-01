#' Create a Standard Project Directory Structure
#'
#' This function creates a predefined directory structure for a computational biology project
#' under a user-specified top-level directory. The structure includes folders for data,
#' source code, QC, clustering, analysis, figures, slides, and paper drafts.
#'
#' @param top_dir Character. The name of the top-level directory where the structure will be created.
#' @return Creates directories on disk. Returns a message upon successful creation.
#' @examples
#' create_project_structure("MyProject")
create_project_structure <- function(top_dir = "NewProject") {
  # Define the relative paths of all subdirectories
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
  
  # Create each directory recursively
  for (d in dirs) {
    dir.create(file.path(top_dir, d), recursive = TRUE, showWarnings = FALSE)
  }
  
  message("✅ Directory structure successfully created under: ", top_dir)
}