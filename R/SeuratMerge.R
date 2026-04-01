#' Merge Multiple Seurat Objects
#'
#' Read multiple Seurat .rds files, preserve RNA assay and optionally include
#' additional assays, merge objects sequentially and join layers.
#'
#' @param seurats Character vector of file paths to Seurat .rds files.
#' @param additional_assays Character or NULL. Name of an additional assay to
#'   preserve (e.g. "CiteSeq"). If NULL only RNA is kept. Default NULL.
#' @return A merged Seurat object with joined layers.
#' @export
#'
#' @examples
#' \dontrun{
#' SeuratMerge(c("s1.rds", "s2.rds"), additional_assays = "CiteSeq")
#' }
SeuratMerge <- function(seurats, additional_assays = NULL) {
  
  # helper: safely retrieve assay names as a plain character vector
  get_assay_names <- function(obj) as.character(Seurat::Assays(obj))
  
  # Read the first Seurat object to seed the merged object
  seu_i <- readRDS(seurats[1])
  seu_n <- CreateSeuratObject(counts = seu_i[["RNA"]]$counts,
                              meta.data = seu_i@meta.data)
  
  # Only add the additional assay to the first object if requested and present
  if (!is.null(additional_assays) &&
      is.character(additional_assays) &&
      additional_assays %in% get_assay_names(seu_i)) {
    seu_n[[additional_assays]] <- CreateAssayObject(counts = seu_i[[additional_assays]]$counts)
  }
  
  # free the raw first object
  rm(seu_i)
  
  # Loop through the remaining files and merge sequentially
  for (i in 2:length(seurats)) {
    seu_i_j <- readRDS(seurats[i])
    seu_i_j_n <- CreateSeuratObject(counts = seu_i_j[["RNA"]]$counts,
                                    meta.data = seu_i_j@meta.data)
    
    # Only add the additional assay if requested and present in this object
    if (!is.null(additional_assays) &&
        is.character(additional_assays) &&
        additional_assays %in% get_assay_names(seu_i_j)) {
      seu_i_j_n[[additional_assays]] <- CreateAssayObject(counts = seu_i_j[[additional_assays]]$counts)
    }
    
    # free raw object before merge to reduce memory pressure
    rm(seu_i_j)
    
    # Merge into accumulating object and join layers
    seu_n <- merge(seu_n, seu_i_j_n) %>% JoinLayers
    rm(seu_i_j_n)
  }
  
  # Final join layers on the fully merged object
  seu_n <- JoinLayers(seu_n)
  
  return(seu_n)
}