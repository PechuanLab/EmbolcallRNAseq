#' Convert Mouse Gene Symbols to Human Gene Symbols
#'
#' This function converts a list of mouse gene symbols to their corresponding human gene symbols using the `biomaRt` package. It queries the Ensembl database to perform the mapping.
#'
#' @param x A character vector of mouse gene symbols (e.g., "mgi_symbol").
#'
#' @return A character vector of unique human gene symbols corresponding to the input mouse gene symbols.
#' @export
#'
#' @examples
#' # Example usage:
#' # Assuming `mouse_genes` is a character vector of mouse gene symbols:
#' mouse_genes <- c("Tnf", "Il6", "Cd4")
#' human_genes <- convertMouseGeneList(mouse_genes)
#' print(human_genes)
#'
convertMouseGeneList <- function(x) {
  # Connect to Ensembl databases for human and mouse
  human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  # Map mouse gene symbols to human gene symbols
  genesV2 = biomaRt::getLDS(
    attributes = c("mgi_symbol"), 
    filters = "mgi_symbol", 
    values = x, 
    mart = mouse, 
    attributesL = c("hgnc_symbol"), 
    martL = human, 
    uniqueRows = TRUE
  )

  # Extract unique human gene symbols
  humanx <- unique(genesV2[, 2])

  # Print the first 6 human gene symbols to the console
  print(head(humanx))

  # Return the unique human gene symbols
  return(humanx)
}