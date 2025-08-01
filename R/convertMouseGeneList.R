#' Title
#'
#' @param x
#'
#' @return Unique gene
#' @export
#'
#' @examples convertMouseGeneList()
convertMouseGeneList <- function(x){

  human =biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  genesV2 = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])

  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
