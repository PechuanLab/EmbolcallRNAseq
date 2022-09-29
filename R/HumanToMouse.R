


#' Title
#'
#' @param x
#'
#' @return Gene vector
#' @export
#'
#' @examples HumanToMouse( )

HumanToMouse <- function(x){

  human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  genesV2 = biomaRt::getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol",
                   values = x , mart = human, attributesL = c("mgi_symbol"),
                   martL = mouse, uniqueRows=T)
  mousex <- unique(genesV2[, 2])

  # Print the first 6 genes found to the screen
  print(head(humanx))
  return
}
