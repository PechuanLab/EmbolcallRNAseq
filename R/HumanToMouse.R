


##' Convert human gene symbols to mouse gene symbols
##'
##' Maps a vector of human HGNC symbols to mouse MGI symbols using biomaRt and the
##' Ensembl datasets for Homo sapiens and Mus musculus.
##'
##' @param x character() Vector of human HGNC gene symbols.
##'
##' @return character() Unique vector of mapped mouse MGI symbols.
##' @export
##'
##' @examples
##' \dontrun{
##'   HumanToMouse(c("TNF", "IL6", "CD4"))
##' }
HumanToMouse <- function(x){

  human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  genesV2 = biomaRt::getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol",
                   values = x , mart = human, attributesL = c("mgi_symbol"),
                   martL = mouse, uniqueRows=TRUE)
  mousex <- unique(genesV2[, 2])

  # Print the first 6 mapped mouse genes
  print(utils::head(mousex))
  return(mousex)
}
