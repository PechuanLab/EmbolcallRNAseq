#' Title
#'
#' @param ExprMat
#' @param Prefix
#' @param Contrastes
#' @param annotations
#' @param Treatment
#' @param species
#' @param logfold
#'
#' @return Master wrapper for multiple functions
#' @export
#'
#' @examples BulkPlots()

BulkPlots <- function(ExprMat,Prefix,Contrastes,annotations,Treatment,species,logfold){

  # Put everything in a folder and be tidy
  dir = getwd()
  onestan = paste(getwd(),paste(Prefix,"DEPlots",sep="_"),sep="/")
  dir.create(onestan)
  setwd(onestan)

  # Loop for all contrasts
  for (i in 1:length(Contrastes)) {
    # Get the contrast create the directory
    contraste = Contrastes[i]
    onestan2 = paste(onestan,contraste,sep="/")
    dir.create(onestan2)
    setwd(onestan2)
    
    # prepare the fitted limma object
    fit3 = limma::treat(fit2,lfc = 0)
    signature = limma::topTreat(fit3, coef=contraste, n=Inf,adjust.method="fdr",sort.by = "P") 
    # save the whole thing in case you would like to run something else
    write.csv(signature, paste0(Prefix,contraste,"allGenes.csv"))
    #Bulk Plot
    try(BulkPlot(signature,Prefix,contraste,species),silent =T)

  }
}
