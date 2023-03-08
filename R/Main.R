#' Main function.
#'
#' This function constructs the harmonized map
#' @param inFile The path for the input dbf file
#' @return  The output is a list with length 4:
#' (1) v_aux: matrice di k+1 elementi (colonna 1:nome v. ausiliarie, colonna 2: 1 se la v. ausiliaria ? selezionata e 0 altrimenti)
#' (2) coeff: numero delle v. ausiliarie selezionalte
#' (3) founded_a: valore di alpha scelto dalla LOOCV
#' (4) estim: coordinate delle unit? areali e relativo valore armonizzato della variabile di interesse. Da estim si può plottare la mappa armonizzata
#' @keywords Akaike, variables selection, harmonization
#' @export
#' @examples
#' h <- Main()

Main <- function(
    inFile = system.file("data", "data.dbf", package = "dataDriven"),
    incProb = "pj",
    depVar = "fj",
    varsToRemove =  c("ID", "IDCLUST50", "Y"),
    coordinates = c("x", "y")
    ) {
  
  AK <- variableSelection(
    inFile = inFile,
    incProb = incProb,
    depVar = depVar,
    varsToRemove =  varsToRemove,
    coordinates = coordinates
  )
  
  h <- harmon(
    inFile = inFile,
    incProb = incProb,
    depVar = depVar,
    varsToRemove =  varsToRemove,
    coordinates = coordinates,
    AK = AK
  )
  
  a <- h[[4]]
  
}