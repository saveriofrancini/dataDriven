#' Main function.
#'
#' This function constructs the harmonized map
#' @param inFile The path for the input dbf file. By default, it is automatically selected a test file provided within the package. To see it use system.file("data", "data.dbf", package = "dataDriven")
#' @param depVar The name of the interest variable column
#' @param varsToRemove The name of the columns to do not consider in the analysis (if there are any)
#' @param coordinates The spatial coordinates of the population areal units
#' @return  The output is a list with length 4:
#' (1) v_aux: matrix of k+1 elements (col 1:name of auxiliary variables, col 2: 1 if the auxiliary variable is selected, 0 otherwise)
#' (2) coeff: number of auxiliary variables selected with the Akaike-type criterion.
#' (3) founded_a: alpha choosen by the LOOCV
#' (4) estim: spatial coordinates of the areal units and the corresponding harmonized values of the interest attribute.
#'            It can be used to plot the harmonized map.
#'
#' @keywords Akaike, variables selection, harmonization
#' @export
#' @examples
#' h <- Main()

Main <- function(inFile = system.file("data", "data.dbf", package = "dataDriven"),
                 depVar = "fj",
                 cluster = "IDCLUST",
                 varsToRemove =  c("ID", "IDCLUST", "Y"),
                 coordinates = c("x", "y")) {
  AK <- variableSelection(
    inFile = inFile,
    depVar = depVar,
    cluster = cluster,
    varsToRemove =  varsToRemove,
    coordinates = coordinates
  )
  
  h <- harmon(
    db = AK[[2]],
    depVar = depVar,
    varsToRemove =  varsToRemove,
    coordinates = coordinates,
    AK = AK[[1]]
  )
  
  return(h[[4]])
  
}