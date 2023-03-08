#' Harmonized map construction.
#'
#' This function constructs the harmonized map
#' @param inFile The path for the input dbf file
#' @param incProb The name of first-order inclusion probability column
#' @param depVar The name of the interest variable column
#' @param varsToRemove The name of the columns to do not consider in the analysis (if there are any)
#' @param coordinates The spatial coordinates of the population areal units
#' @return  The output is a list with length 4:
#' (2) coeff: number of auxiliary variables selected with the Akaike-type criterion.
#' (3) founded_a: alpha choosen by the LOOCV
#' (4) estim: spatial coordinates of the areal units and the corresponding harmonized values of the interest attribute.
#'            It can be used to plot the harmonized map.
#' @keywords Akaike, variables selection, harmonization
#' @export
#' @examples
#' h <- harmon(
#'   inFile = "data.dbf",
#'   incProb = "pj",
#'   depVar = "fj",
#'   varsToRemove =  c("ID", "IDCLUST50", "Y"),
#'   coordinates = c("x", "y"),
#'   AK = AK
#' )

harmon <-
  function(db,
           depVar,
           varsToRemove,
           coordinates,
           AK) {
    library(foreign)
    library(dplyr)
    library(Matrix)
    
    GREG <- function(X, Y, pj) {
      ncamp <- ncol(X)
      A <- a <- 0
      for (i in 1:ncamp) {
        A <- A + ((X[, i] %*% t(X[, i])) / pj[i]) %>% as.matrix()
        a <- a + ((Y[i] * X[, i]) / pj[i]) %>% as.vector()
      }
      output <- list(A, a)
      return(output)
    }
    dist2D <- function(x1, y1, x2 , y2) {
      dist = sqrt((y2 - y1) ^ 2 + (x2 - x1) ^ 2)
      return(dist)
    }
    
    vaux1 <- AK[[1]]
    vaux1U <- AK[[2]]
    X_S = db[is.na(db[, depVar]) == F, colnames(db) %in% c("pj", depVar, varsToRemove, coordinates) == F]
    v_aux <- cbind(c("Intercept", names(X_S)), AK[[3]])
    coeff <- nrow(vaux1)
    Y = db[is.na(db[, depVar]) == F, depVar]
    pj = db[is.na(db[, depVar]) == F, "pj"]
    mat_greg <- GREG(vaux1, Y, pj)
    b_est <- solve(mat_greg[[1]]) %*% mat_greg[[2]]
    S <- db[is.na(db[, depVar]) == F, coordinates]
    U_NS <- db[is.na(db[, depVar]), coordinates]
    
    for (sy in 1:length(Y)) {
      S$ej[sy] <- Y[sy] - t(b_est) %*% vaux1[, sy]
    }
    Tx <- apply(vaux1, 1, sum) + apply(vaux1U, 1, sum)
    sumT <- 0
    for (sb2 in 1:length(Y)) {
      sumT <- sumT + (S$ej[sb2] / pj[sb2])
    }
    Treg <- (t(b_est) %*% Tx + sumT) %>% as.numeric()
    
    ## alpha loocv
    colnames(S) <- c("gx", "gy", "ej")
    DistListS <- lapply(1:nrow(S), function(s) {
      dist2D(S$gx, S$gy, S$gx[s], S$gy[s])
    })
    dist_s1 <- matrix(unlist(DistListS), ncol = nrow(S))
    
    alphaoptimal <- matrix(NA, 19, 1)
    for (alp in 3:21) {
      wij <-
        (dist_s1 ^ -alp) / (apply(
          dist_s1 ^ -alp,
          1,
          FUN = function(x) {
            sum(x[x != Inf])
          }
        ))
      sumtemp <- 0
      for (mi in 1:length(Y)) {
        wijtemp <- wij[mi,-mi]
        Stemp <- S[-mi,]
        error <- wijtemp %*% Stemp$ej
        sumtemp <- sumtemp + (S$ej[mi] - error) ^ 2
      }
      alphaoptimal[alp - 2] <- sumtemp
    }
    opt_alfa <- which.min(alphaoptimal) + 2
    
    colnames(U_NS) <- c("gx", "gy")
    DistList <- lapply(1:nrow(S), function(s) {
      dist2D(U_NS$gx, U_NS$gy, S$gx[s], S$gy[s])
    })
    dist_s <- matrix(unlist(DistList), ncol = nrow(S))
    if (opt_alfa == 21) {
      ###more than one min
      minimum <- apply(dist_s, 1, min)
      M <- NULL
      for (mi in 1:nrow(dist_s)) {
        M[mi] <- list(which(dist_s[mi,] %in% minimum[mi]))
        U_NS$ej[mi] <- (1 / length(M[[mi]])) * sum(S$ej[M[[mi]]])
      }
      ###
    }
    if (opt_alfa != 21) {
      wij <- (dist_s ^ -opt_alfa) / (apply(dist_s ^ -opt_alfa, 1, sum))
      U_NS$ej <- (wij %*% S$ej)
    }
    
    founded_a <- opt_alfa
    for (su in 1:nrow(U_NS)) {
      U_NS$Y[su] <- t(b_est) %*% vaux1U[, su] + U_NS$ej[su]
    }
    U_NS$Y[U_NS$Y < 0] <- 0
    S$Y <- Y
    estim <- rbind(S, U_NS)
    ##Harmonization
    Tidw <- sum(estim$Y)
    estim$Harmon <- estim$Y * (Treg / Tidw)
    estim <- estim[, c("gx", "gy", "Harmon")]
    output <- list(v_aux, coeff, founded_a, estim)
    return(output)
  }