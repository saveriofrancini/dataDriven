#' Akaike variable selection.
#'
#' This function performs the variables selection.
#' @param inFile The path for the input dbf file
#' @param depVar The name of the interest variable column
#' @param varsToRemove The name of the columns to do not consider in the analysis (if there are any)
#' @param coordinates The spatial coordinates of the population areal units
#' @return  The output is a list with length 3:
#' (1) v_aux: values of auxiliary variables selected in the sampled spatial units
#' (2) v_aux_U: values of auxiliary variables selected in the sampled spatial units
#' (3) varprov: vector of k+1 elements (the first is = 1 for default, the others are = 1 if the auxiliary variable is selected and 0 otherwise)
#' @keywords Akaike, variables selection
#' @export
#' @examples
#' AK <- variableSelection(
#'   inFile = "data.dbf",
#'   depVar = "fj",
#'   varsToRemove =  c("ID", "IDCLUST50", "Y"),
#'   coordinates = c("x", "y")
#' )

variableSelection <-
  function(inFile,
           depVar,
           cluster,
           varsToRemove,
           coordinates) {
    library(foreign)
    library(dplyr)
    library(Matrix)
    
    db <- read.dbf(inFile)

    # calculate inclusion probabilities
    clusters <- sort(unique(db[, cluster]))
    if (length(clusters) == 1) {
      db$pj <- 1 / nrow(db)
    } else{
      db$pj <- NA
      for (custer_i in clusters) {
        nSamplesInClusterI <- nrow(db[db[, cluster] == custer_i, ])
        db$pj[db[, cluster] == custer_i] <- 1 / nSamplesInClusterI
      }
    }
    
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
    
    Akaike <- function(X_S, X_NS, Y, pj) {
      var_aux_S <- rbind(1, t(X_S))
      var_aux_U <- rbind(1, t(X_NS))
      v_aux <- var_aux_S[1,]
      v_aux_U <- var_aux_U[1,]
      minMRSS <- matrix(NA, nrow(var_aux_S), 1)
      varprov <- matrix(0, nrow(var_aux_S), 1)
      for (i in 1:nrow(var_aux_S)) {
        MRSS <- matrix(NA, nrow(var_aux_S), 1)
        for (j in 2:nrow(var_aux_S)) {
          v_prov <- rbind(v_aux, var_aux_S[j,])
          if (nrow(v_prov) >= ncol(var_aux_S))
            break
          mat_greg <- GREG(v_prov, Y, pj)
          A <- mat_greg [[1]]
          a <- mat_greg [[2]]
          if (det(A) == 0 || rankMatrix(A) < nrow(A))
            MRSS[j] <- NA
          else{
            b_est <- solve(A) %*% a
            somma <- 0
            for (camp in 1:ncol(var_aux_S)) {
              somma <- somma + (Y[camp] - t(b_est) %*% v_prov[, camp]) ^ 2
            }
            MRSS[j] <-
              somma * ((ncol(var_aux_S) + nrow(v_prov)) / ((ncol(
                var_aux_S
              ) - nrow(v_prov))))
          }
        }
        if (nrow(v_prov) >= ncol(var_aux_S) ||
            all(is.na(MRSS) == T))
          break
        minMRSS[i] <- min(MRSS[!is.na(MRSS)])
        if (i > 1 && minMRSS[i - 1] < minMRSS[i])
          break
        varmin <- which.min(MRSS)
        varprov[varmin] <- 1 + varprov[varmin]
        v_aux <- rbind(v_aux, var_aux_S[varmin,])
        v_aux_U <- rbind(v_aux_U, var_aux_U[varmin,])
      }
      varprov[1] <- 1
      output <- list(v_aux, v_aux_U, varprov)
      return(output)
    }
    
    AK <- Akaike(
      X_S = db[is.na(db[, depVar]) == F, colnames(db) %in% c("pj", depVar, varsToRemove, coordinates) == F],
      X_NS = db[is.na(db[, depVar]), colnames(db) %in% c("pj", depVar, varsToRemove, coordinates) == F],
      Y = db[is.na(db[, depVar]) == F, depVar],
      pj = db[is.na(db[, depVar]) == F, "pj"]
    )
    
    return(list(AK, db))
  }
