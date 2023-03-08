#' Boostrap function.
#'
#' This function applies to bootstrap to construct the error map
#' @param inFile The path for the input dbf file. By default, it is automatically selected a test file provided within the package. To see it use system.file("data", "data.dbf", package = "dataDriven")
#' @param depVar The name of the interest variable column
#' @param varsToRemove The name of the columns to do not consider in the analysis (if there are any)
#' @param coordinates The spatial coordinates of the population areal units
#' @param resampling The number of bootstrap resampling
#' @param cluster The ID strata of the population spatial units. (If simple random sampling is adopted, set cluster = 1)
#' @param outDir The directory to save the output shapefile
#' @return  The output is the predicted map with the associated error map 
#' @keywords Akaike, variables selection, harmonization
#' @export
#' @examples
#' s <- bootStrap()
#'
bootStrap <- function(cluster = "IDCLUST",
                      inFile = system.file("data", "data.dbf", package = "dataDriven"),
                      depVar = "fj",
                      varsToRemove =  c("ID", "Y"),
                      coordinates = c("x", "y"),
                      resampling = 1000,
                      outDir = system.file("data", package = "dataDriven")) {
  library(dataDriven)
  library(doParallel)
  library(foreach)
  library(foreign)
  library(raster)
  
  # identify the sampling used
  db <- read.dbf(inFile)
  db <- db[order(db[, coordinates[1]], db[, coordinates[2]]), ]
  nSamples <- length(db[, depVar][is.na(db[, depVar]) == F])
  clusters <- unique(db[, cluster])
  nClusters <- length(clusters)
  unitsPerCluster <- table(db[, cluster])
  selectedUnits <- lapply(clusters, function(cluster_i) {
    db_i <- db[db[, cluster] == cluster_i,]
    db_i <-
      db_i[order(db_i[, coordinates[1]], db_i[, coordinates[2]]), ]
    which.max(db_i[, depVar])
  })
  selectedUnits <- do.call(c, selectedUnits)
  if (nClusters == 1) {
    sampling <- "random"
  } else{
    if (length(unique(selectedUnits)) == 1) {
      sampling <- "sistematico"
    }
    if (length(unique(selectedUnits)) > 1) {
      sampling <- "opss"
    }
  }
  
  # create initial map
  initialMap <- Main(
    inFile = inFile,
    depVar = depVar,
    cluster = cluster,
    varsToRemove =  c(varsToRemove, cluster),
    coordinates = coordinates
  )
  
  initialMap <- initialMap[order(initialMap$gx, initialMap$gy), ]
  id_originalDb <- paste0(db$x, db$y)
  id_initialMap <- paste0(initialMap$gx, initialMap$gy)
  if (isFALSE(unique(id_originalDb == id_initialMap))) {
    stop("check db sorting")
  }
  
  # perform the bootstrap
  n_cores <- detectCores()
  if (resampling < n_cores)
    n_cores <- n_files
  
  stopImplicitCluster()
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  out <- foreach(
    i = 1:resampling,
    .packages = c("foreign", "dplyr", "Matrix", "dataDriven")
  ) %dopar% {
    db <- read.dbf(inFile)
    db$Harmon <- initialMap$Harmon
    
    if (sampling == "sistematico") {
      sampleFirstCluster <- sample(1:unitsPerCluster[1], 1)
      sampled <- lapply(clusters, function(cluster_i) {
        db_i <- db[db[, cluster] == cluster_i,]
        db_i$Harmon[-sampleFirstCluster] <- NA
        return(db_i)
      })
      # sampled[[1]]$Harmon[sampleFirstCluster]
      db <- do.call(rbind, sampled)
    }
    
    if (sampling == "opss") {
      sampled <- lapply(clusters, function(cluster_i) {
        db_i <- db[db[, cluster] == cluster_i,]
        sampled_i <- sample(x = 1:nrow(db_i), size = 1)
        db_i$Harmon[-sampled_i] <- NA
        return(db_i)
      })
      db <- do.call(rbind, sampled)
    }
    
    if (sampling == "random") {
      sampled <- sample(x = 1:nrow(db), size = nSamples)
      db$Harmon[-sampled] <- NA
    }
    
    file_i <- paste0(tempfile(), i, ".dbf")
    write.dbf(db, file_i)
    rm(db)
    map_i <- Main(
      inFile = file_i,
      depVar = "Harmon",
      cluster = cluster,
      varsToRemove =  c(varsToRemove, depVar, cluster),
      coordinates = coordinates
    )
    
    map_i <- map_i[order(map_i$gx, map_i$gy), ]
    
    return(map_i)
  }
  
  stopCluster(cl)
  gc()
  
  # initialMap[, 1][10]
  cords <- out[[1]][, 1:2]
  
  out <- lapply(out, function(x) {
    (x$Harmon - initialMap$Harmon) ^ 2
  })
  
  p <- do.call(cbind, out)
  p <- apply(p, 1, function(x) {
    sqrt(mean(x))
  })
  initialMap$error <- p
  
  inShp <- paste0(substr(inFile, 1, nchar(inFile) - 4), ".shp")
  inShp <- shapefile(inShp)
  
  inShp <- inShp[order(inShp@data$x, inShp@data$y), ]
  
  inShp@data$error <- initialMap$error
  inShp@data$Harmon <- initialMap$Harmon

  shapefile(inShp, file.path(outDir, "out.shp"), overwrite=TRUE)
  # inShp@data <- inShp@data[,]
  
  return(inShp)
}