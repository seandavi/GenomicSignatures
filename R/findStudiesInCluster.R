#' Find the studies contributing each cluster
#' 
#' @param PCAmodel GenomicSignatures object contianing a PCA model
#' @param ind a numeric vector containing the index of PC clusters you want to find related studies
#' 
findStudiesInCluster = function(PCAmodel, ind = NULL) {
  k = PCAmodel$k   
  z = matrix(0, ncol = k, nrow = sum(PCAmodel$size))
  for (i in seq_along(PCAmodel$cluster)) {
    z[i, PCAmodel$cluster[i]] <- 1
  }
  colnames(z) = paste0("Cl", k, "_", formatC(1:k, width = 2, format = "d", flag = "0"))
  rownames(z) = names(PCAmodel$cluster)

  if (is.null(ind)) {
    return(z)
    stop
  } else {
    for (i in ind) {
      studies = rownames(z[which(z[,i] == 1),]) 
      studies = lapply(studies, function(x) {unlist(strsplit(x, "\\."))[1]})
      studies = unique(unlist(studies))
    }
    return(studies)
  }
}