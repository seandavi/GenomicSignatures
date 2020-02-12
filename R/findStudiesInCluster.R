#' Find the studies contributing each cluster
#'
#' @param PCAmodel GenomicSignatures object contianing a PCA model
#' @param ind a numeric vector containing the index of PC clusters you want to
#' find related studies. Default is `NULL`.
#'
#' @return If `ind` is not specified, this function will create a binary membership
#' matrix with the PCs in the raw and the clusters in the column. If `ind` is specificed,
#' the ouput will be a list of studies in the specified cluster.
#'
#' @export
findStudiesInCluster = function(PCAmodel, ind = NULL) {
  k = PCAmodel$k
  z = matrix(0, ncol = k, nrow = sum(PCAmodel$size))
  for (i in seq_along(PCAmodel$cluster)) {
    z[i, PCAmodel$cluster[i]] <- 1
  }
  colnames(z) = paste0("Cl", k, "_", formatC(1:k, width = 2, format = "d", flag = "0"))
  rownames(z) = names(PCAmodel$cluster)

  if (is.null(ind)) {
      studies = list()
      for (i in 1:ncol(z)) {
          studies[[i]] = rownames(z[which(z[,i] == 1),])
          studies[[i]] = lapply(studies[[i]], function(x) {unlist(strsplit(x, "\\."))[1]})
          studies[[i]] = unique(unlist(studies[[i]]))
          names(studies)[i] = colnames(z)[i]
      }
  } else {
    for (i in ind) {
        studies = rownames(z[which(z[,i] == 1),])
        studies = lapply(studies, function(x) {unlist(strsplit(x, "\\."))[1]})
        studies = unique(unlist(studies))
    }
  }
  return(studies)
}
