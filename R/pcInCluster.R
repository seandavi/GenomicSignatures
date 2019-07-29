#' Create Jaccard's Coefficient Table
#'
#' @param res output from `kmeans()` clustering
#'
#' @return A data frame of Jaccards' Coefficient. Each column represents cluster
#' and rows represent PCs from different datasets.
#'
#' @export
pcInCluster = function(res) {
    k = length(unique(res$cluster))
    z = matrix(0, ncol=k, nrow=sum(res$size))  # sum(res$size) = (# of ds)*(# of pc)

    for(i in seq_along(res$cluster)) {z[i,res$cluster[i]] <- 1}
    colnames(z) = paste0("k", k, "_", 1:k)
    rownames(z) = names(res$cluster)

    return(z)
}
