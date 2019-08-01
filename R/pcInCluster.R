#' Create Jaccard's Coefficient Table
#'
#' @param x a data frame with the common genes (row) and the principle components
#' (PCs) from different datasets (col). The number of column is `the number of PCs
#' x the number of datasets`
#' @param k the number of clusters (k) for `Kmeans` clustering
#'
#' @return A data frame of Jaccards' Coefficient. Each column represents cluster
#' and rows represent PCs from different datasets.
#'
#' @export
pcInCluster = function(x, k) {
    set.seed(123)
    res = kmeans(t(x), k)

    k = length(unique(res$cluster))
    z = matrix(0, ncol=k, nrow=sum(res$size))  # sum(res$size) = (# of ds)*(# of pc)

    for(i in seq_along(res$cluster)) {z[i,res$cluster[i]] <- 1}
    colnames(z) = paste0("k", k, "_", formatC(1:k, width = 2, format = "d", flag = "0"), " (", res$size, ")")
    rownames(z) = names(res$cluster)

    return(z)
}
