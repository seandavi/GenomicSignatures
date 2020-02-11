#' Create Jaccard's Coefficient Table
#'
#' @param x a data frame with the common genes (row) and the principle components
#' (PCs) from different datasets (col). The number of column is `the number of PCs
#' x the number of datasets`
#' @param k the number of clusters for `kmeans` clustering if you provide raw data
#' @param seed a random seed
#'
#' @return A data frame of Jaccards' Coefficient. Each column represents cluster
#' and rows represent PCs from different datasets.
#'
#' @export
PCsInCluster = function(x, k, seed = NULL) {
    if (is.null(seed)) {
        set.seed(123)
    } else {set.seed(seed)}

    if (class(x) == "kmeans") {
        res = x
    } else {res = kmeans(t(x), k)}

    k = length(unique(res$cluster))
    z = matrix(0, ncol=k, nrow=sum(res$size))  # sum(res$size) = (# of ds)*(# of pc)

    for(i in seq_along(res$cluster)) {z[i,res$cluster[i]] <- 1}
    colnames(z) = paste0("Cl", k, "_", formatC(1:k, width = 2, format = "d", flag = "0"), " (", res$size, ")")
    rownames(z) = names(res$cluster)

    return(z)
}
