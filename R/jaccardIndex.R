#' Create Jaccard Similarity Table
#'
#' Calculate Jaccard Similarity Score between the reference cluster and the testing
#' cluster.
#'
#' @import jaccard
#'
#' @param x the Jaccard coefficient table of reference clustering
#' @param y the Jaccard coefficient table of test clustering. nrow(y) should be
#' same as nrow(x)
#'
#' @return A data frame of Jaccard similarity score. Each column represents clusters
#' subjected to test and the rows represent clusters from the reference clustering.
#'
#' @export
jaccardIndex = function(x, y) {
    z = cbind(x, y)
    z = z[which(rowSums(z) != 0),]

    jct = data.frame(matrix(ncol = ncol(y), nrow = ncol(x)))
    rownames(jct) = colnames(x)
    colnames(jct) = colnames(y)

    for (i in 1:nrow(jct)) {
        for (j in 1:ncol(jct)) {
            set.seed(123)
            jct[i,j] = jaccard(x[,i], y[,j])
        }
    }
    return(jct)
}



# n = 4   # the number of cluster from the refernce clustering
#
# jc_tbl = jc_tbl[which((rowSums(jc_tbl) != 0)),]
#
# jc_ind = data.frame(matrix(ncol = (ncol(jc_tbl)-n), nrow = n))
# rownames(jc_ind) = c(paste0("PCSS", 1:n))
# colnames(jc_ind) = colnames(jc_tbl)[(n+1):ncol(jc_tbl)]
#
# for (i in 1:nrow(jc_ind)) {
#   for (j in 1:ncol(jc_ind)) {
#     set.seed(1234)
#     jc_ind[i,j] = jaccard(jc_tbl[,colnames(jc_tbl)[i]],
#                           jc_tbl[,colnames(jc_tbl)[j+n]])
#   }
# }
