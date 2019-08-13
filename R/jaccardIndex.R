#' Create Jaccard Similarity Table
#'
#' Calculate Jaccard Similarity Score between the reference cluster and the testing
#' cluster.
#'
#' @import jaccard
#'
#' @param ref the Jaccard coefficient table of reference clustering
#' @param test the Jaccard coefficient table of test clustering. nrow(y) should be
#' same as nrow(x)
#'
#' @return A data frame of Jaccard similarity score. Each column represents clusters
#' subjected to test and the rows represent clusters from the reference clustering.
#'
#' @export
jaccardIndex = function(ref, test) {

    # test = test[rownames(test) %in% rownames(ref),]
    common_pc = intersect(rownames(test), rownames(ref))
    ref = ref[common_pc,]
    test = test[common_pc,]

    z = cbind(ref, test)
    z = z[which(rowSums(z) != 0),]

    jct = data.frame(matrix(ncol = ncol(test), nrow = ncol(ref)))
    rownames(jct) = colnames(ref)
    colnames(jct) = colnames(test)

    for (i in 1:nrow(jct)) {
        for (j in 1:ncol(jct)) {
            set.seed(123)
            jct[i,j] = jaccard(ref[,i], test[,j])
        }
    }

    return(jct)
}
