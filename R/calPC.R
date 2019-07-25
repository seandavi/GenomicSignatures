#' PCA on gene expression profile
#'
#' @param x Numeric data matrix or numeric data frame: each row corresponds to genes, and each column corresponds to samples of the study.
#' @return results of principle components analysis (PCA) in an object of class \code{\link[stats]{prcomp}}
#' @export
calPC = function(x) {
    dat.pca = prcomp(t(x))
    return(dat.pca)
}
