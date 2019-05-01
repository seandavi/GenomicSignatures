#' Calculate Principle Components
#'
#' @param x gene expression count matrix whit genes and samples of the study
#' @return results of principle components analysis in an object of class \code{\link[stats]{prcomp}}
#' @export
calPC = function(x) {
    dat.pca = prcomp(t(x))
    return(dat.pca)
}
