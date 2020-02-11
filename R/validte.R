.loadingCor = function(dat, avgLoading) {

    if (class(dat) == "ExpressionSet") {
        count = exprs(dat)
    } else if (class(dat) == "SummarizedExperiment") {
        count = assay(dat)
    } else if (class(dat) == "matrix") {
        count = dat
    }

    count = count[apply(count, 1, function (x) {!any(is.na(x) | (x==Inf) | (x==-Inf))}),]

    gene_common = intersect(rownames(avgLoading), rownames(count))
    prcomRes = prcomp(t(count[gene_common,]))
    loadings = prcomRes$rotation[, 1:8]
    loading_cor = abs(cor(avgLoading[gene_common,], loadings[gene_common,],
                          use = "pairwise.complete.obs",
                          method = "pearson"))
    return(loading_cor)
}


#' Validatin of a new dataset
#'
#' @importFrom SummarizedExperiment assay
#'
#' @param dataset Single or a list of SummarizedExperiment (or ExpressionSet) object(s).
#' Rownames are in 'symbol' format.
#' @param avgLoading output from \code{avgLoading} function - a data frame of avaerage
#' loadings. Each column represents cluster and rows represent genes used for PCA.
#' @param level Defibe how to ouput validation in two different forms, \code{c("max", "all")}.
#' Default is "max", which outputs the matrix containing only the maximum coefficient.
#' To get the coefficient of all 8 PCs, set this argument as "all".
#'
#' @return A matrix containing the maximum pearson correlation coefficient between
#' the top 8 PCs of the dataset(s) and pre-calculated average loadings of training
#' datasets. Each row represents a new dataset for test, and each column represents
#' clusters from training datasets. If \code{level = "all"}, a list containing the matrices
#' of the pearson correlation coefficient between all top 8 PCs of the dataset(s) and
#' avgLoading.
#'
#' @export
validate = function(dataset, avgLoading, level = "max") {
    if (!is.list(dataset)) {
        x = .loadingCor(dataset, avgLoading)
        z = apply(x, 1, max)
        return(t(z))
    } else {
        x = lapply(dataset, .loadingCor, avgLoading)
        if (level == "max") {
            z = sapply(x, function(y) {apply(y, 1, max)})
            return(t(z))
        } else if (level == "all") {
            return(t(x))
        }
    }
}
