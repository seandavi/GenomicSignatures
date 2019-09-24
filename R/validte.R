#' Validatin of a new dataset
#'
#' @importFrom SummarizedExperiment assay
#'
#' @param dataset A list of SummarizedExperiment (or ExpressionSet) objects.
#' Rownames are in 'symbol' format.
#' @param avg.loadings output from `avgLoading` function - a data frame of avaerage
#' loadings. Each column represents cluster and rows represent genes used for PCA.
#'
#' @return A matrix containing the maximum pearson correlation coefficiency between
#' the top 8 PCs of the dataset(s) and pre-calculated average loadings of training
#' datasets. Each row represents a new dataset for test, and each column represents
#' clusters from training datasets.
#'
#' @export
validate = function(dataset, avg.loadings) {
  t(sapply(dataset, function(dat) {
    if (class(dat) == "ExpressionSet") {dat = as(dat, "SummarizedExperiment")}
    count = assay(dat)
    count = count[apply(count, 1, function (x) {!any(is.na(x) | (x==Inf) | (x==-Inf))}),]

    gene_common = intersect(rownames(avg.loadings), rownames(count))
    prcomRes = prcomp(t(count[gene_common,]))
    loadings = prcomRes$rotation[, 1:8]
    loading_cor = abs(cor(avg.loadings[gene_common,], loadings[gene_common,],
                          use = "pairwise.complete.obs",
                          method = "pearson"))
    return(apply(loading_cor, 1, max))
  }))
}
