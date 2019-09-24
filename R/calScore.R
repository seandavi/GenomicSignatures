#' Calculate the score for a new dataset
#'
#' Rownames of dataset and avg.loadings should be gene symbols.
#'
#' @param dataset A list of SummarizedExperiment (or ExpressionSet) objects.
#' Rownames are in 'symbol' format.
#' @param avg.loadings Output from `avgLoading` function - a data frame of avaerage
#' loadings. Each column represents cluster and rows represent genes used for PCA.
#'
#' @return A list containing the score matrices for input datasets. Scores are
#' assigned to each sample (row) on each cluster (column).
#'
#' @export
calScore = function(dataset, avg.loadings) {
    if (!is.list(dataset)) {dataset = list(dataset)}
    lapply(dataset, function(dat) {
        if (class(dat) == "ExpressionSet") {
            count = exprs(dat)
        } else if (class(dat) == "SummarizedExperiment") {
            count = assay(dat)
        } else if (class(dat) == "matrix") {
            count = dat
        }
        # if (class(dat) == "ExpressionSet") {dat = as(dat, "SummarizedExperiment")}
        # count = assay(dat)
        count = count[apply(count, 1, function(x) {!any(is.na(x) | (x==Inf) | (x==-Inf))}),]
        count = apply(count, 1, function(x) {x - mean(x)}) %>% t
        gene_common = intersect(rownames(avg.loadings), rownames(count))

        score = t(count[gene_common,]) %*% apply(avg.loadings[gene_common,], 2,
                                                 function(x) x / sqrt(sum(x^2, na.rm = TRUE)))
        # CRC paper version
        # score = t(count[gene_common,]) %*% as.matrix(avg.loadings[gene_common,])
        # score = (t(score) / apply(score, 2, sd)) %>% t

        colnames(score) = colnames(avg.loadings)
        return(score)
    })
}
