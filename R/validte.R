#' Validatin of a new dataset
#'
#' @importFrom SummarizedExperiment assay
#'
#' @param dataset A list of SummarizedExperiment (or ExpressionSet) objects.
#' Rownames are in 'symbol' format.
#' @param avg.loadings output from \code{avgLoading} function - a data frame of avaerage
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
#' avg.loadings.
#'
#' @export
validate = function(dataset, avg.loadings, level = "max") {
  x = lapply(dataset, function(dat) {

      if (class(dat) == "ExpressionSet") {
          count = exprs(dat)
      } else if (class(dat) == "SummarizedExperiment") {
          count = assay(dat)
      } else if (class(dat) == "matrix") {
          count = dat
      }

    count = count[apply(count, 1, function (x) {!any(is.na(x) | (x==Inf) | (x==-Inf))}),]

    gene_common = intersect(rownames(avg.loadings), rownames(count))
    prcomRes = prcomp(t(count[gene_common,]))
    loadings = prcomRes$rotation[, 1:8]
    loading_cor = abs(cor(avg.loadings[gene_common,], loadings[gene_common,],
                          use = "pairwise.complete.obs",
                          method = "pearson"))
    return(loading_cor)
  })

  if (level == "max") {
    z = sapply(x, function(y) {apply(y, 1, max)})
    return(t(z))
  } else if (level == "all") {
    return(t(x))
  }
}




