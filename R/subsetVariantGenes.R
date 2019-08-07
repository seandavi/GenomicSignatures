#' Subset highly variant genes
#'
#' @param x a data frame. Row represents different genes.
#' @param percent.var a percentage of top variant genes to subset. Only the multiple
#' of 10 is allowed (maximum = 100).
#'
#' @return a subset of input data frame, which contains only the top `percent.var`
#' percentage of genes.
#'
#' @export
subsetVariantGenes = function(x, percent.var) {
    expression_variance = apply(x, 1, var, na.rm = TRUE)
    quant = quantile(expression_variance, probs = seq(0, 1, 0.1))
    variant_genes = expression_variance > quant[(11-(percent.var*0.1))]
    return(x[variant_genes,])
}
