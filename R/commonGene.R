#' Common genes from datasets
#'
#' @import tidyr
#'
#' @param setNames a character vector. Name of the datasets.
#' @return a character vector with the common genes of the datasets.
#' @export
commonGene = function(setNames) {
    lapply(setNames,
           function(x) {
               dat = get(x)
               rownames(exprs(dat))
           }) %>% Reduce(intersect, .)
}
