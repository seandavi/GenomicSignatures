#' Create a table of PCs of datasets
#'
#' Any gene with NA or -Inf/Inf values is removed, followed by centering. Each
#' datasets are subset with their commone genes.
#' Principle component analysis is performed on each dataset assuming each samples as
#'
#' @import magrittr
#'
#' @param setNames a character vector of name of the datasets.
#' @param pc.num the number of top PCs to keep for each dataset after PCA. Default is 20.
#' @param comonGene a character vector of the common genes. Default is 'NULL',
#' in which case, the common genes of the provided datasets will be used.
#'
#' @return A data frame of principle components (PCs). Rows are common genes of
#' datasets and columns represent PCs (defined by \code{pc.num}) from each dataset.
#'
#' @keyword internal
#'
pcTable = function (setNames, pc.num = 20, commonGene = NULL) {

  # remove/check NA and Inf values
  training_dataset = list()

  for (set in setNames) {
    dat = get(set)
    exprs = exprs(dat)
    exprs = exprs[apply(exprs, 1, function(x){!any(is.na(x)|(x==Inf)|(x==-Inf))}),]
    exprs = apply(exprs, 1, function (x) {x - mean(x)}) %>% t
    new_dat = dat[rownames(dat) %in% rownames(exprs)]
    exprs(new_dat) = exprs
    assign(set, new_dat)
  }

  # common genes of input datasets
  if (!is.null(commonGene)) {
    common_gene = commonGene
  } else {
    common_gene = lapply(setNames,
                         function(x) {
                           dat = get(x)
                           rownames(exprs(dat))
                         }) %>% Reduce(intersect, .)
  }

  # subset with common_gene
  for (set in setNames) {
    dat = get(set)
    training_dataset[[set]] = dat[rownames(dat) %in% common_gene]
  }

  # Build a data frame with top pc.num (defaul = 20) PCs of the training datasets
  row_num = nrow(exprs(training_dataset[[1]]))
  pc_df = data.frame(matrix(NA, nrow = row_num, ncol = 0))

  for (i in seq_along(training_dataset)) {
    count = exprs(training_dataset[[i]])
    pc = prcomp(t(count))$rotation[,1:pc.num]
    colnames(pc) = paste0(names(training_dataset)[i], ".", colnames(pc))
    pc_df = cbind(pc_df, pc)
  }

  return(pc_df)
}
