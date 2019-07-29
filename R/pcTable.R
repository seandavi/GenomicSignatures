#' Create a table of PCs of datasets
#'
#' Principle component analysis is performed on each dataset assuming each samples as
#'
#' @import magrittr
#'
#' @param setNames a character vector of name of the datasets.
#' @param pc.num the number of top PCs to keep for each dataset after PCA
#'
#' @return a data frame with the common genes among datasets (row) and the number
#' of PCs defined by pc.num (default = 20) from each dataset (column)
#'
#' @export
pcTable = function (setNames, pc.num = 20, commonGene = NULL) {

  # remove/check NA and Inf values
  training_dataset = list()

  for (set in setNames) {
    dat = get(set)

    # exprs = exprs(dat) %>%
    #   rmNaInf %>%
    #   apply(., 1, function (x) {x - mean(x)}) %>% t

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
    colnames(pc) = paste0("ds", i, ".", colnames(pc))
    pc_df = cbind(pc_df, pc)
  }

  return(pc_df)
}
