#' Create a table of PCs of datasets
#'
#' Principle component analysis is performed on each dataset and the table of top
#' principle components (number defined by `pc.num` argument) are saved.
#'
#' @import magrittr
#'
#' @param datasets a list of data frames. Each data frame represent genes x samples of a study.
#' @param pc.num the number of top PCs to keep for each dataset after PCA. Default is 20.
#' @param comonGene a character vector of the common genes. Default is 'NULL',
#' in which case, the common genes of the provided datasets will be used.
#'
#' @return A data frame of principle components (PCs) from each dataset. Rows are
#' common genes from all datasets and columns represent PCs (defined by \code{pc.num})
#' from each dataset.
#'
#' @export
makeTableWithTopPCs = function (datasets, pc.num = 20, commonGene = NULL) {

  # common genes of input datasets
  if (!is.null(commonGene)) {
    common_gene = commonGene
  } else {
    common_gene = lapply(datasets,
                         function(x) {
                           rownames(x)
                         }) %>% Reduce(intersect, .)
  }

  # subset with common_gene
  training_dataset = lapply(datasets, function(x) {x[common_gene,]})

  # Build a data frame with top PCs of the training datasets
  row_num = length(common_gene)
  pc_df = data.frame(matrix(NA, nrow = row_num, ncol = 0))

  for (i in seq_along(training_dataset)) {
    count = training_dataset[[i]]
    pc = prcomp(t(count))$rotation[,1:pc.num]
    colnames(pc) = paste0(names(training_dataset)[i], ".", colnames(pc))
    pc_df = cbind(pc_df, pc)
  }

  return(as.matrix(pc_df))
}
