# Load OVC datasets
library(curatedOvarianData)
load("data/trainingSets_10ovc.rda")
data(list = trainingSets_10ovc)

pc_df = GenomicSuperSignature::pcTable(setNames)
devtools::use_data(pc_df, PCs_10ovc_training, overwrite = TRUE)

# write.csv(pc_df, "inst/extdata/PCs_10ovc_training.csv")

# # common genes of these datasets
# common_gene = GenomicSuperSignature::commonGene(setNames)
#
# # remove/check NA and Inf values + subset with common_gene
# training_dataset = list()
#
# for (set in setNames) {
#     dat = get(set)
#     exprs = exprs(dat) %>%
#         GenomicSuperSignature::rmNaInf %>%
#         apply(., 1, function (x) {x - mean(x)}) %>% t
#
#     training_dataset[[set]] = dat[rownames(dat) %in% common_gene]
# }
#
# # Build a data frame with top 20 PCs of 8 training datasets
# row_num = nrow(exprs(training_dataset[[1]]))
# pc_df = data.frame(matrix(NA, nrow = row_num, ncol = 0))
#
# for (i in seq_along(training_dataset)) {
#     pc = GenomicSuperSignature::calPC(exprs(training_dataset[[i]]))$rotation[,1:20]
#     colnames(pc) = paste0("ds", i, ".", colnames(pc))
#     pc_df = cbind(pc_df, pc)
# }
#
# # Data frame of 12249 genes x 200 PCs
# # write.csv(pc_df, "/inst/extdata/20PCs_10OVCtraining.csv")
