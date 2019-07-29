## Top 20 PCs from 8 CRC training datasets
# training data used for `pca.R` script (from Siyuan)
# 9336 overlapping genes based on the paper ~ not sure what's the issue with "MEOX2" gene

load("/data/GenomicSuperSignature/inst/extdata/expr.all.mapped.RData")
pc.num = 20   # the number of top PCs to keep for training datasets

# Very frustratingly, the original datasets (?, .RData files) in this package are not matching
# with the intermediate file (expr.all.mapped.RData)

expr.all.mapped = lapply(expr.all.mapped, GenomicSuperSignature::rmNaInf)

common_gene = lapply(expr.all.mapped, function(x) {rownames(x)}) %>%
    Reduce(intersect, .)

expr.all.mapped = lapply(expr.all.mapped, '[', common_gene, )

# Build a data frame with top pc.num (defaul = 20) PCs of the training datasets
row_num = nrow(expr.all.mapped[[1]])
pc_df = data.frame(matrix(NA, nrow = row_num, ncol = 0))

for (i in seq_along(expr.all.mapped)) {
    pc = prcomp(t(expr.all.mapped[[i]]))$rotation[,1:pc.num]
    colnames(pc) = paste0(names(expr.all.mapped)[i], ".", colnames(pc))
    pc_df = cbind(pc_df, pc)
}

write.csv(pc_df, "/data/GenomicSuperSignature/inst/extdata/PCs_8crc_training.csv")
