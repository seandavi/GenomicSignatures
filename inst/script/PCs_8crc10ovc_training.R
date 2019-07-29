## Top 20 PCs from 8 CRC training datasets
load("/data/GenomicSuperSignature/inst/extdata/expr.all.mapped.RData")
pc.num = 20   # the number of top PCs to keep for training datasets

expr.all.mapped = lapply(expr.all.mapped, GenomicSuperSignature::rmNaInf)
common_gene = lapply(expr.all.mapped, function(x) {rownames(x)}) %>%
    Reduce(intersect, .)

##### Add OVC training datasets ################################################
# Pre-processing of 8 CRC training datasets is complicating. So I just added this
# part to the existing script (/inst/script/PCs_8crc_training.R).
load("/data/GenomicSuperSignature/data/trainingSets_10ovc.rda")
setNames = trainingSets_10ovc

for (set in setNames) {
    load(paste0("/data/GenomicSuperSignature/data/", set, ".rda"))
    dat = get(set)
    exprs = exprs(dat)
    exprs = exprs[apply(exprs, 1, function(x){!any(is.na(x)|(x==Inf)|(x==-Inf))}),]
    exprs = apply(exprs, 1, function (x) {x - mean(x)}) %>% t
    new_dat = dat[rownames(dat) %in% rownames(exprs)]
    exprs(new_dat) = exprs
    assign(set, new_dat)
}

# common genes of input datasets
common_gene_ovc = lapply(setNames,
                     function(x) {
                         dat = get(x)
                         rownames(exprs(dat))
                     }) %>% Reduce(intersect, .)

# combine CRC and OVC training datasets
common_gene = intersect(common_gene, common_gene_ovc)
expr.all.mapped = lapply(expr.all.mapped, '[', common_gene, )
CrcOv_trainingDatasets = list()

# subset with common_gene
for (set in setNames) {
    dat = get(set)
    CrcOv_trainingDatasets[[set]] = exprs(dat[rownames(dat) %in% common_gene])
}

CrcOv_trainingDatasets = c(CrcOv_trainingDatasets, expr.all.mapped)

# Build a data frame with top pc.num (defaul = 20) PCs of the training datasets
row_num = nrow(CrcOv_trainingDatasets[[1]])
pc_df = data.frame(matrix(NA, nrow = row_num, ncol = 0))

for (i in seq_along(CrcOv_trainingDatasets)) {
    pc = prcomp(t(CrcOv_trainingDatasets[[i]]))$rotation[,1:pc.num]
    colnames(pc) = paste0(names(CrcOv_trainingDatasets)[i], ".", colnames(pc))
    pc_df = cbind(pc_df, pc)
}

write.csv(pc_df, "/data/GenomicSuperSignature/inst/extdata/PCs_8crc10ov_training.csv")
