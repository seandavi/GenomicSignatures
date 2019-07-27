## Top 20 PCs from 8 OVC training datasets
load("/data/trainingSets_8ovc.rda")

for (set in setNames) {
    load(paste0("/data/", set, ".rda"))
}

# common genes of these datasets
common_gene = GenomicSuperSignature::commonGene(setNames)

# remove/check NA and Inf values + subset with common_gene
training_dataset = list()

for (set in setNames) {
    dat = get(set)
    exprs = exprs(dat) %>%
        GenomicSuperSignature::rmNaInf %>%
        apply(., 1, function (x) {x - mean(x)}) %>% t

    exprs(dat) = exprs
    training_dataset[[set]] = dat[rownames(dat) %in% common_gene]
}

# Build a data frame with top 20 PCs of 8 training datasets
row_num = nrow(exprs(training_dataset[[1]]))
pc_df = data.frame(matrix(NA, nrow = row_num, ncol = 0))

for (i in seq_along(training_dataset)) {
    pc = GenomicSuperSignature::calPC(exprs(training_dataset[[i]]))$rotation[,1:20]
    colnames(pc) = paste0("ds", i, ".", colnames(pc))
    pc_df = cbind(pc_df, pc)
}

# Data frame of 13104 genes x 160 PCs
# write.csv(pc_df, "/inst/extdata/20PCs_8OVCtraining.csv")
