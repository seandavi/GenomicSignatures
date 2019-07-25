## Top 20 PCs from 10 OVC training datasets
# 10 training datasets from curatedOvarianData
library(curatedOvarianData)
setNames = c("GSE20565_eset", "GSE2109_eset", "GSE26193_eset", "GSE26712_eset",
             "GSE6008_eset", "GSE9891_eset", "PMID17290060_eset", "TCGA_eset")
agilent_setNames = c("GSE17260_eset", "GSE32062.GPL6480_eset")

setNames = c(setNames, agilent_setNames)
data(list = setNames)

# remove/check NA and Inf values
exprs(GSE20565_eset) = exprs(GSE20565_eset) %>% rmNaInf
exprs(GSE2109_eset) = exprs(GSE2109_eset) %>% rmNaInf
exprs(GSE26193_eset) = exprs(GSE26193_eset) %>% rmNaInf
exprs(GSE26712_eset) = exprs(GSE26712_eset) %>% rmNaInf
exprs(GSE6008_eset) = exprs(GSE6008_eset) %>% rmNaInf
exprs(GSE9891_eset) = exprs(GSE9891_eset) %>% rmNaInf
exprs(PMID17290060_eset) = exprs(PMID17290060_eset) %>% rmNaInf
exprs(TCGA_eset) = exprs(TCGA_eset) %>% rmNaInf
exprs(GSE17260_eset) = exprs(GSE17260_eset) %>% rmNaInf 
exprs(GSE32062.GPL6480_eset) = exprs(GSE32062.GPL6480_eset) %>% rmNaInf 

common_gene = Reduce(intersect, list(rownames(exprs(GSE20565_eset)), rownames(exprs(GSE2109_eset)),
                                     rownames(exprs(GSE26193_eset)), rownames(exprs(GSE26712_eset)),
                                     rownames(exprs(GSE6008_eset)), rownames(exprs(GSE9891_eset)),
                                     rownames(exprs(PMID17290060_eset)), rownames(exprs(TCGA_eset)),
                                     rownames(exprs(GSE17260_eset)), rownames(exprs(GSE32062.GPL6480_eset))
))

# List of training datasets subset with common_gene
training_dataset = list()
training_dataset[["GSE20565_eset"]] = GSE20565_eset[rownames(GSE20565_eset) %in% common_gene]
training_dataset[["GSE2109_eset"]] = GSE2109_eset[rownames(GSE2109_eset) %in% common_gene]
training_dataset[["GSE26193_eset"]] = GSE26193_eset[rownames(GSE26193_eset) %in% common_gene]
training_dataset[["GSE26712_eset"]] = GSE26712_eset[rownames(GSE26712_eset) %in% common_gene]
training_dataset[["GSE6008_eset"]] = GSE6008_eset[rownames(GSE6008_eset) %in% common_gene]
training_dataset[["GSE9891_eset"]] = GSE9891_eset[rownames(GSE9891_eset) %in% common_gene]
training_dataset[["PMID17290060_eset"]] = PMID17290060_eset[rownames(PMID17290060_eset) %in% common_gene]
training_dataset[["TCGA_eset"]] = TCGA_eset[rownames(TCGA_eset) %in% common_gene]
training_dataset[["GSE17260_eset"]] = GSE17260_eset[rownames(GSE17260_eset) %in% common_gene]
training_dataset[["GSE32062.GPL6480_eset"]] = GSE32062.GPL6480_eset[rownames(GSE32062.GPL6480_eset) %in% common_gene]

# Build a data frame with top 20 PCs of 10 training datasets
row_num = nrow(exprs(training_dataset[[1]]))
pc_df = data.frame(matrix(NA, nrow = row_num, ncol = 0))

library(GenomicSuperSignature)

for (i in seq_along(training_dataset)) {
  pc = calPC(exprs(training_dataset[[i]]))$rotation[,1:20]
  colnames(pc) = paste0("ds", i, ".", colnames(pc))
  pc_df = cbind(pc_df, pc)
}

# Data frame of 12249 genes x 200 PCs
write.csv(pc_df, "/data/GenomicSuperSignature/inst/extdata/20PCs_10OVCtraining.csv")
