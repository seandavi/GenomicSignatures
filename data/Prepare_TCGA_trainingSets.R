# Pre-processing of TCGA RNAseq datasets

### 4 cancer datasets from curatedTCGAData
load("/data/GenomicSuperSignature/data/testSets_tcga.rda")
setNames = testSets_tcga
dataset = list()

for (set in setNames) {
    load(paste0("/data/GenomicSuperSignature/data/", set, ".rda"))
    dataset[[set]] = get(set) %>%
        EnrichmentBrowser::idMap(., org = "hsa", from = "ENTREZID", to = "SYMBOL")

    # log2 transformation of raw counts
    assay(dataset[[set]]) = log2(assay(dataset[[set]]) + 1)
}
names(dataset) = c("BRCA", "COAD", "LUAD", "READ")

### TCGA_OV RNAseq data from curatedOvarianData
load("/data/GenomicSuperSignature/data/TCGA.RNASeqV2_eset.rda")
x = as(TCGA.RNASeqV2_eset, "SummarizedExperiment")

# assay(x) = assay(x)[rowSums(assay(x)) > 2,]
rs <- rowSums(assay(x) > 2)
keep <-  rs >= ncol(x) / 2
x <- x[keep,]
dataset[["OV"]] = x

TCGA_trainingSets = dataset
save(TCGA_trainingSets, file = "/data/GenomicSuperSignature/data/TCGA_trainingSets.rda")
