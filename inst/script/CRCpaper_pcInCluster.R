ds = 8
pc = 20

x = data.frame(matrix(nrow = (ds*pc), ncol = 0))

row_names = c()
for (i in c(1:ds)) {
    n = paste0(rep(paste0("ds", i), pc), ".PC",1:pc)
    row_names = c(row_names, n)
}
rownames(x) = row_names

ref_cl = c("ds1.PC1"=1, "ds2.PC2"=1, "ds3.PC1"=1, "ds4.PC1"=1, "ds5.PC2"=1, "ds6.PC1"=1, "ds6.PC2"=1, "ds7.PC1"=1, "ds8.PC7"=1,
           "ds1.PC2"=2, "ds2.PC3"=2, "ds3.PC2"=2, "ds3.PC3"=2, "ds4.PC2"=2, "ds5.PC4"=2, "ds7.PC2"=2, "ds8.PC8"=2,
           "ds1.PC3"=3, "ds2.PC5"=3, "ds3.PC4"=3, "ds4.PC4"=3, "ds6.PC6"=3, "ds7.PC4"=3,
           "ds1.PC4"=4, "ds2.PC4"=4, "ds4.PC3"=4, "ds5.PC1"=4, "ds6.PC3"=4, "ds7.PC3"=4)

for (i in 1:nrow(x)) {
    if (rownames(x)[i] %in% names(ref_cl)[1:9]) {
        x$PCSS1[i] = 1L
    } else {x$PCSS1[i] = 0L}
    if (rownames(x)[i] %in% names(ref_cl)[10:17]) {
        x$PCSS2[i] = 1L
    } else {x$PCSS2[i] = 0L}
    if (rownames(x)[i] %in% names(ref_cl)[18:23]) {
        x$PCSS3[i] = 1L
    } else {x$PCSS3[i] = 0L}
    if (rownames(x)[i] %in% names(ref_cl)[24:29]) {
        x$PCSS4[i] = 1L
    } else {x$PCSS4[i] = 0L}
}

# add the real dataset names instead of ds1-8
new_rownames = lapply(trainingSets_8crc, function(x) {paste0(x, paste0(".PC", 1:20))}) %>% unlist(.)
rownames(x) = new_rownames

# add the number of PCs in each cluster
colnames(x) = paste0(colnames(x), " (", c(9,8,6,6), ")")

write.csv(x, "/data/GenomicSuperSignature/inst/extdata/CRCpaper_pcInCluster.csv")


