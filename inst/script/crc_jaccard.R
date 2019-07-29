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
        x$crc_cl1[i] = 1L
    } else {x$crc_cl1[i] = 0L}
    if (rownames(x)[i] %in% names(ref_cl)[10:17]) {
        x$crc_cl2[i] = 1L
    } else {x$crc_cl2[i] = 0L}
    if (rownames(x)[i] %in% names(ref_cl)[18:23]) {
        x$crc_cl3[i] = 1L
    } else {x$crc_cl3[i] = 0L}
    if (rownames(x)[i] %in% names(ref_cl)[24:29]) {
        x$crc_cl4[i] = 1L
    } else {x$crc_cl4[i] = 0L}
}

write.csv(x, "/inst/data/crc_jaccard.csv")
