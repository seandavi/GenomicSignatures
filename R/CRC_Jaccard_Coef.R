#' CRC Jaccard's Coefficient Table
#'
#' Create Jaccard's Coefficient table represents the binary principle components
#' membership of 8 CRC datasets used in the previous work by Ma *et al.*(2018)
#' (Genome Biol. 2018 Sep 25;19(1):142. doi: 10.1186/s13059-018-1511-4).
#'
#' @param res Output from `kmeans()` clustering
#' @param ds the number of datasets used for clutering
#' @param pc the number of PCs kept from each dataset
#' @return A data frame of Jaccards' coefficients between your dataset and CRC
#' data. Each column represents cluster (crc_cl1 ~ crc_cl4, which represent PCSS1
#' ~ PCSS4) and the rows represent PCs from different datasets used in clustering.
#'
#' @export
CRC_Jaccard_Coef = function(res, ds, pc) {
    num_rows = ds * pc
    x = data.frame(matrix(nrow = num_rows, ncol = 0))
    rownames(x) = names(res$cluster)

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

    return(x)
}
