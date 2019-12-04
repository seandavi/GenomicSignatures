#' Calculate average loadings of each cluster
#'
#' Name of the input vector is principle component (pc) from different datasets (ds)
#' and the integer is the cluster number assigned to each ds.pc.
#'
#' @param df A data frame. Each row represents principle components from different
#' training datasets. Each column represents genes used for PCA analysis.
#' @param k The number of clusters used for k-means clustering
#' @return A data frame of average loadings. Each column represents cluster
#' and rows represent genes used for PCA.
#' @export
avgLoading = function(df, k) {

    # Kmeans clustering
    set.seed(123)
    res = kmeans(df, centers = k)
    x = table(res$cluster) %>% as.data.frame()

    # Separate the PC table into each cluster
    for (i in 1:nrow(x)) {
        assign(paste0("Cl", nrow(x), "_", formatC(i, width = 2, format = "d", flag = "0")),
               df[res$cluster == i,] %>% t)
    }

    # the number of unique datasets in each cluster
    unique_sets = c()
    for (i in 1:nrow(x)) {
        ind = which(res$cluster == i)
        dataSetName = gsub("_LV\\d+$", "", names(res$cluster)[ind])
        uniqueDataSetName = length(unique(dataSetName))
        unique_sets = c(unique_sets, uniqueDataSetName)
    }

    # Calculate the average of gene expressions in each cluster
    cl_n = formatC(1:nrow(x), width = 2, format = "d", flag = "0")
    cl_ls = mget(paste0("Cl", nrow(x), "_", cl_n))
    names(cl_ls) = paste0(names(cl_ls), " (", res$size, "/", unique_sets, ")")
    avg.loadings = data.frame(lapply(cl_ls, function(cl) {apply(cl, 1, mean, na.rm=FALSE)}),
                              row.names = rownames(cl_ls[[1]]), check.names = FALSE)
    rownames(avg.loadings) = colnames(df)
    return(avg.loadings)
}


