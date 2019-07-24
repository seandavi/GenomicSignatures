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
    set.seed(123)
    res = kmeans(df, centers = k)

    x = table(res$cluster) %>% as.data.frame()
    cl_n = formatC(1:nrow(x), width = 2, format = "d", flag = "0")
    names(x)[1] = "cl"
    x$cl = c(paste0("Cl_", cl_n))

    for (i in 1:nrow(x)) {
        assign(paste0("Cl", nrow(x), "_", formatC(i, width = 2, format = "d", flag = "0")),
               df[res$cluster == i,] %>% t)
    }

    cl_ls = mget(paste0("Cl", nrow(x), "_", cl_n))
    avg.loadings = data.frame(lapply(cl_ls, function(cl) {apply(cl, 1, mean, na.rm=FALSE)}),
                              row.names = rownames(cl_ls[[1]]))
    return(avg.loadings)
}


