#' Calculate average loadings of each cluster
#'
#' Name of the input vector is principle component (pc) from different datasets (ds)
#' and the integer is the cluster number assigned to each ds.pc. This function includes
#' kmenas clustering, so `set.seed()` before using this function for reproducible
#' results.
#'
#' @param df A data frame. Each row represents principle components from different
#' training datasets. Each column represents genes used for PCA analysis.
#' @param k The number of clusters used for k-means clustering
#' @param n The number of top principle components from each datasets used for model buildling
#' @param seed a random seed
#'
#' @return A list of kmeans clustering results and a data frame of average loadings.
#' Each column represents cluster and rows represent genes used for PCA.
#'
#' @export
buildAvgLoading = function(df, k, n, seed = NULL) {

    # Kmeans clustering
    if (!is.null(seed)) {
        set.seed(123)
    }

    res = kmeans(df, centers = k)
    x = table(res$cluster) %>% as.data.frame()

    # Separate the PC table into each cluster
    for (i in 1:k) {
        assign(paste0("Cl", k, "_", formatC(i, width = 2, format = "d", flag = "0")),
               df[res$cluster == i,] %>% t)
    }

    # the number of unique datasets in each cluster
    unique_sets = c()
    for (i in 1:k) {
        ind = which(res$cluster == i)
        dataSetName = gsub(".PC\\d+$", "", names(res$cluster)[ind])
        uniqueDataSetName = length(unique(dataSetName))
        unique_sets = c(unique_sets, uniqueDataSetName)
    }

    # Calculate the average of loadings in each cluster
    cl_n = formatC(1:k, width = 2, format = "d", flag = "0")
    cl_ls = mget(paste0("Cl", k, "_", cl_n))
    names(cl_ls) = paste0(names(cl_ls), " (", res$size, "/", unique_sets, ")")
    avg.loadings = data.frame(lapply(cl_ls, function(cl) {apply(cl, 1, mean, na.rm=FALSE)}),
                              row.names = rownames(cl_ls[[1]]), check.names = FALSE)
    rownames(avg.loadings) = colnames(df)

    # Save kmeans clustering outputs and avgloading
    res$avgLoading = avg.loadings
    res$k = k
    res$n = n
    return(res)
}


