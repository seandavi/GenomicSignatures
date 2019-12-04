#' Plot heatmap of the sample scores
#'
#' @import ComplexHeatmap
#' @import circlize
#' @import grid
#'
#' @param dat a matrix with samples (row) and PrcompClusters (column)
#' @param dataName a character string containing the name of the dataset
#' @param avgLoading a character string containing the average loading used for scoring
#' @param cluster_rows a logical value, that controls whether to make cluster on rows. Default is TRUE.
#' @param cluster_columns a logical value, that controls whether to make cluster on columns. Default is TRUE.
#' @prarm `...` any additional argument for `ComplexHeatmap::Heatmap`
#'
#' @return a heatmap version of the sample score
#'
#' @export
sampleScoreHeatmap = function(dat, dataName, avgLoading,
                              cluster_rows = TRUE, cluster_columns = TRUE,
                              show_row_names = TRUE,
                              show_column_names = TRUE, ...) {
    ComplexHeatmap::Heatmap(dat,
                            cluster_rows = cluster_rows,
                            cluster_columns = cluster_columns,
                            row_title = dataName,
                            row_title_gp = gpar(fontsize = 11, fontface = "bold"),
                            column_title = avgLoading,
                            column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                            row_names_gp = gpar(fontsize = 0.7),
                            column_names_gp = gpar(fontsize = 5),
                            show_row_names = show_row_names,
                            show_column_names = show_column_names,
                            heatmap_width = unit(10, "cm"),
                            heatmap_legend_param = list(title = "score"), ...)
}
