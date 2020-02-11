#' Plot heatmap of the sample scores
#'
#' @import ComplexHeatmap
#' @importfrom circlize colorRamp2
#' @importfrom grid gpar
#'
#' @param score a matrix with samples (row) and PrcompClusters (column). If it
#' is a simple vector, it will be converted to a one-column matrix.
#' @param dataName a character string containing the name of the dataset
#' @param modelName a character string containing the average loading used for scoring
#' @param cluster_rows a logical value, that controls whether to make cluster on rows. Default is TRUE.
#' @param cluster_columns a logical value, that controls whether to make cluster on columns. Default is TRUE.
#' @prarm `...` any additional argument for `ComplexHeatmap::Heatmap`
#'
#' @return a heatmap version of the sample score
#'
#' @export
sampleScoreHeatmap = function(score, dataName, modelName,
                              cluster_rows = TRUE,
                              cluster_columns = TRUE,
                              show_row_names = TRUE,
                              show_column_names = TRUE,
                              row_names_gp = 0.7,
                              column_names_gp = 5, ...) {
    ComplexHeatmap::Heatmap(score,
                            cluster_rows = cluster_rows,
                            cluster_columns = cluster_columns,
                            row_title = dataName,
                            row_title_gp = gpar(fontsize = 11, fontface = "bold"),
                            column_title = modelName,
                            column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                            row_names_gp = gpar(fontsize = row_names_gp),
                            column_names_gp = gpar(fontsize = column_names_gp),
                            show_row_names = show_row_names,
                            show_column_names = show_column_names,
                            heatmap_width = unit(10, "cm"),
                            heatmap_legend_param = list(title = "score"), ...)
}
