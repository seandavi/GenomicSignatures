#' Color the matrix in a heatmap-stype
#'
#' @import ComplexHeatmap
#' @import circlize
#' @import grid
#'
#' @param dat a matrix subjected to add background colors in each slot correlated
#' to its value - similar to heatmap
#' @param title a character string containing the title of the table
#' @param color a numeric vector of length 3. Number represents the values assigned
#' to "grey", "white", and "red", respectively. Default is `c(0, 0.4, 0.8)`.
#'
#' @return a heatmap version of the input matrix
#'
#' @export
heatmapTable = function(dat, title, color) {

    if (missing(color)) {
        color = c(0, 0.4, 0.8)
    } else {color = color}

    ComplexHeatmap::Heatmap(dat, col = colorRamp2(color, c("grey", "white", "red")), name = "Corr",
                          cluster_rows = FALSE, cluster_columns = FALSE,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(sprintf("%.2f", dat[i, j]), x, y, gp = gpar(fontsize = 8))
                          },
                          row_names_gp = gpar(fontsize = 10),
                          row_names_max_width = unit(0.5, "cm"),
                          column_title = title,
                          column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                          column_names_side = "bottom",
                          column_names_gp = gpar(fontsize = 10),
                          column_names_rot = 90,
                          column_names_max_height = unit(1, "cm"),
                          heatmap_width = unit(ncol(dat), "cm"),
                          heatmap_height = unit(nrow(dat) + 2, "cm"),
                          show_heatmap_legend = FALSE
  )}
