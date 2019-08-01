#' Color the matrix in a heatmap-stype
#'
#' @import ComplexHeatmap
#' @import circlize
#' @import grid
#'
#' @param dat a matrix subjected to add background colors in each slot correlated
#' to its value - similar to heatmap
#' @param row_title a character string containing the title of the row
#' @param column_title a character string containing the title of the column
#' @param breaks a numeric vector of length 3. Number represents the values assigned
#' to three colors. Default is `c(0, 0.5, 1)`.
#' @param colors a character vector of length 3. Each represents the color assigned
#' to three breaks. Default is `c("grey", "white", "red")`.
#'
#' @return a heatmap version of the input matrix
#'
#' @export
heatmapTable = function(dat, column_title, row_title, breaks, colors) {

  if (missing(breaks)) {
    breaks = c(0, 0.5, 1)
  } else {breaks = breaks}

  if (missing(colors)) {
    colors = c("grey", "white", "red")
  } else {colors = colors}

  if (missing(column_title)) {
    column_title = character(0)
    hmHeight = nrow(dat)
  } else {
    column_title = column_title
    hmHeight = nrow(dat) + 2
  }

  if (missing(row_title)) {
    row_title = character(0)
    hmWidth = ncol(dat)
  } else {
    row_title = row_title
    hmWidth = ncol(dat) + 2
  }

    ComplexHeatmap::Heatmap(dat, col = colorRamp2(breaks, colors), name = "Corr",
                          cluster_rows = FALSE, cluster_columns = FALSE,
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(sprintf("%.2f", dat[i, j]), x, y, gp = gpar(fontsize = 8))
                          },
                          row_title = row_title,
                          row_title_gp = gpar(fontsize = 11, fontface = "bold"),
                          row_names_gp = gpar(fontsize = 10),
                          row_names_max_width = unit(0.5, "cm"),
                          column_title = column_title,
                          column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                          column_names_side = "bottom",
                          column_names_gp = gpar(fontsize = 10),
                          column_names_rot = 90,
                          column_names_max_height = unit(1, "cm"),
                          heatmap_width = unit(hmWidth, "cm"),
                          heatmap_height = unit(hmHeight, "cm"),
                          show_heatmap_legend = FALSE
  )}
