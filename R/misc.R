#' Misc Functions for GenomicSuperSignature package
#'
#' Remove missing / Inf values in an expression
#'
#' @param x a matrix with the expression values
#' @return the updated input matrix where NA and Inf values are removed
#' @export
rmNaInf <- function(x) {
  if(!is.matrix(x)) stop("x must be a matrix (expression values)!")
  x <- x[apply(x, 1, function(row) {
    !any(is.na(row) | row %in% c(Inf, -Inf))
  }), ]
  return(x)
}

# gg color function
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }

#' Filling in zero values for classification
#'
#' @param exprs an ExpressionSet object
#' @param names gene names
#' @export
fillInZero <- function(exprs, names) {
  if (!all(rownames(exprs) %in% names)) stop("Gene names don't match!")
  exprs_new <- matrix(0, length(names), ncol(exprs))
  dimnames(exprs_new) <- list(names, colnames(exprs))
  exprs_new[rownames(exprs), ] <- exprs
  return(exprs_new)
}

#' Calculate silhouette width
#'
#' @param distance a matrix like object of distance
#' @param lvls levels
#' @return data frame with claculated silhouette width
#' @export
calcSilWidth <- function(distance, lvls) {
  lvl.factor <- factor(lvls)
  distance.mat <- as.matrix(distance)
  df.return <- data.frame(lvl1 = NULL, lvl2 = NULL, avg.sil.width = NULL)
  for(i in 1:(nlevels(lvl.factor) - 1)) {
    lvl1 <- levels(lvl.factor)[i]
    if(sum(lvl.factor %in% lvl1) == 0) next
    for(j in (i + 1):nlevels(lvl.factor)) {
      lvl2 <- levels(lvl.factor)[j]
      if(sum(lvl.factor %in% lvl2) == 0) next
      distance.mat.tmp <- distance.mat[lvl.factor %in% c(lvl1, lvl2),
                                       lvl.factor %in% c(lvl1, lvl2)]
      lvl.factor.tmp <- lvl.factor[lvl.factor %in% c(lvl1, lvl2)]
      avg.sil.width <- summary(silhouette(lvl.factor.tmp %>% as.numeric,
                                          dist=distance.mat.tmp %>% as.dist) )$avg.width
      df.return <- rbind(df.return,
                         data.frame(lvl1 = lvl1, lvl2 = lvl2, avg.sil.width = avg.sil.width))
    }
  }
  avg.sil.width <- summary(silhouette(lvl.factor %>% as.numeric,
                                      dist=distance.mat %>% as.dist) )$avg.width
  df.return <- rbind(df.return,
                     data.frame(lvl1 = "all", lvl2 = "all", avg.sil.width = avg.sil.width))
  return(df.return)
}
