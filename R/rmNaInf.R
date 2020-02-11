#' Remove missing and Inf values in an expression 
#' 
#' @export
rmNaInf = function(x) {
  
  if (!is.matrix(x)) {
    stop("x must be a matrix (expression values)!")
  }
  
  x = x[apply(x, 1, function(row) {
    !any(is.na(row) | row %in% c(Inf, -Inf))
  }), ]
  return(x)
}