#' Calculate the number of PC to return
#'
#' This function estimates the number of 'significant' principle components for
#' the SVD decomposition. Output from this function is the minimum `k` for `PLIER::PLIER()`.
#'
#' @param data the same data as to be used for PLIER (z-score recommended) or alternatively the result of an svd calculation
#' @param method Either "eblow" (fast) or "permutation" (slower, but less heuristic)
#' @param B number of permutations
#' @param seed seed for reproducibility
#'
#' @return
#' @example
#'
#' @note This function is adopted from `PLIER::num.pc()`.
#'
#' @export
numpc = function (data, method = "elbow", B = 20, seed = NULL) {

  method = match.arg(method, c("elbow", "permutation"))
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if ((class(data)!="list") & (class(data)!="rsvd")) {
    message("Computing svd")
    n = ncol(data)
    m = nrow(data)
    data = rowNorm(data)

    if (n < 500) {
      k = n
    } else {
      k = max(200, n/4)
    }

    if (k == n) {
      uu = svd(data)
    } else {
      set.seed(123456)
      uu = rsvd(data, k, q = 3)
    }
  }

  else if (!is.null(data[["d"]])) {   # If data is the result from svd calculation
    if (method == "permutation") {
      message("Original data is needed for permutation method.\nSetting method to elbow")
      method = "elbow"
    }
    uu = data
  }


  if (method == "permutation") {
    dstat = uu$d[1:k]^2 / sum(uu$d[1:k]^2)
    dstat0 = matrix(0, nrow = B, ncol = k)

    for (i in 1:B) {
      dat0 = t(apply(data, 1, sample, replace = FALSE))
      if (k == n) {
        uu0 = svd(dat0)
      } else {
        set.seed(123456)
        uu0 = rsvd(dat0,k, q=3)
      }
      dstat0[i,] = uu0$d[1:k]^2 / sum(uu0$d[1:k]^2)
    }

    psv = rep(1, k)

    for (i in 1:k) {
      psv[i] = mean(dstat0[, i] >= dstat[i])
    }

    for (i in 2:k) {
      psv[i] = max(psv[(i - 1)], psv[i])
    }

    nsv = sum(psv <= 0.1)
  }

  else if (method == "elbow") {
    xraw = abs(diff(diff(uu$d)))
    x = smooth(xraw, twiceit = TRUE)
    nsv = which(x <= quantile(x, 0.5))[1] + 1
  }
  return(nsv)
}

