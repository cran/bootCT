#' Create matrix of lagged variables
#'
#' This function lags a set of variables in a matrix, each with a separate index. It is also possible to retain only the last lag order.
#'
#' @param X numeric matrix whose columns are subject to lagging
#' @param k vector of lag orders
#' @param last.only If TRUE only the k-th order lag will be computed, otherwise all lags from 1 to k
#'
#' @return a matrix whose columns are the original variables and the k-th order lagged variables. Column name suffix ".lx".
#' @examples
#'
#' data(ger_macro)
#'
#' lag_mts(X = ger_macro, k = 3, last.only = FALSE)
#'
#' @export
lag_mts =
  function(X,
           k,
           last.only = F) {

    X = as.data.frame(X)
    n = nrow(X)
    d = ncol(X)
    if(length(k) < d){
      warning("length of lags mismatched, taking the first element")
      k = rep(k[1],d)
    }
    if(length(k) == 1){
      k = rep(k[1],d)
    }
    na.head = list()
    new.ts = list()
    cnames = list()
    t = 1
    if (!last.only) {
      for (j in 1:d) {
        if (k[j] > 0) {
          if (last.only) {
            lk = k[j]
          } else{
            lk = 1:k[j]
          }
          for (i in 1:k[j]) {
            na.head[[t]] = rep(NA, lk[i])
            new.ts[[t]] = c(na.head[[t]], X[(1:(n - lk[i])), j])
            cnames[[t]] = paste0(colnames(X)[j], ".l", lk[i], sep = "")
            t = t + 1
          }
        }
      }
    } else{
      for (j in 1:d) {
        na.head[[t]] = rep(NA, k[j])
        new.ts[[t]] = c(na.head[[t]], X[(1:(n - k[j])), j])
        cnames[[t]] = paste0(colnames(X)[j], ".l", k[j], sep = "")
        t = t + 1
      }
    }
    XL = data.frame(do.call(cbind, new.ts))
    colnames(XL) = unlist(cnames)
    return(XL)
  }
