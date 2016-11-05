#' Compute periodically correlated DPCA filters and scores
#'
#' @param X multivariate time series
#' @param PC series of filters returned from pcdpca
#' @keywords pcdpca
#' @export
#' @examples
#' pcdpca.scores(X)
pcdpca.scores <- function(X, XI)
{
  Y.est = list()
  for (d in 1:T){
    Y.est[[d]] = XI[[d]] %c% X  # applies the filter
  }
  Y.res = Y.est[[1]]

  if (T > 1){
    for (d in 2:T){
      idx = T*(0:n) + d
      idx = idx[idx <= n]
      Y.res[idx,] = Y.est[[d]][idx,]
    }
  }
  Y.res
}

