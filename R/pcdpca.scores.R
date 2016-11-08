#' Compute periodically correlated DPCA scores, given the filters XI
#'
#' @param X multivariate time series
#' @param XI series of filters returned from pcdpca
#' @keywords pcdpca
#' @export
#' @examples
#' pcdpca.scores(X)
pcdpca.scores <- function(X, XI)
{
  period = XI$period
  Y = XI %c% pc2stat(X,period=period)
  stat2pc(Y,period=period,n=nrow(X))
}


