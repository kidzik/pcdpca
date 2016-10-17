#' For given scores \code{Y} and dynamic principal components \code{XI}
#' retrive a series from which scores \code{Y} were calculated.
#' This procedure should be seen as the inverse of \code{\link{pcdpca.scores}}.
#'
#' @title Retrieve a process from given scores
#' @param Y scores process
#' @param XI principal components series
#' @return Retrived process X
#' @seealso \code{\link{pcdpca.scores}}, \code{\link{pcdpca}}
#' @references Kidzinski, Kokoszka, Jouzdani
#' Dynamic principal components of periodically correlated functional time series
#' Research report, 2016
#' @export
pcdpca.inverse = function(Y,XI){
  T = length(XI)
  n = nrow(Y)

  X.est = list()

  Y.est = list()
  for (d in 1:T){
    X.est[[d]] = t(rev(XI[[1]])) %c% Y
  }
  X.res = X.est[[1]]

  for (d in 2:T){
    idx = T*(0:n) + d
    idx = idx[idx <= n]
    X.res[idx,] = X.est[[d]][idx,]
  }
  X.res
}
