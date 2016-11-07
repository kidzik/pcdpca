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
  s = dim(XI[[1]]$operators)[2]
  n = nrow(Y)
  d = ncol(Y)
  T = s / d

  X.est = list()

  Y.est = list()

  ii0 = seq(2,4,by = 2)
  otmp = XI[[1]]$operators
  XI[[1]]$operators[ii0,,] = XI[[2]]$operators[ii0,,]
  XI[[2]]$operators[ii0,,] = otmp[ii0,,]

  for (d in 1:T){
    X.est[[d]] = t(rev(XI[[d]])) %c% Y
  }
  X.res = X.est[[1]]

  if (T>1){
    for (d in 2:T){
      idx = T*(0:n) + d
      idx = idx[idx <= n]
      X.res[idx,] = X.est[[d]][idx,]
    }
  }
  X.res
}
