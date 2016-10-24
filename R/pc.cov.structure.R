#' Given two multivariate periodically correlated time series \eqn{X_t} and \eqn{Y_t} with period \eqn{T}
#' estimates the covariances of \eqn{Cov(X_{d+k} Y_d)} for \eqn{k \in [-q,q]} and some index \eqn{d}.
#' \code{\link{lagged.cov}} is used for the estimation at each lag.
#'
#' @title Estimate the periodically correlated covariance structure within a given window \eqn{k \in [-q,q]}
#' @param X first process
#' @param Y second process, if null then autocovariance of \code{X} is computed
#' @param q size of the window (covariances from \code{-q} to \code{q} will be computed)
#' @param d index in period between \eqn{0} and \eqn{T-1}
#' @param T period
#' @return a time domain operator
#' @export
#' @examples
#' X = rar(100)
#' Y = rar(100)
#' pc.cov.structure(X,Y,3,0,3)
pc.cov.structure = function(X,Y,q,j1,j2,T){
  if (is.null(Y))
    Y = X

  if (!is.matrix(X) || !is.matrix(Y))
    stop("X and Y must be matrices")
  if (dim(X)[1] != dim(Y)[1])
    stop("Number of observations must be equal")

  nbasisX = dim(X)[2]
  nbasisY = dim(Y)[2]
  n = dim(X)[1]

  Ch = array(0,c(2*q+1,nbasisX,nbasisY))

  for (h in (-q):q)
    Ch[h+q+1,,] = pc.lagged.cov(X,Y,h,j1,j2,T)
  Ch

  timedom(Ch,-q:q)
}
