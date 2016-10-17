#' For a given periodically correlated multivariate process \code{X} eigendecompose it's spectral density
#' and use an inverse fourier transform to get coefficients of the optimal filters.
#'
#' @title Compute periodically correlacted DPCA filter coefficients
#' @param X multivariate stationary time series
#' @param V correlation structure between coefficients of vectors (default diagonal)
#' @param lags requested filter coefficients
#' @param q window for spectral density estimation as in \code{\link{spectral.density}}
#' @param weights as in \code{\link{spectral.density}}
#' @param freq frequency grid to estimate on as in \code{\link{spectral.density}}
#' @return principal components series
#' @references Kidzinski, Kokoszka, Jouzdani
#' Dynamic principal components of periodically correlated functional time series
#' Research report, 2016
#' @seealso \code{\link{dpca.inverse}}, \code{\link{dpca.scores}}
#' @export
pcdpca = function(X,V=NULL,lags=-10:10,q=NULL,weights=NULL,freq=NULL,T=2){
  if (T < 2)
    stop("T must be greater or equal 2")
  if (!is.matrix(X))
    stop("X must be a matrix")
  if (!is.vector(lags) || any(!is.positiveint(abs(lags))))
    stop("lags must be a vector of integers")
  if (is.null(V))
    V = diag(dim(X)[2])
  if (is.null(q))
    q = 10

  pcdpcas = list()
  for (d in 1:T){
    Ch = pc.cov.structure(X,X,q,d-1,T)
    #    Ch = cov.structure(X,X,q)
    SD = spectral.density(X,q=q,weights=weights,freq=freq,Ch=Ch)
    E = freqdom.eigen(SD)

    nbasis = dim(E$vectors)[2]

    XI = array(0,c(length(lags),nbasis,nbasis))

    for (component in 1:nbasis)
      XI[,component,] = t(exp(-(SD$freq %*% t(lags)) * 1i)) %*% E$vectors[,,component] / length(SD$freq)

    pcdpcas[[d]] = timedom(Re(XI[length(lags):1,,]),lags)
  }
  pcdpcas
}
