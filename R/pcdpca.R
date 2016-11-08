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
#' @param period period of the periodic time series
#' @return principal components series
#' @references Kidzinski, Kokoszka, Jouzdani
#' Dynamic principal components of periodically correlated functional time series
#' Research report, 2016
#' @seealso \code{\link{pcdpca.inverse}}, \code{\link{pcdpca.scores}}
#' @export
pcdpca = function(X,V=NULL,lags=-10:10,q=NULL,weights=NULL,freq=NULL,period=NULL){
  if (is.null(period))
    stop("You have to specify the period.")
  if (period <= 1)
    stop("Period must be greater or equal 2. Otherwise use the freqdom package.")
  if (!is.null(V)){
    V = kronecker(diag(period),V)
  }
  XI = dprcomp(pc2stat(X,period=period),V=V,lags=lags,q=q,weights=weights,freq=freq)
  XI$period = period
  XI
}
