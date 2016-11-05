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
  if (!is.matrix(X))
    stop("X must be a matrix")
  if (!is.vector(lags) || any(!is.positiveint(abs(lags))))
    stop("lags must be a vector of integers")
  if (is.null(V))
    V = diag(dim(X)[2])
  if (is.null(q))
    q = 10

  pcdpcas = list()
  CVStructure = c()
  for (j1 in 0:(T-1)){
    Row = c()
    for (j2 in 0:(T-1)){
      Ch = pc.cov.structure(X,X,q,j1,j2,T)
      if (j2 == 0)
        Row = Ch
      else
        Row = timedom.cbind(Row,Ch)
    }
    if (j1 == 0)
      CVStructure = Row
    else
      CVStructure = timedom.rbind(CVStructure,Row)
  }

  SD = spectral.density(X,q=q,weights=weights,freq=freq,Ch=CVStructure)
  E = freqdom.eigen(SD)

  nbasis = dim(E$vectors)[2]

  XI = array(0,c(length(lags),nbasis,nbasis))

  for (component in 1:nbasis)
    XI[,component,] = t(exp(-(SD$freq %*% t(lags)) * 1i)) %*% E$vectors[,,component] / length(SD$freq)

  PC = timedom(Re(XI[length(lags):1,,]),lags)

  s = dim(PC$operators)[2]
  n = nrow(X)
  d = ncol(X)
  T = s / d

  XI = list()
  lags = PC$lags
  nlags = length(lags)

  XI[[1]] = array(0,c(nlags,T,T))
  XI[[2]] = array(0,c(nlags,T,T))

  ll0 = seq(1,nlags,by=2)
  ll1 = seq(2,nlags-1,by=2)

  #  print(dim(XI[[1]][ll0,,]))
  #  print(dim(PC$operators[seq(6,16,by=1), 1:2, 1:2]))

  shift = (nlags-1)/2
  midlag = 1 + shift

  # print(ll0)
  # print(seq(midlag - shift/2,midlag + shift/2,by=1))
  # print(ll1)
  # print(seq( 1+ midlag - shift/2,midlag + shift/2,by=1))

  XI[[1]][ll0,,] = PC$operators[seq(midlag - shift/2,midlag + shift/2,by=1), 1:2, 1:2]
  XI[[1]][ll1,,] = PC$operators[seq(1 + midlag - shift/2,midlag + shift/2,by=1), 3:4, 1:2]

  XI[[2]][ll1,,] = PC$operators[seq(midlag - shift/2,midlag + shift/2 - 1,by=1), 1:2, 3:4]
  XI[[2]][ll0,,] = PC$operators[seq(midlag - shift/2,midlag + shift/2,by=1), 3:4, 3:4]

  XI[[1]] = timedom(Re(XI[[1]][1:length(lags),,]),PC$lags)
  XI[[2]] = timedom(Re(XI[[2]][1:length(lags),,]),PC$lags)
  XI
}
