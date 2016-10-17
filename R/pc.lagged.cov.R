#' For a given multivariate stationary periodically correlated time series estimates a covarianve matrix
#' \eqn{C_{XY}^k_d = Cov(X_{d+k},Y_d)} using the formula
#' \deqn{\hat C_{XY}^k_d =  \frac{1}{[n/T]} \sum_{i=1}^{[n/T]} X_{k+i*T+d} Y_{k+d}'. }
#'
#' @title Compute cross covariance with a given lag for periodically correlated time series
#' @param X first process
#' @param Y second process, if null then autocovariance of X is computed
#' @param lag the lag that we are interested in
#' @param d index in period
#' @param T period
#' @return Covariance matrix
#' @export
#' @examples
#' X = rar(100)
#' Y = rar(100)
#' pc.lagged.cov(X,Y)
pc.lagged.cov = function(X,Y,lag,d,T)
{
  if (is.null(Y))
		Y = X

  if (dim(X)[1] != dim(Y)[1])
    stop("Number of observations must be equal")
  if (!is.matrix(X) || !is.matrix(Y))
    stop("X and Y must be matrices")

	n = dim(X)[1]
	h = abs(lag)

  if (n - 1 <= h)
	  stop(paste("Too little observations to compute lagged covariance with lag",h))

	idxX = 1:(n-h)
	idxY = 1:(n-h)+h

	subidx = T*(0:n) + d + 1
	subidx = subidx[subidx <= n-h]
	subidx = subidx[subidx >= 0]
	idxX = idxX[subidx]
	idxY = idxY[subidx]

  M = t(X[idxX,]) %*% (Y[idxY,])/(n)

	if (lag < 0){
	  M = t(Y[idxX,]) %*% (X[idxY,])/(n)
	  M = t(M)
	}
	M
}

