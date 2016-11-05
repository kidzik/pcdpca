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
pc.lagged.cov = function(X,Y,lag,j1,j2,T)
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

	# idxX = 1:(n-h)
	# idxY = 1:(n-h)+h

	idY = T*(0:n) + j1 + 1 + h * T
	idY = idY[idY <= n]
	# print(idY)
	nr = length(idY)

	idX = T*(0:n) + j2 + 1
	idX = idX[1:nr]
	# print(idX)


  M = t(X[idX,]) %*% (Y[idY,])/(nr)

	if (lag < 0){
	  M = t(Y[idX,]) %*% (X[idY,])/(nr)
	  M = t(M)
	}
	t(M)
}

