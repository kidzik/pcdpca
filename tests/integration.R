library(freqdom)
library(MASS)
library(pcdpca)

library(devtools)
install(".")
library(pcdpca)

set.seed(1)
d = 2
n = 1000*d

Sigma = diag(2)
Sigma[1,2] = 0.3
Sigma[2,1] = Sigma[1,2]
Sigma = Sigma

noise = function(n){ mvrnorm(1,c(0,0),Sigma) }

Psi = matrix(rnorm(d^2),d,d)
Psi = Psi %*% t(Psi)
P = 0.3 * Psi / norm(Psi)

X = rar(n, d, Psi = P, noise = noise)

train = 1:(n/2)
test = (n/2 + 1):n

## Dynamic PCA ##
XI.est = dprcomp(X[train,],q=3,lags=-2:2,weights="Bartlett",freq=pi*(-150:150/150))  # finds the optimal filter
Y.est = XI.est %c% X  # applies the filter
Y.est[,-1] = 0 # forces the use of only one component
Xdpca.est = t(rev(XI.est)) %c% Y.est    # deconvolution
sum(Y.est^2)

## Periodically correlated PCA ##
XI.est.pc = pcdpca(X[train,],q=3,lags=-2:2,weights="Bartlett",freq=pi*(-150:150/150),T=d)  # finds the optimal filter
Y.est.pc = pcdpca.scores(X, XI.est.pc)  # applies the filter
Y.est.pc[,-1] = 0 # forces the use of only one component
Xpcdpca.est = pcdpca.inverse(Y.est.pc, XI.est.pc)  # deconvolution
sum(Y.est.pc^2)

cat("NMSE DPCA = ")
cat(MSE(X[test,],Xdpca.est[test,]) / MSE(X[test,],0))
cat("\nNMSE PCDPCA = ")
cat(MSE(X[test,],Xpcdpca.est[test,]) / MSE(X[test,],0))
cat("\n")
