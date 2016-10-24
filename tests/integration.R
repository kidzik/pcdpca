library(freqdom)
library(MASS)
library(pcdpca)

library(devtools)
install(".")
library(pcdpca)

set.seed(1)
d = 2
n = 1000*d
Psi = matrix(rnorm(d^2),d,d)
Psi = Psi %*% t(Psi)
Psi = 0 * Psi / norm(Psi)
X = rar(n,d,Psi = Psi)# + sign(sin(2 * pi * matrix(rep(c(1:d), d * n), n, d) / 2))

train = 1:(n/2)
test = (n/2 + 1):n

## Dynamic PCA ##
XI.est = dprcomp(X[train,],q=10,weights="Bartlett",freq=pi*(-150:150/150))  # finds the optimal filter
Y.est = XI.est %c% X  # applies the filter
Y.est[,-1] = 0 # forces the use of only one component
Xdpca.est = t(rev(XI.est)) %c% Y.est    # deconvolution

## Periodically correlated PCA ##
XI.est = pcdpca(X[train,],q=10,weights="Bartlett",freq=pi*(-150:150/150),T=d)  # finds the optimal filter
Y.est = pcdpca.scores(X, XI.est)  # applies the filter
Y.est[,-1] = 0 # forces the use of only one component
Xpcdpca.est = pcdpca.inverse(Y.est, XI.est)  # deconvolution

cat("NMSE DPCA = ")
cat(MSE(X[test,],Xdpca.est[test,]) / MSE(X[test,],0))
cat("\nNMSE PCDPCA = ")
cat(MSE(X[test,],Xpcdpca.est[test,]) / MSE(X[test,],0))
cat("\n")
