library(freqdom)
library(MASS)
library(pcdpca)

library(devtools)
install(".")
library(pcdpca)

set.seed(1)
d = 2
n = 100*d
Psi = matrix(rnorm(d^2),d,d)
Psi = Psi %*% t(Psi)
Psi = 0 * Psi / norm(Psi)
X = rar(n,d,Psi = Psi)


R = pc.lagged.cov(X, X, 0, 0, 1, 2)
pc.lagged.cov(X, X, -10, 0, 0, 2)

lagged.cov(X, X, 1) - t(lagged.cov(X, X, -1))

pc.lagged.cov(X, X, 1, 0, 0, 2) - t(pc.lagged.cov(X, X, -1, 0, 0, 2))

# pc.lagged.cov(X, X, 0, 0, 4) - lagged.cov(X, X, 0)

