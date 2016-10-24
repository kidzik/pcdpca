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
X = rar(n,d,Psi = Psi)

pc.cov.structure(X,X,10,1,0,2)

cov.structure(X,X,10) # should be close
# sum(timedom.norms(A)$norms)
# sum(timedom.norms(B)$norms)
# sum(timedom.norms(A - B)$norms)
