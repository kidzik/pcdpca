library(freqdom)
library(MASS)
library(pcdpca)

library(devtools)
install(".")
library(pcdpca)

set.seed(1)
d = 2
n = 10000*d

Psi = matrix(rnorm(d^2),d,d)
Psi = Psi %*% t(Psi)
X = rar(n,d,Psi = 0 * Psi / norm(Psi))

pc.cov.structure(X,X,2,1,0,2)

A = cov.structure(X,X,2)
B = pc.cov.structure(X,X,2,0,0,2) # should be close
C = pc.cov.structure(X,X,2,1,0,2)
D = pc.cov.structure(X,X,2,1,1,2)
sum(timedom.norms(A - B)$norms)
sum(timedom.norms(A - D)$norms)

X = rar(n,d,Psi = 0.5 * Psi / norm(Psi))
A = cov.structure(X,X,2)
C = pc.cov.structure(X,X,2,0,0,2)
D = pc.cov.structure(X,X,2,1,1,2)
sum(timedom.norms(C - D)$norms)

C = pc.cov.structure(X,X,1,1,0,2)
D = pc.cov.structure(X,X,1,0,1,2)
sum(timedom.norms(C - t(D) )$norms)


# sum(timedom.norms(A)$norms)
# sum(timedom.norms(B)$norms)
# sum(timedom.norms(A - B)$norms)
