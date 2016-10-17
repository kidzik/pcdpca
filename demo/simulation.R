library(freqdom)
library(MASS)
library(pcdpca)

set.seed(1)
n = 2001
d = 5
Psi = matrix(rnorm(d^2),d,d)
Psi = Psi %*% t(Psi)
Psi = 0.2 * Psi / norm(Psi)
X = rar(n,d,Psi = Psi) + matrix(rep(c(1:d),d*n),n,d)

## Static PCA ##
PR = prcomp(X)
Y1 = PR$x
Y1[,-1] = 0
Xpca.est = Y1 %*% t(PR$rotation)

## Dynamic PCA ##
XI.est = dprcomp(X,q=4,weights="Bartlett",freq=pi*(-150:150/150))  # finds the optimal filter
Y.est = XI.est %c% X  # applies the filter
Y.est[,-1] = 0 # forces the use of only one component
Xdpca.est = t(rev(XI.est)) %c% Y.est    # deconvolution

## Periodically correlated PCA ##
XI.est = pcdpca(X,q=3,weights="Bartlett",freq=pi*(-150:150/150),T=2)  # finds the optimal filter
Y.est = pcdpca.scores(X, XI.est)  # applies the filter
Y.est[,-1] = 0 # forces the use of only one component
Xpcdpca.est = pcdpca.inverse(Y.est, XI.est)  # deconvolution

# Write down results
ind = 1:n
cat("NMSE PCA =  ")
cat(MSE(X[ind,],Xpca.est[ind,]) / MSE(X[ind,],0))
cat("\nNMSE DPCA = ")
cat(MSE(X[ind,],Xdpca.est[ind,]) / MSE(X[ind,],0))
cat("\nNMSE PCDPCA = ")
cat(MSE(X[ind,],Xpcdpca.est[ind,]) / MSE(X[ind,],0))
cat("\n")
