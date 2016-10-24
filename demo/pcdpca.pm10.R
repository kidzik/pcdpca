library(freqdom)
library(MASS)
library(fda)
data(pm10)

library(devtools)
# install(".")
library(pcdpca)

## Static PCA ##
PR = prcomp(t(X$coef))
Y1 = PR$x
Y1[,-1] = 0
Xpca.est = Y1 %*% t(PR$rotation)

## Dynamic PCA ##
XI.est = dprcomp(t(X$coef),q=4,weights="Bartlett",freq=pi*(-150:150/150))  # finds the optimal filter
Y.est = XI.est %c% t(X$coef)  # applies the filter
Y.est[,-1] = 0 # forces the use of only one component
Xdpca.est = t(rev(XI.est)) %c% Y.est    # deconvolution

## Periodically correlated PCA ##
XI.est = pcdpca(t(X$coef),q=3,weights="Bartlett",freq=pi*(-150:150/150),T=7)  # finds the optimal filter
Y.est = pcdpca.scores(t(X$coef), XI.est)  # applies the filter
Y.est[,-1] = 0 # forces the use of only one component
Xpcdpca.est = pcdpca.inverse(Y.est, XI.est)  # deconvolution

# Creates functional objects
Xpcdpca.est.fd = fd(t(Re(Xpcdpca.est)),basis=X$basis)
Xdpca.est.fd = fd(t(Re(Xdpca.est)),basis=X$basis)
Xpca.fd = fd(t(Xpca.est),basis=X$basis)

# Write down results
n = dim(X$coef)[2]
ind = 1:n
cat("NMSE PCA =  ")
cat(MSE(t(X$coef)[ind,],Xpca.est[ind,]) / MSE(t(X$coef)[ind,],0))
cat("\nNMSE DPCA = ")
cat(MSE(t(X$coef)[ind,],Xdpca.est[ind,]) / MSE(t(X$coef)[ind,],0))
cat("\nNMSE PCDPCA = ")
cat(MSE(t(X$coef)[ind,],Xpcdpca.est[ind,]) / MSE(t(X$coef)[ind,],0))
cat("\n")
