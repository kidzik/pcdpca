library(freqdom)
library(MASS)
library(fda)
library(pcdpca)

data(pm10)

n = dim(X$coef)[2]

## Static PCA ##
PR = prcomp(t(X$coef))
Y1 = PR$x
Y1[,-1] = 0
Xpca.est = Y1 %*% t(PR$rotation)

## Dynamic PCA ##
XI.est = dprcomp(t(X$coef),q=8,lags=-10:10,weights="Bartlett",freq=pi*(-150:150/150))  # finds the optimal filter
Y.est = XI.est %c% t(X$coef)  # applies the filter
Y.est[,-1] = 0 # forces the use of only one component
Xdpca.est = t(rev(XI.est)) %c% Y.est    # deconvolution

## Periodically correlated PCA ##
XI.est = pcdpca(t(X$coef),q=2,lags=-2:2,weights="Bartlett",freq=pi*(-150:150/150),period=7)  # finds the optimal filter
Y.est = pcdpca.scores(t(X$coef), XI.est)  # applies the filter
Y.est[,-1] = 0 # forces the use of only one component
Xpcdpca.est = pcdpca.inverse(Y.est, XI.est)  # deconvolution

# Creates functional objects
Xpcdpca.est.fd = fd(t(Re(Xpcdpca.est)),basis=X$basis)
Xdpca.est.fd = fd(t(Re(Xdpca.est)),basis=X$basis)
Xpca.fd = fd(t(Xpca.est),basis=X$basis)

# Write down results
ind = 5:(n-5)
cat("NMSE PCA =  ")
cat(MSE(t(X$coef)[ind,],Xpca.est[ind,]) / MSE(t(X$coef)[ind,],0))
cat("\nNMSE DPCA = ")
cat(MSE(t(X$coef)[ind,],Xdpca.est[ind,]) / MSE(t(X$coef)[ind,],0))
cat("\nNMSE PCDPCA = ")
cat(MSE(t(X$coef)[ind,],Xpcdpca.est[ind,]) / MSE(t(X$coef)[ind,],0))
cat("\n")
