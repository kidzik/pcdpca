if (!requireNamespace("fda", quietly = TRUE)) {
  stop("fda package is needed for this demo to work. Please install it.",
       call. = FALSE)
}
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("MASS package is needed for this demo to work. Please install it.",
       call. = FALSE)
}

library(MASS)
library(fda)
library(freqdom)
library(pcdpca)
data(pm10)

n = dim(X$coef)[2]

rev.freqdom = function(XI){
  XI$freq = rev(XI$freq)
  XI
}

## Static PCA ##
PR = prcomp(t(X$coef))
Y1 = PR$x
Y1[,-1] = 0
Xpca = Y1 %*% t(PR$rotation)

## Dynamic PCA ##
XI.est = dprcomp(t(X$coef),q=20,weights="Bartlett",freq=pi*(-150:150/150))  # finds the optimal filter
Y.est = XI.est %c% t(X$coef)  # applies the filter
Y.est[,-1] = 0 # forces the use of only one component
Xdpca.est = t(rev(XI.est)) %c% Y.est    # deconvolution

## Periodically correlated PCA ##
XI.est = pcdpca(t(X$coef),q=2,lags=-2:2,weights="Bartlett",freq=pi*(-150:150/150),period=7)  # finds the optimal filter
Y.pcest = pcdpca.scores(t(X$coef), XI.est)  # applies the filter
Y.pcest[,-1] = 0 # forces the use of only one component
Xpcdpca.est = pcdpca.inverse(Y.pcest, XI.est)  # deconvolution

# Creates functional objects
Xpcdpca.est.fd = fd(t(Re(Xpcdpca.est)),basis=X$basis)
Xdpca.est.fd = fd(t(Re(Xdpca.est)),basis=X$basis)
Xpca.fd = fd(t(Xpca),basis=X$basis)

# Write down results
ind = 1:n
cat("NMSE DPCA =  ")
cat(MSE(t(X$coef)[ind,],Xdpca.est[ind,]) / MSE(t(X$coef)[ind,],0))
cat("\nNMSE PCA = ")
cat(MSE(t(X$coef)[ind,],Xpca[ind,]) / MSE(t(X$coef)[ind,],0))
cat("\nNMSE PCDPCA = ")
cat(MSE(t(X$coef)[ind,],Xpcdpca.est[ind,]) / MSE(t(X$coef)[ind,],0))
cat("\n")

# Figure: 10 observations reconstructed from the first component
ind = 1:10 + 20
par(mfrow=c(1,4),ps = 12, cex = 1, cex.main = 1)
plot(X[ind],ylim=c(-5,3),xlab="Intraday time", ylab="Sqrt(PM10)",lwd=2)
title("Original curves")
plot(Xpca.fd[ind],ylim=c(-5,3),xlab="Intraday time", ylab="Sqrt(PM10)",lwd=2)
title("PCA curves")
plot(Xdpca.est.fd[ind],ylim=c(-5,3),xlab="Intraday time", ylab="Sqrt(PM10)",lwd=2)
title("DPCA curves")
plot(Xpcdpca.est.fd[ind],ylim=c(-5,3),xlab="Intraday time", ylab="Sqrt(PM10)",lwd=2)
title("PCDPCA curves")
par(mfrow=c(1,1))
