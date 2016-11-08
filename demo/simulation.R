library("freqdom")
library(devtools)
library(fda)
install(".")
library(pcdpca)

## Prepare the process
d = 7
n = 100
A = t(t(matrix(rnorm(d*n),ncol=d,nrow=n))*7:1)
B = t(t(matrix(rnorm(d*n),ncol=d,nrow=n))*7:1)
C = t(t(matrix(rnorm(d*n),ncol=d,nrow=n))*7:1)

X = matrix(0,ncol=d,nrow=3*n)
X[3*(1:n) - 1,] = A
X[3*(1:n) - 2,] = A + B
X[3*(1:n) ,] = 2*A - B + C

basis = create.fourier.basis(nbasis=7)
X.fd = fd(t(Re(X)),basis=basis)
plot(X.fd)

## Hold out some datapoints
train = 1:(50*3)
test = (50*3) : (3*n)

## Static PCA ##
PR = prcomp(as.matrix(X[train,]))
Y1 = as.matrix(X) %*% PR$rotation
Y1[,-1] = 0
Xpca.est = Y1 %*% t(PR$rotation)

## Dynamic PCA ##
XI.est = dprcomp(as.matrix(X[train,]),q=3,lags=-2:2,weights="Bartlett",freq=pi*(-150:150/150))  # finds the optimal filter
Y.est = XI.est %c% X  # applies the filter
Y.est[,-1] = 0 # forces the use of only one component
Xdpca.est = t(rev(XI.est)) %c% Y.est    # deconvolution

## Periodically correlated PCA ##
XI.est.pc = pcdpca(as.matrix(X[train,]),q=3,lags=-2:2,weights="Bartlett",freq=pi*(-150:150/150),period=3)  # finds the optimal filter
Y.est.pc = pcdpca.scores(X, XI.est.pc)  # applies the filter
Y.est.pc[,-1] = 0 # forces the use of only one component
Xpcdpca.est = pcdpca.inverse(Y.est.pc, XI.est.pc)  # deconvolution

## Results
cat("NMSE PCA = ")
r0 = MSE(X[test,],Xpca.est[test,]) / MSE(X[test,],0)
cat(r0)
cat("\nNMSE DPCA = ")
r1 = MSE(X[test,],Xdpca.est[test,]) / MSE(X[test,],0)
cat(r1)
cat("\nNMSE PCDPCA = ")
r2 = MSE(X[test,],Xpcdpca.est[test,]) / MSE(X[test,],0)
cat(r2)
cat("\n")
