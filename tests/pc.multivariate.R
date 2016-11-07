library("freqdom")
library(devtools)
install(".")
library(pcdpca)

d = 2
n = 10000
A = t(c(1,0.01)*t(matrix(rnorm(d*n),ncol=d,nrow=n)))
B = t(c(1,0.01)*t(matrix(rnorm(d*n),ncol=d,nrow=n)))
A[,2] = 0
B[,2] = 0
X = matrix(0,ncol=2,nrow=d*n)

X[2*(1:n) - 1,] = A
X[2*(1:n),] = B
X = A

cov11 = pc.cov.structure(X,X,2,1,1,2)
cov10 = pc.cov.structure(X,X,2,1,0,2)
cov00 = pc.cov.structure(X,X,2,0,0,2)
cov01 = pc.cov.structure(X,X,2,1,0,2)
cov = cov.structure(X,X,2)

train = 1:(n/2)
test = (1+n/2):n

## Dynamic PCA ##
XI.est = dprcomp(as.matrix(X[train,]),q=10,lags=-2:2,weights="Bartlett",freq=pi*(-150:150/150))  # finds the optimal filter
Y.est = XI.est %c% X  # applies the filter
Y.est[,-1] = 0 # forces the use of only one component
Xdpca.est = t(rev(XI.est)) %c% Y.est    # deconvolution
sum(Y.est^2)

## Periodically correlated PCA ##
XI.est.pc = pcdpca(as.matrix(X[train,]),q=1,lags=-2:2,weights="Bartlett",freq=pi*(-50:50/50),T=d)  # finds the optimal filter
Y.est.pc = pcdpca.scores(X, XI.est.pc)  # applies the filter
Y.est.pc[,-1] = 0 # forces the use of only one component
Xpcdpca.est = pcdpca.inverse(Y.est.pc, XI.est.pc)  # deconvolution
sum(Y.est.pc^2)

cat("NMSE DPCA = ")
r1 = MSE(X[test,],Xdpca.est[test,]) / MSE(X[test,],0)
cat(r1)
cat("\nNMSE PCDPCA = ")
r2 = MSE(X[test,],Xpcdpca.est[test,]) / MSE(X[test,],0)
cat(r2)
cat("\n")
RES = rbind(RES, c(r1,r2))
