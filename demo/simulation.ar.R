library("freqdom")
library(pcdpca)

## Prepare the process
period = 2
d = 5
n = 2000

RES = c()

for (run in 1:100){
Psi = c()
for (i in 1:period){
  PsiTmp = c()
  for (j in 1:period){
    P = matrix(rnorm(d*d, sd = exp(-(1:d)/d) ),d)
    P = P %*% t(P)
    PsiTmp = cbind(PsiTmp,0.5 * P / norm(P)) * sign(runif(1,-1,1))
  }
  Psi = rbind(Psi,PsiTmp)
}
Xd = rar(n = n,Psi = Psi,noise = function(n) { rnorm(n, exp( - (1:n) / n))  } )
X = t(matrix(t(Xd),d))

## Hold out some datapoints
train = 1:(nrow(X)/2)
test = (nrow(X)/2) : (nrow(X))

## Static PCA ##
PR = prcomp(as.matrix(X[train,]))
Y1 = as.matrix(X) %*% PR$rotation
Y1[,-1] = 0
Xpca.est = Y1 %*% t(PR$rotation)

sqn = floor(sqrt(n))
lags = -sqn:sqn

## Dynamic PCA ##
XI.est = dprcomp(as.matrix(X[train,]),q=sqn,lags=-lags,weights="Bartlett",freq=pi*(-150:150/150))  # finds the optimal filter
Y.est = XI.est %c% X  # applies the filter
Y.est[,-1] = 0 # forces the use of only one component
Xdpca.est = t(rev(XI.est)) %c% Y.est    # deconvolution

## Periodically correlated PCA ##
XI.est.pc = pcdpca(as.matrix(X[train,]),q=sqn,lags=lags,weights="Bartlett",freq=pi*(-150:150/150),period=period)  # finds the optimal filter
Y.est.pc = pcdpca.scores(X, XI.est.pc)  # applies the filter
Y.est.pc[,-1] = 0 # forces the use of only one component
Xpcdpca.est = pcdpca.inverse(Y.est.pc, XI.est.pc)  # deconvolution

## Results
r0 = MSE(X[test,],Xpca.est[test,]) / MSE(X[test,],0)
r1 = MSE(X[test,],Xdpca.est[test,]) / MSE(X[test,],0)
r2 = MSE(X[test,],Xpcdpca.est[test,]) / MSE(X[test,],0)
row = c(r0,r1,r2)
print(row)
RES = rbind(RES,row)
}

colnames(RES) = c("PCA","DPCA","PCDPCA")
df = data.frame(RES,row.names = NULL)

colMeans(df)
summary(df)
apply(df, 2, sd)

t.test(df$DPCA - df$PCDPCA)
boxplot(df)