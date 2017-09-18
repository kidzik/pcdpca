library("ggplot2")
library("reshape2")
library('latex2exp')
library("freqdom")
library(pcdpca)

MSE = function(X,Y) { mean((X-Y)**2) / mean((X)**2) }
my.dpca = function (X, q=floor(ncol(X)^{1/3}), L = 30, freq = (-1000:1000/1000) * pi, Ndpc = dim(X)[2])
{
  if (!is.matrix(X))
    stop("X must be a matrix")
  res = list()
  res$spec.density = freqdom::spectral.density(X, freq = freq, q = q)
  res$filters = freqdom::dpca.filters(res$spec.density, q = L, Ndpc = Ndpc)
  res
}

## Prepare the process
period = 2
L = 6
d = 5
n = 1000
sqn = floor(n^(1/3))

RES = c()
all.times = c()

for (run in 1:100){
  times = c()
Psi = c()
for (i in 1:period){
  PsiTmp = c()
  for (j in 1:period){
    P = matrix(rnorm(d*d, sd = exp(-(1:d)/d) ),d)
    PsiTmp = cbind(PsiTmp,0.9 * P / norm(P))
  }
  Psi = rbind(Psi,PsiTmp)
}
Xd = rar(n = n,Psi = Psi )
X = t(matrix(t(Xd),d))

## Hold out some datapoints
train = 1:(nrow(X)/2)
test = (nrow(X)/2 + 1) : (nrow(X))

## Static PCA ##
start_time <- Sys.time()
PR = prcomp(as.matrix(X[train,]))
Y1 = as.matrix(X) %*% PR$rotation
Y1[,-1] = 0
Xpca.est = Y1 %*% t(PR$rotation)
end_time <- Sys.time()
times = cbind(times, end_time - start_time)

# ## Results
r0 = MSE(X[test,],Xpca.est[test,]) / MSE(X[test,],0)

## Dynamic PCA ##
start_time <- Sys.time()
XI.est = my.dpca(as.matrix(X[train,]),L = L/3, q=sqn, freq = pi*(-150:150/150), Ndpc = 1)  # finds the optimal filter
Y.est = freqdom::filter.process(X, XI.est$filters )
Xdpca.est = freqdom::filter.process(Y.est, t(rev(XI.est$filters)) )    # deconvolution
end_time <- Sys.time()
times = cbind(times, end_time - start_time)

r1 = MSE(X[test,],Xdpca.est[test,]) / MSE(X[test,],0)

## Dynamic PCA ##
start_time <- Sys.time()
XI.est = my.dpca(as.matrix(X[train,]),L = L, q=sqn, freq = pi*(-150:150/150), Ndpc = 1)  # finds the optimal filter
Y.est = freqdom::filter.process(X, XI.est$filters )
Xdpca.est = freqdom::filter.process(Y.est, t(rev(XI.est$filters)) )    # deconvolution
end_time <- Sys.time()
times = cbind(times, end_time - start_time)

r1 = c(r1, MSE(X[test,],Xdpca.est[test,]) / MSE(X[test,],0))

Ls = c(3,6)
qs = c(floor(sqn/2),sqn,2*sqn)

r2 = c()
for (L in Ls){
  for (q in qs){
    ## Periodically correlated PCA ##
    start_time <- Sys.time()
    XI.est.pc = pcdpca(as.matrix(X[train,]), q = q, freq = pi*(-150:150/150), period = period, L = L)  # finds the optimal filter
    Y.est.pc = pcdpca.scores(X, XI.est.pc)  # applies the filter
    Y.est.pc[,-1] = 0 # forces the use of only one component
    Xpcdpca.est = pcdpca.inverse(Y.est.pc, XI.est.pc)  # deconvolution
    end_time <- Sys.time()
    times = cbind(times, end_time - start_time)

    # cat("\nNMSE PCDPCA = ")
    r2 = c(r2, MSE(X[test,],Xpcdpca.est[test,]) / MSE(X[test,],0))
  }
}

# cat(r2)
# cat("\n")
row = c(r0,r1,r2)
print(row)
print(times)
RES = rbind(RES,row)
all.times = rbind(all.times, times)
}

## PLOTTING ##
nms = c(paste("PC-DPCA \n q=",qs,"; L=",Ls[1],sep=""),
        paste("PC-DPCA \n q=",qs,"; L=",Ls[2],sep=""))

nms = c("PCA","DPCA\nL=3","DPCA",nms)
colnames(RES) = nms
colnames(all.times) = nms
cols = c(c(1,3),4:9)
#cols = 1:length(nms)
df = data.frame(RES[,cols],row.names = NULL)

colMeans(RES)
apply(RES, 2, sd)

colMeans(all.times)
apply(all.times, 2, sd)

df.long = melt(df)
ggplot(df.long, aes(variable, value)) + theme_bw() + theme(text = element_text(size=30), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_x_discrete(labels=nms[cols]) + xlab("") + ylim(c(0.35,0.9)) + ylab("Normalized mean squared error") + labs(title = "Simulation study 2") +
  geom_boxplot()
