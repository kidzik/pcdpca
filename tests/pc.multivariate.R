library("freqdom")
library(devtools)
install(".")
library(pcdpca)

n = 100
A = rar(n,2,matrix(0.3*c(1,0.3,0.3,1),nrow = 2)) + c(-1,2)
B = rar(n,2,matrix(0.3*c(1,0.9,0.1,1),nrow = 2)) + c(3,1)
C = matrix(0,ncol=2,nrow=2*n)

C[2*(1:n),] = A
C[2*(1:n) - 1,] = B + 3*A
