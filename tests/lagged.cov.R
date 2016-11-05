library(freqdom)
library(MASS)
library(pcdpca)

library(devtools)
install(".")
library(pcdpca)

set.seed(1)
d = 2
n = 100*d
Psi = matrix(rnorm(d^2),d,d)
Psi = Psi %*% t(Psi)

X = rar(n,d,Psi = 0 * Psi / norm(Psi))


R = pc.lagged.cov(X, X, 0, 0, 1, 2)
pc.lagged.cov(X, X, -10, 0, 0, 2)

lagged.cov(X, X, 1) - t(lagged.cov(X, X, -1))

pc.lagged.cov(X, X, 1, 0, 0, 2) - t(pc.lagged.cov(X, X, -1, 0, 0, 2))
pc.lagged.cov(X, X, 1, 0, 0, 2) - t(pc.lagged.cov(X, X, -1, 1, 1, 2))
lagged.cov(X, X, 0) - t(pc.lagged.cov(X, X, 0, 1, 1, 2))

# pc.lagged.cov(X, X, 0, 0, 4) - lagged.cov(X, X, 0)


# VAR(1)
# X_(t+1) = Psi(X_t) + e_(t+1)
# e \sim N(0,Simga)
# V = Var(X_t) = \Sum_{i=0}^\infty Psi^i Sigma t(Psi^i)
# Cov(X_{t+n},X_t) = Psi^n V
Sigma = diag(2)
Sigma[1,2] = 0.3
Sigma[2,1] = Sigma[1,2]
Sigma = Sigma

noise = function(n){ mvrnorm(1,c(0,0),Sigma) }

Psi = matrix(rnorm(d^2),d,d)
Psi = Psi %*% t(Psi)
P = 0.3 * Psi / norm(Psi)

# Theoretical
V = Sigma
pow = P
for (i in 1:20){
  V = V + pow %*% Sigma %*% t(pow)
  pow = pow %*% P
}

# Empirical
X = rar(100000, d, Psi = P, noise = noise)
lagged.cov(X)

norm(t(P %*% V) - lagged.cov(X, X, 1))
norm(t(P %*% V) - pc.lagged.cov(X, X, -1, 0, 1, 2))
norm(P %*% P %*% V - pc.lagged.cov(X, X, 1, 0, 1, 3))
