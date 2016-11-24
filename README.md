# pcdpca

Implementation of "Dynamic principal components of periodically correlated functional time series".

Two examples in `demo` directory:

  - pm10 data from Graz (comparizon with DFPCA paper)
  - simplation with parametrized periodicity

## Installation

    library("devtools")
    install_github("kidzik/pcdpca"")

## Running a demo

    library("pcdpca")
    demo("simulation")
    demo("pcdpca.pm10")

## Usage

Let `X` be a multivariate time series, a matrix with `n` observations and `d` covariates.
Let `period` be the period. Then

    XI.est = pcdpca(X,q=3,weights="Bartlett",freq=pi*(-150:150/150),period=2)  # finds the optimal filter
    Y.est = pcdpca.scores(X, XI.est)  # applies the filter
    Y.est[,-1] = 0 # forces the use of only one component
    Xpcdpca.est = pcdpca.inverse(Y.est, XI.est)  # deconvolution
    cat(MSE(X,Xpcdpca.est[ind,]) / MSE(X[ind,],0)) # variance explained
