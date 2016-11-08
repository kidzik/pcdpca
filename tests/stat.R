X = t(matrix(1:120,nrow=6))

Y = stat2pc(X,2)
assertthat::are_equal(Y[2,1],2)
assertthat::are_equal(Y[3,2],9)
assertthat::are_equal(dim(Y),c(40,3))

Y = pc2stat(X,2)
assertthat::are_equal(Y[2,1],13)
assertthat::are_equal(Y[3,2],26)
assertthat::are_equal(dim(Y),c(10,12))
