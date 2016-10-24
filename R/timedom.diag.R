#' @export
timedom.diag = function(A,T)
{
  res = list()
  d = dim(A$operators)[3]
  dr = floor(d / T)


  for (i in 0:(T-1)){
    C = A
    idx = 1:dr + dr * i
    C$operators = A$operators[,idx,idx]
    res[[i+1]] = C
  }
  res
}

