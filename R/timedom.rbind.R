#' @export
timedom.rbind = function(A,B)
{
  C = A
  nlags = length(A$lags)
  dims = c(nlags, dim(A$operators)[2] + dim(B$operators)[2], dim(A$operators)[3])
  C$operators = array(0, dims)
  for (i in 1:nlags){
    C$operators[i,,] = rbind(A$operators[i,,],B$operators[i,,])
  }
  C
}

