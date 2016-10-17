inClass = function(X,cls){
  !is.null(oldClass(X)) && oldClass(X) == cls
}

is.positiveint = function (n){
  is.numeric(n) && n > 0
}
