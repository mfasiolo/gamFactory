.contract2 <- function(X, v, ind){
  d <- length(v)
  Xv <- matrix(nrow = nrow(X), ncol = d)
  for(ii in 1:d){
    Xv[ , ii] <- X[ , drop(ind[ii, ]), drop = FALSE] %*% v
  }
  return( Xv )
}