.contract3_vect <- function(x, v, ind){
  d <- length(v)
  xv <- numeric((d+1)*d/2)
  zz <- 1
  for(jj in 1:d){
    for(kk in jj:d){
      xv[zz] <- crossprod(x[drop(ind[jj, kk, ])], v)
      zz <- zz + 1
    }
  }
  return( xv )
}

# Matrix version (not needed at the moment)
# .contract3 <- function(X, v, ind){
#   d <- length(v)
#   Xv <- matrix(nrow = nrow(X), ncol = (d+1)*d/2)
#   zz <- 1
#   for(jj in 1:d){
#     for(kk in jj:d){
#       Xv[ , zz] <- X[ , drop(ind[jj, kk, ]), drop = FALSE] %*% v
#       zz <- zz + 1
#     }
#   }
#   return( Xv )
# }