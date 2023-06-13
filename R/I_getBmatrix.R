#################
# Reparametrisation matrix to work with diagonal penalty
#
# r is rank of penalty matrix P
.getBmatrix <- function(P, r){

  d <- ncol( P )
  ei <- eigen( P )
  v <- ei$values

  if(r < d){
    v[ (r+1):d ] <- 1
  }

  B <- ei$vectors %*% diag(1/sqrt(v), d)

  return( B )
}
