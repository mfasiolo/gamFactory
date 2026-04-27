#################
# Reparametrisation matrix to work with diagonal penalty
#
# r is rank of penalty matrix P
.getBmatrix <- function(P, r){

  d <- ncol( P )
  
  # if P is very close to identity, return identity matrix
  # It would work even without this special cause but this getting a B 
  # that is a permutation matrix. 
  if( max(abs(P - diag(d))) < 1e-6 ){ 
    return( diag(d) )
  }
  
  ei <- eigen( P )
  v <- ei$values

  if(r < d){
    v[ (r+1):d ] <- 1
  }

  B <- ei$vectors %*% diag(1/sqrt(v), d)

  return( B )
}
