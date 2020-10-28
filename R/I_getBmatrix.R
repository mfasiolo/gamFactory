#################
# Reparametrisation matrix to work with diagonal penalty
# 
.getBmatrix <- function(P, r){
  
  d <- ncol( P )
  ei <- eigen( P )
  v <- ei$values
  
  if(r < d){
    v[ (r+1):d ] <- 1
  }
  
  B <- ei$vectors * 1 / sqrt(v)
  
  return( B )
}
