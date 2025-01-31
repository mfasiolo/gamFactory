# P-splines penalty matrix
.Ppen <- function(o, k){
  if(o == 0){ return(diag(1, k)) }
  P <- diff(diag(k), differences = o)
  return( t(P) %*% P )
}