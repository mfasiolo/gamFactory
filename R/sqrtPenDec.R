#'
#' Penalise null space of penalty
#' 
#' @name sqrtPenDec
#' @rdname sqrtPenDec
#' @export sqrtPenDec
#'
sqrtPenDec <- function(P, r, eps = 1e-5){
  
  d <- ncol( P )
  ei <- eigen( P )
  v <- ei$values
  
  if( r < d ){ 
    v[ -(1:r) ] <- max(v) * eps
  }
  
  B <- t( t(ei$vectors) * sqrt(v) )
  
  return( B )
}
