#'
#' Wrapper for PsplineDesign
#' 
#' @name constrSplineDes
#' @rdname constrSplineDes
#' @export constrSplineDes
#'
constrSplineDes <- function(k, m, lim, B){
  
  force(k); force(m); force(lim); force(B)
  
  out <- function(x, deriv){
    PsplineDesign(x = x, k = k, m = m, lim = lim, B = B, deriv = deriv)
  }
  
  return( out )
}