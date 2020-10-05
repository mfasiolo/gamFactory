#'
#' Wrapper for PsplineDesign
#' 
#' @name constrSplineDes
#' @rdname constrSplineDes
#' @export constrSplineDes
#'
constrSplineDes <- function(k, m, lim, B, NS){
  
  force(k); force(m); force(lim); force(B); force(NS)
  
  out <- function(x, deriv){
    PsplineDesign(x = x, k = k, m = m, lim = lim, B = B, NS = NS, deriv = deriv)
  }
  
  return( out )
}