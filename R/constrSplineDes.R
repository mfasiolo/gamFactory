#'
#' Wrapper for PsplineDesign
#' 
#' @name constrSplineDes
#' @rdname constrSplineDes
#' @export constrSplineDes
#'
constrSplineDes <- function(k, m, lim){
  
  force(k); force(m); force(lim)

  out <- function(x, deriv){
    PsplineDesign(x = x, k = k, m = m, lim = lim, deriv = deriv)
  }
  
  return( out )
}