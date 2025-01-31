#'
#' Wrapper to fdDeriv() 
#' 
#' @description This is a wrapper useful for checking derivative of families with the fdDeriv() function. 
#' @name createGPDs
#' @param obj for instance the output of createGPD.
#' @param dropH if TRUE the Hessian will be transformed to a vector using flattenSymMat.
#' @rdname wrap_der_fun
#' @export wrap_der_fun
#'
wrap_der_fun <- function(obj, dropH = FALSE){
  
  d0 <- function(param){
    return( obj(param, deriv = 0)$d0 )
  }
  d1 <- function(param){
    return( obj(param, deriv = 1)$d1 )
  }
  d2 <- function(param){
    out <- obj(param, deriv = 2)$d2
    if( dropH && is.matrix(out) ) { out <- symmat2vec(out) }
    return( out )
  }
  d3 <- function(param){
    return( obj(param, deriv = 3)$d3 )
  }
  
  return( list("d0" = d0, "d1" = d1, "d2" = d2, "d3" = d3) )
  
}
