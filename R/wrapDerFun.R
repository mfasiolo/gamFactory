#'
#' Wrapper to fdDeriv() 
#' 
#' @description This is a wrapper useful for checking derivative of families with the fdDeriv() function. 
#' @name createGPDs
#' @param obj for instance the output of createGPD.
#' @param dropHessian if TRUE the Hessian will be transformed to a vector using flattenSymMat.
#' @rdname wrapDerFun
#' @export wrapDerFun
#'
wrapDerFun <- function(obj, dropHessian = FALSE){
  
  d0 <- function(param){
    return( obj(param, deriv = 0)$d0 )
  }
  d1 <- function(param){
    return( obj(param, deriv = 1)$d1 )
  }
  d2 <- function(param){
    out <- obj(param, deriv = 2)$d2
    if( dropHessian && is.matrix(out) ) { out <- flattenSymMat(out) }
    return( out )
  }
  d3 <- function(param){
    return( obj(param, deriv = 3)$d3 )
  }
  
  return( list("d0" = d0, "d1" = d1, "d2" = d2, "d3" = d3) )
  
}
