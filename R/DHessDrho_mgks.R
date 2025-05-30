#'
#' Derivative of MGKS Hessian w.r.t smoothing parameters
#' 
#' @name DHessDrho.mgks
#' @rdname DHessDrho.mgks
#' @export 
#'
DHessDrho.mgks <- function(o, llk, DbDr){
  
  DHessDrho.nexpsm(o = o, llk = llk, DbDr = DbDr)
  
}
  
