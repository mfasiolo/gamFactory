#'
#' Derivative of si_nexpsm Hessian w.r.t smoothing parameters
#' 
#' @name DHessDrho.si_nexpsm
#' @rdname DHessDrho.si_nexpsm
#' @export 
#'
DHessDrho.si_nexpsm <- function(o, llk, DbDr){
  
  DHessDrho.nexpsm(o = o, llk = llk, DbDr = DbDr)
  
}
