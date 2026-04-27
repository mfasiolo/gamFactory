#'
#' Derivative of si_posi Hessian w.r.t smoothing parameters
#' 
#' @name DHessDrho.si_posi
#' @rdname DHessDrho.si_posi
#' @export 
#'
DHessDrho.si_posi <- function(o, llk, DbDr){
  
  DHessDrho.nexpsm(o = o, llk = llk, DbDr = DbDr)
  
}