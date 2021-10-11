#'
#' Derivative of standard Hessian w.r.t smoothing parameters
#' 
#' @name DHessDrho.stand
#' @rdname DHessDrho.stand
#' @export 
#'
DHessDrho.stand <- function(o, llk, DbDr){
  
  X <- o$store$X
  m <- ncol( DbDr )        # Number of smoothing parameters
  DeDr <- X %*% DbDr       # Deta/Drho
  V <- llk$d3 * DeDr

  d1H <- list()
  for(ii in 1:m) {
    d1H[[ii]] <- crossprod(X, V[ , ii] * X)
  } 

  return(d1H)
}







