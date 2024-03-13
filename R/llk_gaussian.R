##################
#'
#' Log-likelhood of a Gaussian model
#' 
#' @description XXX.
#' @param param XXX.
#' @param deriv XXX.
#' @rdname llk_gaussian
#' @examples 
#' library(gamFactory)
#' n <- 10
#' y <- rnorm(n)
#' param <- c(0.5, 1.5) # mu and 1/sigma
#' 
#' # Derivatives of Gaussian log-lik up to order 3
#' llk_gaussian(y = y, param = param, deriv = 3)
#' 
#' # Wrap derivatives for compatibility with gamFactory::checkDeriv
#' obj <- list(
#'   "d0" = function(param){
#'     sum(llk_gaussian(y = y, param = param, deriv = 0)$d0)
#'   },
#'   "d1" = function(param){
#'     colSums(do.call("cbind", llk_gaussian(y = y, param = param, deriv = 1)$d1))
#'     
#'   },
#'   "d2" = function(param){
#'     colSums(do.call("cbind", llk_gaussian(y = y, param = param, deriv = 2)$d2))
#'     
#'   },
#'   "d3" = function(param){
#'     colSums(do.call("cbind", llk_gaussian(y = y, param = param, deriv = 3)$d3))
#'   })
#' 
#' check_deriv(obj = obj, param = param, ord = 1:3)
#' @export
#' 
llk_gaussian <- function(y, param, deriv = 0, ...) {
  
  if (is.list(param) ) param <- do.call("cbind", param)
  if (is.vector(param)) param <- matrix(param, nrow = 1)
  if (ncol(param) != 2) stop("Wrong number of parameters provided")
  
  p <- ncol( param )
  mu <- param[ , 1, drop = TRUE]
  tau <- param[ , 2, drop = TRUE] # 1 / sigma
  tau2 <- tau^2                   # 1 / sigma^2
  n <- length(y)
  
  if (length(mu) == 1) {
    mu <- rep(mu, n)
    tau <- rep(tau, n)
    tau2 <- rep(tau2, n)
  }
  
  ymu <- y - mu        
  ymu2 <- ymu ^ 2
    
  d0 <- - .5 * log(2 * pi) + log(tau) - .5 * tau2 * ymu2

  out <- list()
  out$d0 <- d0
  
  if( deriv > 0 )
  {
    d1 <- tau2 * ymu
    d2 <- 1/tau - tau*ymu2
    out[["d1"]] <- list(d1, d2)
    
    if( deriv > 1 ){
      d11 <- - tau2
      d12 <- 2 * d1 / tau
      d22 <- - ymu2 - 1/tau2
      out[["d2"]] <- list(d11, d12, d22) 
      
      if( deriv > 2 ){
        zeros <- rep(0, n)
        d111 <- zeros
        d112 <- - 2 * tau
        d122 <- 2 * ymu
        d222 <- 2 / tau^3
        out[["d3"]] <- list(d111, d112, d122, d222) 
        
        if( deriv > 3){
          d1111 <- d1112 <- d1222 <- zeros
          d1122 <- rep(-2, n)
          d2222 <- -6 / tau2^2
          out[["d4"]] <- list(d1111, d1112, d1122, d1222, d2222)
        }
      }
    }
  }
    
  return( out )
  
}