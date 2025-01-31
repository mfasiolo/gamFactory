
##################
#'
#' Log-likelhood of a Gaussian model
#' 
#' @description XXX.
#' @param param XXX.
#' @param deriv XXX.
#' @rdname llk_binomial
#' @export llk_binomial
#' @examples 
#' library(gamFactory)
#' nobs <- 10
#' n <- 50
#' y <- rbinom(nobs, size = n, prob = 0.3)
#' param <- c(0.4)
#' 
#' # Derivatives of Gaussian log-lik up to order 3
#' llk_binomial(y = y, param = param, n = n, deriv = 3)
#' 
#' # Wrap derivatives for compatibility with gamFactory::checkDeriv
#' obj <- list(
#'   "d0" = function(param){
#'     sum(llk_binomial(y = y, param = param, n = n, deriv = 0)$d0)
#'   },
#'   "d1" = function(param){
#'     colSums(do.call("cbind", llk_binomial(y = y, param = param, n = n, deriv = 1)$d1))
#'     
#'   },
#'   "d2" = function(param){
#'     colSums(do.call("cbind", llk_binomial(y = y, param = param, n = n, deriv = 2)$d2))
#'     
#'   },
#'   "d3" = function(param){
#'     colSums(do.call("cbind", llk_binomial(y = y, param = param, n = n, deriv = 3)$d3))
#'   })
#' 
#' check_deriv(obj = obj, param = param, ord = 1:3)
#' 
llk_binomial <- function(y, param, n, deriv = 0, ...) {
  
  if (is.list(param) ) param <- do.call("cbind", param)
  if (is.vector(param)) param <- matrix(param, ncol = 1)
  if (ncol(param) != 1) stop("Wrong number of parameters provided")
  
  
  p <- ncol( param )
  mu <- param[ , 1, drop = TRUE] # rate
  nobs <- length(y)
  
  if(length(n) != nobs){
    if(length(n) == 1){
      n <- rep(n, nobs)
    } else {
      stop("length(n) != length(y)")
    }
  }
  
  if (length(mu) == 1) {
    mu <- rep(mu, nobs)
  }
  
  lmu <- log(mu)
  nmy <- n - y
  omm <- 1 - mu
  
  d0 <- lchoose(n = n, k = y) + y * lmu + nmy * log1p(-mu) 
  
  out <- list()
  out$d0 <- d0
  
  if( deriv > 0 )
  {
    d1 <- y / mu - nmy / omm
    out[["d1"]] <- list(d1)
    
    if( deriv > 1 ){
      d11 <- - y / mu^2 - nmy / omm^2
      out[["d2"]] <- list(d11) 
      
      if( deriv > 2 ){
        d111 <- 2 * (y / mu^3 - nmy / omm^3)
        out[["d3"]] <- list(d111) 
        
        if( deriv > 3){
          d1111 <- 6 * (- y / mu^4 - nmy / omm^4)
          out[["d4"]] <- list(d1111)
        }
      }
    }
  }
  
  return( out )
  
}