
##################
#'
#' Log-likelhood of a Gaussian model
#' 
#' @description XXX.
#' @param param XXX.
#' @param deriv XXX.
#' @rdname llk_poisson
#' @export llk_poisson
#' @examples 
#' library(gamFactory)
#' n <- 10
#' y <- rpois(n, 15)
#' param <- c(12.5) # mu and 1/sigma
#' 
#' # Derivatives of Gaussian log-lik up to order 3
#' llk_poisson(y = y, param = param, deriv = 3)
#' 
#' # Wrap derivatives for compatibility with gamFactory::checkDeriv
#' obj <- list(
#'   "d0" = function(param){
#'     sum(llk_poisson(y = y, param = param, deriv = 0)$d0)
#'   },
#'   "d1" = function(param){
#'     colSums(do.call("cbind", llk_poisson(y = y, param = param, deriv = 1)$d1))
#'     
#'   },
#'   "d2" = function(param){
#'     colSums(do.call("cbind", llk_poisson(y = y, param = param, deriv = 2)$d2))
#'     
#'   },
#'   "d3" = function(param){
#'     colSums(do.call("cbind", llk_poisson(y = y, param = param, deriv = 3)$d3))
#'   })
#' 
#' check_deriv(obj = obj, param = param, ord = 1:3)
#'
#' 
llk_poisson <- function(y, param, deriv = 0, ...) {
  
  if (is.list(param) ) param <- do.call("cbind", param)
  if (is.vector(param)) param <- matrix(param, ncol = 1)
  if (ncol(param) != 1) stop("Wrong number of parameters provided")
  
  p <- ncol( param )
  mu <- param[ , 1, drop = TRUE] # rate
  n <- length(y)
  
  if (length(mu) == 1) {
    mu <- rep(mu, n)
  }
  
  lmu <- log(mu)
  
  d0 <- y * lmu - mu - lfactorial(y)
  
  out <- list()
  out$d0 <- d0
  
  if( deriv > 0 )
  {
    d1 <- y / mu - 1
    out[["d1"]] <- list(d1)
    
    if( deriv > 1 ){
      d11 <- - y / mu^2
      out[["d2"]] <- list(d11) 
      
      if( deriv > 2 ){
        d111 <- 2 * y / mu^3
        out[["d3"]] <- list(d111) 
        
        if( deriv > 3){
          d1111 <- - 6 * y / mu^4
          out[["d4"]] <- list(d1111)
        }
      }
    }
  }
  
  return( out )
  
}