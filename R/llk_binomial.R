
##################
#'
#' Log-likelihood of a Binomial model
#'
#' @description Log-likelihood of the Binomial distribution, and its derivatives with
#'              respect to \code{param = mu}, the success probability (\code{y} counts
#'              successes out of \code{n} trials, both on the response, natural scale).
#'              Returned in the list format used by [`gamFactory::llk_gaussian`]
#'              and friends.
#' @param y a vector of observed counts (number of successes out of \code{n} trials).
#' @param param a matrix (or list) with 1 column (element), the success probability
#'              \code{mu}.
#' @param n number of trials, either a single value (recycled for every observation) or
#'          a vector of the same length as \code{y}.
#' @param deriv integer between 0 and 4 indicating the maximum derivative order to
#'              return: 0 only returns \code{d0} (the log-density itself), while 1-4
#'              additionally return \code{d1}-\code{d4}.
#' @rdname llk_binomial
#' @export llk_binomial
#' @examples 
#' library(gamFactory)
#' nobs <- 10
#' n <- 50
#' y <- rbinom(nobs, size = n, prob = 0.3)
#' param <- c(0.4)
#' 
#' # Derivatives of Binomial log-lik up to order 4
#' llk_binomial(y = y, param = param, n = n, deriv = 4)
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
#'   },
#'   "d4" = function(param){
#'     colSums(do.call("cbind", llk_binomial(y = y, param = param, n = n, deriv = 4)$d4))
#'   })
#'
#' check_deriv(obj = obj, param = param, ord = 1:4)
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