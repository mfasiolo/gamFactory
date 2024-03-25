#'
#' Log-likelihood of Generalized Pareto Distribution
#' 
#' @description XXX.
#' @param np XXX.
#' @name llk_gpd
#' @rdname llk_gpd
#' @export llk_gpd
#' @examples
#' library(gamFactory)
#' n <- 5
#' param <- c(1.5, 0.3)
#' y <- rexp(n)
#' 
#' # Derivatives of Gaussian log-lik up to order 3
#' llk_gpd(y = y, param = param, deriv = 3)
#' 
#' # Wrap derivatives for compatibility with gamFactory::checkDeriv
#' obj <- list(
#'   "d0" = function(param){
#'     sum(llk_gpd(y = y, param = param, deriv = 0)$d0)
#'   },
#'   "d1" = function(param){
#'     colSums(do.call("cbind", llk_gpd(y = y, param = param, deriv = 1)$d1))
#' 
#'   },
#'   "d2" = function(param){
#'     colSums(do.call("cbind", llk_gpd(y = y, param = param, deriv = 2)$d2))
#' 
#'   },
#'   "d3" = function(param){
#'     colSums(do.call("cbind", llk_gpd(y = y, param = param, deriv = 3)$d3))
#'   })
#' 
#' check_deriv(obj = obj, param = param, ord = 1:3)
#' 
llk_gpd <- function(y, param, deriv = 0, ...) {
  
  if (is.list(param) ) param <- do.call("cbind", param)
  if (is.vector(param)) param <- matrix(param, nrow = 1)
  if (ncol(param) != 2) stop("Wrong number of parameters provided")
  
  phi <- param[ , 1, drop = TRUE]
  xi <- param[ , 2, drop = TRUE]
  
  tol <- 1e-7
  if( any(abs(xi) < tol) ){ 
    xi[xi < tol & xi > 0] <- tol
    xi[xi > -tol & xi <= 0] <- - tol
  }
  
  Cxi <- xi * (1 + xi) * y
  
  out <- list()
  
  out[["d0"]] <- log1p(xi) - log(phi) - (1/xi+1) * log1p(Cxi / phi)
  
  if( deriv > 0 ){
    d1 <- (-phi+y+xi*y)/(phi*(phi+Cxi))
    d2 <- 1/(1+xi)+1/xi^2*log1p(Cxi/phi) - (1/xi+1)*(2*xi+1)*y/(phi+Cxi)
    out[["d1"]] <- list(d1, d2)
    
    if( deriv > 1 ){
      d11 <- (phi^2-2*phi*(1+xi)*y-xi*(1+xi)^2*y^2)/(phi^2*(phi+Cxi)^2)
      d12 <- -(1+xi)*y*(-2*phi+y+xi*y)/(phi*(phi+Cxi)^2)
      d22 <- - 1/(1+xi)^2-2/xi^3*log1p(Cxi/phi) + 
        2/xi^2*(2*xi+1)*y/(phi+Cxi) -
        (1/xi+1)*(2*y*(phi+Cxi)-(2*xi+1)^2*y^2)/(phi+Cxi)^2
      out[["d2"]] <- list(d11, d12, d22) 
      
      if( deriv > 2 ){
        A <- phi+xi*(1+xi)*y
        term_1 <- 2/(1+xi)^3-2*(1+xi)*(y+2*xi*y)^3/(xi*A^3)
        term_2 <- 6*(1+xi)*(1+2*xi)*y^2/(xi*A^2)-3*(y+2*xi*y)^2/(xi^2*A^2)
        term_3 <- 6*y/(xi^2*A)-6*(y+2*xi*y)/(xi^3*A)
        term_4 <- 6*log1p(Cxi/phi)/xi^4
        d111 <- (-2*phi^3+6*phi^2*(1+xi)*y+6*phi*Cxi^2/xi+2*Cxi^3/xi)/(phi^3*(phi+Cxi)^3)
        d112 <-(1+xi)*y*(-4*phi^2+3*phi*(1+xi)*y+xi*(1+xi)^2*y^2)/(phi^2*(phi+Cxi)^3)
        d122 <- 2*y*(phi^2-3*phi*(1+xi)^2*y+(1+xi)^3*y^2)/(phi*(phi+Cxi)^3)
        d222 <- term_1 + term_2 + term_3 + term_4
        out[["d3"]] <- list(d111, d112, d122, d222) 
      }
    }
  }
  
  return( out )
}


