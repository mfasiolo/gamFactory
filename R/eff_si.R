#' Build Single Index Effect Evaluator
#' 
#' @description This function acts as a constructor for single index 
#' smooth effects. It bundles the covariate data and a basis expansion method into an
#' object that can evaluate the effect and its analytical derivatives dynamically.
#' 
#' @param Xi A numeric design matrix of dimensions \code{n} by \code{na} containing 
#' the covariates to be projected down into a single dimension.
#' @param basis A structure or list containing a basis evaluation function 
#' \code{evalX(x, deriv)}. This function must accept a vector of projected values 
#' and return a list containing basis matrices (\code{X0}, \code{X1}, \code{X2}, \code{X3}) 
#' corresponding to the specified derivative order.
#' @param a0 An optional numeric offset vector of length \code{na} added directly 
#' to the single index coefficients \code{alpha}. If \code{NULL}, it defaults 
#' to a vector of zeros matching the column count of \code{Xi}.
#' 
#' @details The function maps covariates using linear transformation. 
#' The matrix \code{Xi} is combined with parameters to construct a univariate index, 
#' which is subsequently transformed via an outer smooth basis expansion.
#' 
#' For a given combined parameter vector \code{param = c(alpha, beta)}:
#' \enumerate{
#'   \item The single index projection step is calculated as:
#'         \deqn{\eta_i = \mathbf{X}_{i \cdot} (\bm{\alpha} + \mathbf{a}_0)}
#'   \item The outer smooth function value is calculated via linear combination with the basis matrix:
#'         \deqn{f(\eta_i) = \mathbf{B}_0(\eta_i) \bm{\beta}}
#' }
#' 
#' When \code{deriv >= 1}, partial derivatives of the spline basis with respect to 
#' the index projection (\eqn{\eta}) are computed analytically via:
#' \deqn{f^{(m)}(\eta_i) = \mathbf{B}_m(\eta_i) \bm{\beta}, \qquad \text{for } m \in \{1, 2, 3\}}
#' 
#' @return An object of class \code{c("si", "nested")} containing:
#' \itemize{
#'   \item{\code{eval(param, deriv = 0)}}{ A function designed to accept a combined 
#'         numeric parameter vector \code{param} (where the first \code{ncol(Xi)} elements 
#'         represent the single index weights \code{alpha}, and the remainder represent 
#'         the spline coefficients \code{beta}) and return a fully evaluated state object.}
#' }
#' @seealso [trans_linear], [s_nest]
#' @name eff_si
#' @rdname eff_si
#' @export eff_si
#'
eff_si <- function(Xi, basis, a0 = NULL){
  
  force(Xi); force(basis); force(a0)
  
  eval <- function(param, deriv = 0){
    
    na <- ncol( Xi )
    nb <- length(param) - na 
    
    # Single index and spline coefficients
    alpha <- param[ 1:na ]
    beta <- param[ -(1:na) ]
    
    if( is.null(a0) ){
      a0 <- alpha * 0
    }
    
    # Project covariates on single index vector 
    ax <- drop( Xi %*% (alpha + a0) )
    
    # Build P-spline basis and its derivatives
    # The error is probably due to the fact that no observations falls within range
    store <- basis$evalX(x = ax, deriv = deriv)
    store$Xi <- Xi
    if( deriv >= 1 ){
      store$f1 <- drop( store$X1 %*% beta )
      store$g1 <- Xi
      if( deriv >= 2 ){
        store$f2 <- drop( store$X2 %*% beta )
        if( deriv >= 3 ){
          store$f3 <- drop( store$X3 %*% beta )
        }
      }
    }
    
    o <- eff_si(Xi = Xi, basis = basis)
    o$f <- drop( store$X0 %*% beta )
    o$param <- param
    o$a0 <- a0
    o$na <- na
    o$store <- store
    o$deriv <- deriv
    return( o )
    
  }
  
  out <- structure(list("eval" = eval), class = c("si", "nested"))
  
  return( out )
  
}









