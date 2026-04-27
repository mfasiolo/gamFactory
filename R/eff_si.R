#'
#' Build single index effect
#' 
#' @param Xi matrix to be projected via single index vector \code{alpha}.
#' @param basis function which takes \code{si = Xi\%*\%alpha} as input and returns model
#'                  matrix and its derivatives w.r.t. \code{si}.
#' @param positive_si logical, whether to constrain the single index to be positive.
#' @name eff_si
#' @rdname eff_si
#' @export eff_si
#'
eff_si <- function(Xi, basis, a0 = NULL, positive_si = FALSE){
  positive_si <- isTRUE(positive_si)
  force(Xi); force(basis); force(a0); force(positive_si)

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
    if (positive_si) {
      ax <- drop( Xi %*% exp(alpha + a0) )
    } else {
      ax <- drop( Xi %*% (alpha + a0) )
    }
    
    # Build P-spline basis and its derivatives
    # The error is probably due to the fact that no observations falls within range
    store <- basis$evalX(x = ax, deriv = deriv)
    store$Xi <- Xi
 
    # derivatives w.r.t ax
    if (deriv >= 1) store$f1 <- drop(store$X1 %*% beta)
    if (deriv >= 2) store$f2 <- drop(store$X2 %*% beta)
    if (deriv >= 3) store$f3 <- drop(store$X3 %*% beta)
    
    # derivatives of g w.r.t alpha. 
      if (!positive_si) {
        # Further derivative g2 and g3 are zero because the transformation is linear
        store$g1 <- Xi
      } else {
        store$g <- ax
        g_diag <- t(t(Xi) * as.vector(exp(alpha + a0)))
        store$g1 <- g_diag
        
        if (deriv >= 2) {
          n_lower <- na * (na + 1) / 2
          store$g2 <- matrix(0, nrow = nrow(Xi), ncol = n_lower)
          idx2 <- if (na == 1) 1 else cumsum(c(1, na:2))
          store$g2[, idx2] <- g_diag
        }
        
        if (deriv >= 3) {
          n_tensor <- na * (na + 1) * (na + 2) / 6
          store$g3 <- matrix(0, nrow = nrow(Xi), ncol = n_tensor)
          idx3 <- if (na == 1) 1 else cumsum(c(1, (na:2 * (na:2 + 1)) / 2))
          store$g3[, idx3] <- g_diag
        }
      }
    
    o <- eff_si(Xi = Xi, basis = basis, a0 = a0, positive_si = positive_si)
    o$f <- drop( store$X0 %*% beta )
    o$param <- param
    o$a0 <- a0
    o$na <- na
    o$store <- store
    o$deriv <- deriv

    return( o )
    
  }
  
  
  if (positive_si) {
    # We need to return the class "si_posi" to trigger the "general" formulas
    # for derivative computation later on (on those for the standard "si" case).
    out <- structure(list("eval" = eval), class = c("si_posi", "nested"))
  } else {
    out <- structure(list("eval" = eval), class = c("si", "nested"))
  }
  return( out )
  
}









