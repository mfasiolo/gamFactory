#'
#' Build single index effect
#' 
#' @param Xi matrix to be projected via single index vector \code{alpha}.
#' @param splineDes function which takes \code{si = Xi\%*\%alpha} as input and returns model
#'                  matrix and its derivatives w.r.t. \code{si}.
#' @name eff_si
#' @rdname eff_si
#' @export eff_si
#'
eff_si <- function(Xi, splineDes){
  
  force(Xi); force(splineDes)

  eval <- function(param, deriv = 0){
    
    na <- ncol( Xi )
    nb <- length(param) - na 
    
    # Single index and spline coefficients
    alpha <- param[ 1:na ]
    beta <- param[ -(1:na) ]
    
    # Project covariates on single index vector 
    ax <- drop( Xi %*% alpha )
    
    # Build P-spline basis and its derivatives
    # The error is probably due to the fact that no observations falls within range
    store <- splineDes(x = ax, deriv = deriv)
    store$Xi <- Xi
    if( deriv >= 1 ){
      store$f1 <- drop( store$X1 %*% beta )
      if( deriv >= 2 ){
        store$f2 <- drop( store$X2 %*% beta )
        if( deriv >= 3 ){
          store$f3 <- drop( store$X3 %*% beta )
        }
      }
    }
    
    o <- eff_si(Xi = Xi, splineDes = splineDes)
    o$f <- drop( store$X0 %*% beta )
    o$param <- param
    o$na <- na
    o$store <- store
    o$deriv <- deriv

    return( o )
    
  }
  
  out <- structure(list("eval" = eval), class = c("singleIndex"))
  
  return( out )
  
}









