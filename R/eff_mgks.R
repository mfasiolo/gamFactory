#'
#' Build nested adaptive exponential smooth effect
#' 
#' @name eff_mgks
#' @rdname eff_mgks
#' @export eff_mgks
#'
eff_mgks <- function(y, X, Xi, splineDes){
  
  force(y); force(X); force(Xi); force(splineDes);
  
  eval <- function(param, deriv = 0){
    
    na <- ncol( Xi ) + 1
    nb <- length(param) - na 
    
    # Scale, inner and outer coefficients
    a0 <- exp( param[1] )
    alpha <- param[ 2:na ]
    beta <- param[ -(1:na) ]
    
    inner <- mgks(y = y, X = X, X0 = Xi, beta = alpha, deriv = deriv)
    
    # Build P-spline basis and its derivatives
    # We are also getting the derivatives of the inner linear predictor w.r.t.
    # the scale parameter a0 and the coefficients alpha
    store <- splineDes(x = a0 * inner$d0, deriv = deriv)
    store$g <- a0 * inner$d0
    store$Xi <- Xi
    if( deriv >= 1 ){
      store$f1 <- drop( store$X1 %*% beta )
      store$g1 <-  cbind(store$g, a0 * inner$d1)
      if( deriv >= 2 ){
        store$f2 <- drop( store$X2 %*% beta )
        store$g2 <- cbind(store$g1, a0 * inner$d2)
        if( deriv >= 3 ){
          store$f3 <- drop( store$X3 %*% beta )
          store$g3 <- cbind(store$g2, a0 * inner$d3)
        }
      }
    }
    
    o <- eff_mgks(y = y, X = X, Xi = Xi, splineDes = splineDes)
    o$f <- drop( store$X0 %*% beta )
    o$param <- param
    o$na <- na
    o$store <- store
    o$deriv <- deriv
    
    return( o )
    
  }
  
  out <- structure(list("eval" = eval), class = c("mgks", "nested"))
  
  return( out )
  
}

