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
    
    na <- ncol( Xi )
    nb <- length(param) - na 
    
    # Inner and outer coefficients
    alpha <- param[ 1:na ]
    beta <- param[ -(1:na) ]
    
    inner <- mgks(y = y, X = X, X0 = Xi, beta = alpha, deriv = deriv)
    
    # Build P-spline basis and its derivatives
    # The error is probably due to the fact that no observations falls within range
    store <- splineDes(x = inner$d0, deriv = deriv)
    store$g <- inner$d0
    store$Xi <- Xi
    if( deriv >= 1 ){
      store$f1 <- drop( store$X1 %*% beta )
      store$g1 <- inner$d1
      if( deriv >= 2 ){
        store$f2 <- drop( store$X2 %*% beta )
        store$g2 <- inner$d2
        if( deriv >= 3 ){
          store$f3 <- drop( store$X3 %*% beta )
          store$g3 <- inner$d3
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

